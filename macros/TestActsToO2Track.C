#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "O2Interface/ActsToO2Track.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"

#include <TCanvas.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TPad.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <random>
#include <sstream>
#include <string>

namespace
{

constexpr int kActsColor = kRed + 1;
constexpr int kO2Color = kBlue + 1;
constexpr int kNumberOfDrawPoints = 100;

template <typename Position, typename Momentum>
TGraph* makeProjection(const Position& position, const Momentum& momentum, const bool drawRZ, const double lengthCm)
{
  const double norm = std::sqrt(momentum[0] * momentum[0] + momentum[1] * momentum[1] + momentum[2] * momentum[2]);
  auto* graph = new TGraph(kNumberOfDrawPoints);
  graph->SetBit(TObject::kCanDelete);
  for (int point = 0; point < kNumberOfDrawPoints; ++point) {
    const double path = lengthCm * point / (kNumberOfDrawPoints - 1.);
    const double x = position[0] + path * momentum[0] / norm;
    const double y = position[1] + path * momentum[1] / norm;
    const double z = position[2] + path * momentum[2] / norm;
    graph->SetPoint(point, drawRZ ? std::hypot(x, y) : x, drawRZ ? z : y);
  }
  return graph;
}

void styleTrack(TGraph& graph, const int color, const int lineStyle, const int lineWidth)
{
  graph.SetLineColor(color);
  graph.SetLineStyle(lineStyle);
  graph.SetLineWidth(lineWidth);
}

void drawStartMarker(const double x, const double y, const int color, const int style)
{
  auto* marker = new TMarker(x, y, style);
  marker->SetMarkerColor(color);
  marker->SetMarkerSize(1.1);
  marker->SetBit(TObject::kCanDelete);
  marker->Draw();
}

template <typename ActsPosition, typename ActsMomentum,
          typename O2Position, typename O2Momentum>
void drawTrackComparison(TCanvas& canvas,
                         const ActsPosition& actsPosition,
                         const ActsMomentum& actsMomentum,
                         const O2Position& o2Position,
                         const O2Momentum& o2Momentum,
                         const std::size_t trackNumber,
                         const double lengthCm)
{
  canvas.Clear();
  canvas.Divide(2, 1);

  std::ostringstream title;
  title << "Track " << trackNumber << " (straight tangent, L=" << lengthCm << " cm)";

  canvas.cd(1);
  gPad->SetGrid();
  auto* actsXY = makeProjection(actsPosition, actsMomentum, false, lengthCm);
  auto* o2XY = makeProjection(o2Position, o2Momentum, false, lengthCm);
  styleTrack(*actsXY, kActsColor, 1, 4);
  styleTrack(*o2XY, kO2Color, 7, 3);
  actsXY->SetTitle((title.str() + ";x (cm);y (cm)").c_str());
  actsXY->Draw("AL");
  o2XY->Draw("L SAME");
  drawStartMarker(actsPosition[0], actsPosition[1], kActsColor, 20);
  drawStartMarker(o2Position[0], o2Position[1], kO2Color, 24);
  auto* legend = new TLegend(0.12, 0.76, 0.45, 0.90);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(actsXY, "ACTS", "l");
  legend->AddEntry(o2XY, "converted O2", "l");
  legend->SetBit(TObject::kCanDelete);
  legend->Draw();

  canvas.cd(2);
  gPad->SetGrid();
  auto* actsRZ = makeProjection(actsPosition, actsMomentum, true, lengthCm);
  auto* o2RZ = makeProjection(o2Position, o2Momentum, true, lengthCm);
  styleTrack(*actsRZ, kActsColor, 1, 4);
  styleTrack(*o2RZ, kO2Color, 7, 3);
  actsRZ->SetTitle((title.str() + ";r (cm);z (cm)").c_str());
  actsRZ->Draw("AL");
  o2RZ->Draw("L SAME");
  drawStartMarker(std::hypot(actsPosition[0], actsPosition[1]),
                  actsPosition[2], kActsColor, 20);
  drawStartMarker(std::hypot(o2Position[0], o2Position[1]),
                  o2Position[2], kO2Color, 24);

  auto* note = new TLatex(
    0.12, 0.94, "Lines should overlap; O2 is dashed for visibility");
  note->SetNDC();
  note->SetTextSize(0.035);
  note->SetBit(TObject::kCanDelete);
  note->Draw();

  canvas.Modified();
  canvas.Update();
}

} // namespace

bool TestActsToO2Track(const Acts::BoundTrackParameters& actsTrack,
                       const Acts::GeometryContext& geometryContext,
                       const double relativeTolerance = 1.e-5,
                       const double eigenvalueTolerance = 1.e-8,
                       TCanvas* comparisonCanvas = nullptr,
                       const std::size_t trackNumber = 0,
                       const double drawingLengthCm = 50.)
{
  try {
    const auto converted = standalone_o2::convertActsToO2TrackWithCovariance(actsTrack, geometryContext);
    standalone_o2::TrackParametrization::dim3_t convertedPosition{};
    standalone_o2::TrackParametrization::dim3_t convertedMomentum{};
    converted.getXYZGlo(convertedPosition);
    const bool validMomentum = converted.getPxPyPzGlo(convertedMomentum);
    const Acts::Vector3 actsPositionNative = actsTrack.position(geometryContext);
    const Acts::Vector3 actsMomentumNative = actsTrack.momentum();
    const std::array<double, 3> actsPosition{actsPositionNative[0] / Acts::UnitConstants::cm, actsPositionNative[1] / Acts::UnitConstants::cm, actsPositionNative[2] / Acts::UnitConstants::cm};
    const std::array<double, 3> actsMomentum{actsMomentumNative[0] / Acts::UnitConstants::GeV, actsMomentumNative[1] / Acts::UnitConstants::GeV, actsMomentumNative[2] / Acts::UnitConstants::GeV};
    const int actsCharge = static_cast<int>(std::lround(actsTrack.charge() / Acts::UnitConstants::e));
    const auto close = [relativeTolerance](const double lhs, const double rhs) { const double scale = std::max({1., std::abs(lhs), std::abs(rhs)}); return std::abs(lhs - rhs) <= relativeTolerance * scale; };
    bool kinematicsOK = validMomentum && converted.getCharge() == actsCharge;
    std::cout << std::setprecision(9) << "component       ACTS             converted        delta\n";
    constexpr std::array<const char*, 6> labels{"x [cm]", "y [cm]", "z [cm]", "px [GeV]", "py [GeV]", "pz [GeV]"};
    for (int component = 0; component < 3; ++component) {
      const double reference = actsPosition[component];
      const double value = convertedPosition[component];
      kinematicsOK = kinematicsOK && close(reference, value);
      std::cout << std::setw(11) << labels[component] << "  " << std::setw(16) << reference << "  " << std::setw(16) << value << "  " << std::setw(12) << value - reference << '\n';
    }
    for (int component = 0; component < 3; ++component) {
      const double reference = actsMomentum[component];
      const double value = convertedMomentum[component];
      kinematicsOK = kinematicsOK && close(reference, value);
      std::cout << std::setw(11) << labels[component + 3] << "  " << std::setw(16) << reference << "  " << std::setw(16) << value << "  " << std::setw(12) << value - reference << '\n';
    }
    Eigen::Matrix<double, standalone_o2::kNParams, standalone_o2::kNParams> covariance;
    bool covarianceFinite = true;
    for (std::size_t row = 0; row < standalone_o2::kNParams; ++row) {
      for (std::size_t column = 0; column < standalone_o2::kNParams; ++column) {
        covariance(row, column) = converted.getCovariance(row, column);
        covarianceFinite = covarianceFinite && std::isfinite(covariance(row, column));
      }
    }
    const double asymmetry = (covariance - covariance.transpose()).cwiseAbs().maxCoeff();
    const Eigen::SelfAdjointEigenSolver<decltype(covariance)> solver{covariance};
    const bool eigensolverOK = solver.info() == Eigen::Success;
    const double minimumEigenvalue = eigensolverOK ? solver.eigenvalues().minCoeff() : -std::numeric_limits<double>::infinity();
    const double covarianceScale = std::max(1., covariance.cwiseAbs().maxCoeff());
    const bool covarianceOK = covarianceFinite && asymmetry <= relativeTolerance * covarianceScale && eigensolverOK && minimumEigenvalue >= -eigenvalueTolerance * covarianceScale;
    std::cout << "charge: ACTS=" << actsCharge << ", converted=" << converted.getCharge() << '\n'
              << "covariance max asymmetry: " << asymmetry << '\n'
              << "covariance minimum eigenvalue: " << minimumEigenvalue << '\n';
    if (comparisonCanvas != nullptr && validMomentum) {
      drawTrackComparison(*comparisonCanvas, actsPosition, actsMomentum, convertedPosition, convertedMomentum, trackNumber, drawingLengthCm);
    }
    const bool passed = kinematicsOK && covarianceOK;
    std::cout << (passed ? "SANITY CHECK PASSED" : "SANITY CHECK FAILED") << '\n';
    return passed;
  } catch (const std::exception& error) {
    std::cerr << "SANITY CHECK FAILED: " << error.what() << '\n';
    return false;
  }
}
/// drawTracks=false keeps the original non-graphical behavior.
/// If outputPdf is non-empty, every drawn track is saved as one PDF page.
bool TestActsToO2Track(const std::size_t nTracks = 100,
                       const bool drawTracks = true,
                       const char* outputPdf = "ActsToO2TrackComparison.pdf",
                       const double drawingLengthCm = 50.)
{
  std::mt19937 random{0xAC7502};
  std::uniform_real_distribution<double> localPosition{-100., 100.};
  std::uniform_real_distribution<double> azimuth{-std::numbers::pi, std::numbers::pi};
  std::uniform_real_distribution<double> polarAngle{0.1 * std::numbers::pi, 0.9 * std::numbers::pi};
  std::uniform_real_distribution<double> momentum{0.1, 100.};
  std::bernoulli_distribution positiveCharge{0.5};
  const auto geometryContext = Acts::GeometryContext::dangerouslyDefaultConstruct();
  const auto surface = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3::Zero());
  static std::unique_ptr<TCanvas> comparisonCanvas;
  const bool savePdf = drawTracks && outputPdf != nullptr && outputPdf[0] != '\0';
  if (drawTracks) {
    comparisonCanvas = std::make_unique<TCanvas>("actsToO2TrackCanvas", "ACTS to O2 track conversion", 1400, 650);
    if (savePdf) {
      comparisonCanvas->Print((std::string{outputPdf} + "[").c_str());
    }
  }
  bool allPassed = true;
  for (std::size_t trackNumber = 0; trackNumber < nTracks; ++trackNumber) {
    Acts::BoundVector parameters = Acts::BoundVector::Zero();
    parameters[Acts::eBoundLoc0] = localPosition(random) * Acts::UnitConstants::cm;
    parameters[Acts::eBoundLoc1] = localPosition(random) * Acts::UnitConstants::cm;
    parameters[Acts::eBoundPhi] = azimuth(random);
    parameters[Acts::eBoundTheta] = polarAngle(random);
    const double charge = positiveCharge(random) ? 1. : -1.;
    parameters[Acts::eBoundQOverP] = charge * Acts::UnitConstants::e / (momentum(random) * Acts::UnitConstants::GeV);
    Acts::BoundMatrix covariance = Acts::BoundMatrix::Zero();
    covariance.diagonal() << 1.e-4, 1.e-4, 1.e-6, 1.e-6, 1.e-8, 1.e-4;
    const Acts::BoundTrackParameters track{surface, parameters, covariance, Acts::ParticleHypothesis::pion()};
    std::cout << "Random track " << trackNumber + 1 << '/' << nTracks << '\n';
    const bool passed = TestActsToO2Track(track, geometryContext, 1.e-5, 1.e-8, drawTracks ? comparisonCanvas.get() : nullptr, trackNumber + 1, drawingLengthCm);
    if (drawTracks && savePdf) {
      comparisonCanvas->Print(outputPdf);
    }
    if (!passed) {
      allPassed = false;
      break;
    }
  }
  if (drawTracks && savePdf) {
    comparisonCanvas->Print((std::string{outputPdf} + "]").c_str());
    std::cout << "Track comparison pages written to " << outputPdf << '\n';
  }
  if (allPassed) {
    std::cout << "All " << nTracks << " random tracks passed the sanity check\n";
  }
  return allPassed;
}