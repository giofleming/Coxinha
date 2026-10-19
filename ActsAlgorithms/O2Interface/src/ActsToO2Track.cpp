#include "ActsToO2Track.hpp"

#include "Acts/Definitions/Units.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace
{

constexpr float kPi = 3.14159265358979323846f;
constexpr float kTwoPi = 2.f * kPi;
constexpr float kHalfPi = 0.5f * kPi;
constexpr int kNSectors = 18;
constexpr float kSectorSpan = kTwoPi / static_cast<float>(kNSectors);
constexpr float kSafe = 1.e-5f;

float bringToPlusMinusPi(float angle) noexcept
{
  while (angle >= kPi) {
    angle -= kTwoPi;
  }
  while (angle < -kPi) {
    angle += kTwoPi;
  }
  return angle;
}

float angleToSectorAlpha(const float angle) noexcept
{
  float normalized = angle;
  while (normalized < 0.f) {
    normalized += kTwoPi;
  }
  while (normalized >= kTwoPi) {
    normalized -= kTwoPi;
  }
  const int sector = static_cast<int>(normalized / kSectorSpan) % kNSectors;
  return bringToPlusMinusPi((static_cast<float>(sector) + 0.5f) * kSectorSpan);
}

void rotateZ(standalone_o2::TrackParametrization::dim3_t& vector,
             const float angle) noexcept
{
  const float sine = std::sin(angle);
  const float cosine = std::cos(angle);
  const float x = vector[0];
  const float y = vector[1];
  vector[0] = x * cosine - y * sine;
  vector[1] = x * sine + y * cosine;
}

std::size_t covarianceIndex(std::size_t row, std::size_t column)
{
  if (row >= standalone_o2::kNParams ||
      column >= standalone_o2::kNParams) {
    throw std::out_of_range("O2 covariance index is outside [0,4]");
  }
  if (column > row) {
    std::swap(row, column);
  }
  return row * (row + 1) / 2 + column;
}

void requireFinite(const Acts::Vector3& vector, const char* name)
{
  if (!std::isfinite(vector[0]) || !std::isfinite(vector[1]) ||
      !std::isfinite(vector[2])) {
    throw std::invalid_argument(
      std::string{"Cannot convert ACTS track: non-finite "} + name);
  }
}

int convertCharge(const double actsCharge)
{
  const double chargeInE = actsCharge / Acts::UnitConstants::e;
  if (!std::isfinite(chargeInE)) {
    throw std::invalid_argument("Cannot convert ACTS track: non-finite charge");
  }

  const double roundedCharge = std::round(chargeInE);
  if (std::abs(chargeInE - roundedCharge) > 1.e-6) {
    throw std::invalid_argument(
      "Cannot convert ACTS track: standalone O2 track requires an integer charge");
  }
  if (roundedCharge < static_cast<double>(std::numeric_limits<int>::min()) ||
      roundedCharge > static_cast<double>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument("Cannot convert ACTS track: charge is outside the int range");
  }
  return static_cast<int>(roundedCharge);
}

struct ConvertedKinematics {
  standalone_o2::TrackParametrization::dim3_t position;
  standalone_o2::TrackParametrization::dim3_t momentum;
  int charge;
};

ConvertedKinematics convertKinematics(
  const Acts::BoundTrackParameters& actsTrack,
  const Acts::GeometryContext& geometryContext)
{
  const Acts::Vector3 position = actsTrack.position(geometryContext);
  const Acts::Vector3 momentum = actsTrack.momentum();
  requireFinite(position, "position");
  requireFinite(momentum, "momentum");

  const double px = momentum[0] / Acts::UnitConstants::GeV;
  const double py = momentum[1] / Acts::UnitConstants::GeV;
  if (std::hypot(px, py) <= 0.) {
    throw std::invalid_argument(
      "Cannot convert ACTS track: standalone O2 track requires non-zero transverse momentum");
  }

  return {
    {static_cast<float>(position[0] / Acts::UnitConstants::cm),
     static_cast<float>(position[1] / Acts::UnitConstants::cm),
     static_cast<float>(position[2] / Acts::UnitConstants::cm)},
    {static_cast<float>(px), static_cast<float>(py),
     static_cast<float>(momentum[2] / Acts::UnitConstants::GeV)},
    convertCharge(actsTrack.charge())};
}

template <typename T, typename = void>
struct HasFixedMomentumHypothesis : std::false_type {
};

template <typename T>
struct HasFixedMomentumHypothesis<
  T, std::void_t<decltype(std::declval<const T&>().hasMomentumHypothesis())>>
  : std::true_type {
};

double momentumDerivative(const Acts::BoundTrackParameters& actsTrack,
                          const double momentum,
                          const double qOverP)
{
  const auto& hypothesis = actsTrack.particleHypothesis();
  if constexpr (HasFixedMomentumHypothesis<
                  std::decay_t<decltype(hypothesis)>>::value) {
    if (hypothesis.hasMomentumHypothesis()) {
      return 0.;
    }
  }
  return -momentum / qOverP;
}

standalone_o2::TrackParametrizationWithError::labCovMat_t
  convertCovariance(const Acts::BoundTrackParameters& actsTrack,
                    const Acts::GeometryContext& geometryContext)
{
  const auto& optionalCovariance = actsTrack.covariance();
  if (!optionalCovariance.has_value()) {
    throw std::invalid_argument(
      "Cannot convert ACTS track covariance: input covariance is absent");
  }

  const Acts::BoundMatrix& boundCovariance = *optionalCovariance;
  for (unsigned int row = 0; row < Acts::eBoundSize; ++row) {
    for (unsigned int column = 0; column < Acts::eBoundSize; ++column) {
      if (!std::isfinite(boundCovariance(row, column))) {
        throw std::invalid_argument(
          "Cannot convert ACTS track covariance: non-finite matrix element");
      }
    }
  }

  const Acts::Vector3 position = actsTrack.position(geometryContext);
  const Acts::Vector3 direction = actsTrack.direction();
  const double momentum = actsTrack.absoluteMomentum();
  const double qOverP = actsTrack.qOverP();
  requireFinite(position, "position");
  requireFinite(direction, "direction");
  if (!std::isfinite(momentum) || momentum <= 0. ||
      !std::isfinite(qOverP) || qOverP == 0.) {
    throw std::invalid_argument(
      "Cannot convert ACTS track covariance: invalid momentum or q/p");
  }

  const Acts::BoundToFreeMatrix boundToFree =
    actsTrack.referenceSurface().boundToFreeJacobian(
      geometryContext, position, direction);
  const double dpDqOverP =
    momentumDerivative(actsTrack, momentum, qOverP);

  std::array<std::array<double, Acts::eBoundSize>, 6> jacobian{};
  for (unsigned int bound = 0; bound < Acts::eBoundSize; ++bound) {
    for (unsigned int component = 0; component < 3; ++component) {
      jacobian[component][bound] =
        boundToFree(Acts::eFreePos0 + component, bound) /
        Acts::UnitConstants::cm;
      jacobian[component + 3][bound] =
        (momentum * boundToFree(Acts::eFreeDir0 + component, bound) +
         direction[component] * dpDqOverP *
           boundToFree(Acts::eFreeQOverP, bound)) /
        Acts::UnitConstants::GeV;
    }
  }

  double globalCovariance[6][6]{};
  for (unsigned int row = 0; row < 6; ++row) {
    for (unsigned int column = 0; column < 6; ++column) {
      for (unsigned int i = 0; i < Acts::eBoundSize; ++i) {
        for (unsigned int j = 0; j < Acts::eBoundSize; ++j) {
          globalCovariance[row][column] +=
            jacobian[row][i] * boundCovariance(i, j) *
            jacobian[column][j];
        }
      }
    }
  }

  standalone_o2::TrackParametrizationWithError::labCovMat_t packed{};
  std::size_t index = 0;
  for (unsigned int row = 0; row < 6; ++row) {
    for (unsigned int column = 0; column <= row; ++column) {
      const double value = globalCovariance[row][column];
      if (!std::isfinite(value)) {
        throw std::invalid_argument(
          "Cannot convert ACTS track covariance: non-finite transformed element");
      }
      packed[index++] = static_cast<float>(value);
    }
  }
  return packed;
}

} // namespace

namespace standalone_o2
{

TrackParametrization::TrackParametrization(
  const dim3_t& xyz, const dim3_t& pxpypz, const int charge,
  const bool sectorAlpha, const PID pid)
  : mPID{pid}
{
  const float radialPosition2 = xyz[0] * xyz[0] + xyz[1] * xyz[1];
  float alpha = (sectorAlpha || radialPosition2 < 1.f)
                  ? std::atan2(pxpypz[1], pxpypz[0])
                  : std::atan2(xyz[1], xyz[0]);
  if (sectorAlpha) {
    alpha = angleToSectorAlpha(alpha);
  }

  float sine = std::sin(alpha);
  float cosine = std::cos(alpha);
  if (cosine * pxpypz[0] + sine * pxpypz[1] < 0.f) {
    alpha = std::atan2(pxpypz[1], pxpypz[0]);
    if (sectorAlpha) {
      alpha = angleToSectorAlpha(alpha);
    }
    sine = std::sin(alpha);
    cosine = std::cos(alpha);
  }

  if (std::abs(sine) < 2.f * kSafe) {
    if (alpha > 0.f) {
      alpha += alpha < kHalfPi ? 2.f * kSafe : -2.f * kSafe;
    } else {
      alpha += alpha > -kHalfPi ? -2.f * kSafe : 2.f * kSafe;
    }
  } else if (std::abs(cosine) < 2.f * kSafe) {
    if (alpha > 0.f) {
      alpha += alpha > kHalfPi ? 2.f * kSafe : -2.f * kSafe;
    } else {
      alpha += alpha > -kHalfPi ? 2.f * kSafe : -2.f * kSafe;
    }
  }

  dim3_t localPosition = xyz;
  dim3_t localMomentum = pxpypz;
  rotateZ(localPosition, -alpha);
  rotateZ(localMomentum, -alpha);

  const float inversePt =
    1.f / std::hypot(localMomentum[0], localMomentum[1]);
  mX = localPosition[0];
  mAlpha = alpha;
  mParameters[kY] = localPosition[1];
  mParameters[kZ] = localPosition[2];
  mParameters[kSnp] = localMomentum[1] * inversePt;
  mParameters[kTgl] = localMomentum[2] * inversePt;
  mAbsCharge = std::abs(charge);
  mParameters[kQ2Pt] = charge != 0 ? charge * inversePt : inversePt;

  if (std::abs(1.f - mParameters[kSnp]) < kSafe) {
    mParameters[kSnp] = 1.f - kSafe;
  } else if (std::abs(-1.f - mParameters[kSnp]) < kSafe) {
    mParameters[kSnp] = -1.f + kSafe;
  }
}

int TrackParametrization::getCharge() const noexcept
{
  if (mAbsCharge == 0) {
    return 0;
  }
  return getQ2Pt() < 0.f ? -mAbsCharge : mAbsCharge;
}

TrackParametrization::value_t TrackParametrization::getPt() const noexcept
{
  const float effectiveCharge =
    mAbsCharge == 0 ? 1.f : static_cast<float>(mAbsCharge);
  return effectiveCharge / std::abs(getQ2Pt());
}

TrackParametrization::value_t TrackParametrization::getP() const noexcept
{
  return getPt() * std::sqrt(1.f + getTgl() * getTgl());
}

void TrackParametrization::getXYZGlo(dim3_t& xyz) const noexcept
{
  xyz = {mX, getY(), getZ()};
  rotateZ(xyz, mAlpha);
}

bool TrackParametrization::getPxPyPzGlo(dim3_t& pxpypz) const noexcept
{
  if (!std::isfinite(getQ2Pt()) || std::abs(getQ2Pt()) <= 0.f ||
      std::abs(getSnp()) >= 1.f) {
    return false;
  }
  const float pt = getPt();
  const float localPx =
    pt * std::sqrt((1.f - getSnp()) * (1.f + getSnp()));
  pxpypz = {localPx, pt * getSnp(), pt * getTgl()};
  rotateZ(pxpypz, mAlpha);
  return true;
}

TrackParametrizationWithError::TrackParametrizationWithError(
  const dim3_t& xyz, const dim3_t& pxpypz,
  const labCovMat_t& globalCovariance, const int charge,
  const bool sectorAlpha, const PID pid)
  : TrackParametrization{xyz, pxpypz, charge, sectorAlpha, pid}
{
  float fullGlobalCovariance[6][6]{};
  std::size_t packedIndex = 0;
  for (std::size_t row = 0; row < 6; ++row) {
    for (std::size_t column = 0; column <= row; ++column) {
      fullGlobalCovariance[row][column] =
        fullGlobalCovariance[column][row] =
          globalCovariance[packedIndex++];
    }
  }

  dim3_t localMomentum = pxpypz;
  rotateZ(localMomentum, -getAlpha());
  const float sine = std::sin(getAlpha());
  const float cosine = std::cos(getAlpha());
  const float pt2 = localMomentum[0] * localMomentum[0] +
                    localMomentum[1] * localMomentum[1];
  const float inversePt = 1.f / std::sqrt(pt2);
  const float inversePt3 = inversePt / pt2;
  const float effectiveCharge =
    charge != 0 ? static_cast<float>(charge) : 1.f;

  float jacobian[kNParams][6]{};
  jacobian[kY][0] = -sine;
  jacobian[kY][1] = cosine;
  jacobian[kZ][2] = 1.f;

  const float u = localMomentum[0];
  const float v = localMomentum[1];
  const float w = localMomentum[2];
  const float dSnpDu = -u * v * inversePt3;
  const float dSnpDv = u * u * inversePt3;
  const float dTglDu = -w * u * inversePt3;
  const float dTglDv = -w * v * inversePt3;
  const float dQ2PtDu = -effectiveCharge * u * inversePt3;
  const float dQ2PtDv = -effectiveCharge * v * inversePt3;

  jacobian[kSnp][3] = dSnpDu * cosine - dSnpDv * sine;
  jacobian[kSnp][4] = dSnpDu * sine + dSnpDv * cosine;
  jacobian[kTgl][3] = dTglDu * cosine - dTglDv * sine;
  jacobian[kTgl][4] = dTglDu * sine + dTglDv * cosine;
  jacobian[kTgl][5] = inversePt;
  jacobian[kQ2Pt][3] = dQ2PtDu * cosine - dQ2PtDv * sine;
  jacobian[kQ2Pt][4] = dQ2PtDu * sine + dQ2PtDv * cosine;

  for (std::size_t row = 0; row < kNParams; ++row) {
    for (std::size_t column = 0; column <= row; ++column) {
      float value = 0.f;
      for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
          value += jacobian[row][i] * fullGlobalCovariance[i][j] *
                   jacobian[column][j];
        }
      }
      mCovariance[covarianceIndex(row, column)] = value;
    }
  }
}

TrackParametrizationWithError::value_t
  TrackParametrizationWithError::getCovariance(
    const std::size_t row, const std::size_t column) const
{
  return mCovariance[covarianceIndex(row, column)];
}

TrackParametrization convertActsToO2Track(
  const Acts::BoundTrackParameters& actsTrack,
  const Acts::GeometryContext& geometryContext,
  const bool sectorAlpha, const PID pid)
{
  const ConvertedKinematics converted =
    convertKinematics(actsTrack, geometryContext);
  return TrackParametrization{converted.position, converted.momentum,
                              converted.charge, sectorAlpha, pid};
}

TrackParametrizationWithError convertActsToO2TrackWithCovariance(
  const Acts::BoundTrackParameters& actsTrack,
  const Acts::GeometryContext& geometryContext,
  const bool sectorAlpha, const PID pid)
{
  const ConvertedKinematics converted =
    convertKinematics(actsTrack, geometryContext);
  const TrackParametrizationWithError::labCovMat_t covariance =
    convertCovariance(actsTrack, geometryContext);
  return TrackParametrizationWithError{
    converted.position, converted.momentum, covariance,
    converted.charge, sectorAlpha, pid};
}

} // namespace standalone_o2
