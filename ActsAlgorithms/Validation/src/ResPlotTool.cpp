// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ResPlotTool.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cmath>
#include <format>
#include <iomanip>
#include <sstream>

using namespace ActsExamples;

namespace AliceActsTrk
{

ResPlotTool::ResPlotTool(const ResPlotTool::Config& cfg,
                         Acts::Logging::Level lvl)
  : m_cfg(cfg), m_logger(Acts::getDefaultLogger("ResPlotTool", lvl))
{
  const auto& etaAxis = m_cfg.varBinning.at("Eta");
  const auto& phiAxis = m_cfg.varBinning.at("Phi");
  const auto& ptAxis = m_cfg.varBinning.at("Pt");
  const auto& pullAxis = m_cfg.varBinning.at("Pull");

  ACTS_DEBUG("Initialize the histograms for residual and pull plots");
  for (unsigned int parID = 0; parID < Acts::eBoundSize; parID++) {
    std::string parName = m_cfg.paramNames.at(parID);

    std::string parResidual = "Residual_" + parName;
    const auto& residualAxis = m_cfg.varBinning.at(parResidual);

    // residual distributions
    m_res.emplace(parName, Histogram1(
                             std::format("res_{}", parName),
                             std::format("Residual of {}", parName),
                             std::array{residualAxis}));

    // residual vs eta scatter plots
    m_resVsEta.emplace(parName,
                       Histogram2(
                         std::format("res_{}_vs_eta", parName),
                         std::format("Residual of {} vs eta", parName),
                         std::array{etaAxis, residualAxis}));

    // residual vs phi scatter plots
    m_resVsPhi.emplace(parName,
                       Histogram2(
                         std::format("res_{}_vs_phi", parName),
                         std::format("Residual of {} vs phi", parName),
                         std::array{phiAxis, residualAxis}));

    // residual vs pT scatter plots
    m_resVsPt.emplace(parName, Histogram2(
                                 std::format("res_{}_vs_pT", parName),
                                 std::format("Residual of {} vs pT", parName),
                                 std::array{ptAxis, residualAxis}));

    // pull distributions
    m_pull.emplace(
      parName, Histogram1(
                 std::format("pull_{}", parName),
                 std::format("Pull of {}", parName), std::array{pullAxis}));

    // pull vs eta scatter plots
    m_pullVsEta.emplace(parName, Histogram2(
                                   std::format("pull_{}_vs_eta", parName),
                                   std::format("Pull of {} vs eta", parName),
                                   std::array{etaAxis, pullAxis}));

    // pull vs eta scatter plots
    m_pullVsPhi.emplace(parName, Histogram2(
                                   std::format("pull_{}_vs_phi", parName),
                                   std::format("Pull of {} vs phi", parName),
                                   std::array{phiAxis, pullAxis}));

    // pull vs pT scatter plots
    m_pullVsPt.emplace(parName, Histogram2(
                                  std::format("pull_{}_vs_pT", parName),
                                  std::format("Pull of {} vs pT", parName),
                                  std::array{ptAxis, pullAxis}));
  }

  // Add momentum resolution
  // residual distributions
  m_res.emplace("res_pt_o_pt", Histogram1("res_pt_o_pt",
                                          "Residual of pT over pT",
                                          std::array{m_cfg.varBinning.at("Residual_pt_o_pt")}));

  m_resVsEta.emplace("res_pt_o_pt_vs_eta", Histogram2("res_pt_o_pt_vs_eta",
                                                      "Residual of pT over pT",
                                                      std::array{etaAxis, m_cfg.varBinning.at("Residual_pt_o_pt")}));

  m_resVsPhi.emplace("res_pt_o_pt_vs_phi", Histogram2("res_pt_o_pt_vs_phi",
                                                      "Residual of pT over pT",
                                                      std::array{etaAxis, m_cfg.varBinning.at("Residual_pt_o_pt")}));

  m_resVsEtaVsPhi.emplace("res_pt_o_pt_vs_eta_vs_phi", Histogram3("res_pt_o_pt_vs_eta_vs_phi",
                                                                  "Residual of pT over pT",
                                                                  std::array{etaAxis, phiAxis, m_cfg.varBinning.at("Residual_pt_o_pt")}));

  for (unsigned int etaBin = 0; etaBin < m_cfg.etaBins.size() - 1; etaBin++) {

    const auto& residualAxis = m_cfg.varBinning.at("Residual_pt_o_pt");
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << m_cfg.etaBins[etaBin];
    std::string etaBin_low = ss.str();
    ss.str("");

    ss << std::fixed << std::setprecision(2) << m_cfg.etaBins[etaBin + 1];
    std::string etaBin_high = ss.str();
    ss.str("");

    std::string resName = "res_ptopt_" + etaBin_low + "_" + etaBin_high;

    m_resVsPt.emplace(resName, Histogram2(resName,
                                          resName,
                                          std::array{ptAxis, residualAxis}));
  }
}

void ResPlotTool::fill(const Acts::GeometryContext& gctx,
                       const SimParticleState& truthParticle,
                       const Acts::BoundTrackParameters& fittedParamters)
{
  using ParametersVector = Acts::BoundTrackParameters::ParametersVector;
  using Acts::VectorHelpers::eta;
  using Acts::VectorHelpers::perp;
  using Acts::VectorHelpers::phi;
  using Acts::VectorHelpers::theta;

  // get the fitted parameter (at perigee surface) and its error
  Acts::BoundVector trackParameter = fittedParamters.parameters();

  // get the perigee surface
  const Acts::Surface& pSurface = fittedParamters.referenceSurface();

  // get the truth position and momentum
  ParametersVector truthParameter = ParametersVector::Zero();

  // get the truth perigee parameter
  Acts::Intersection3D intersection =
    pSurface
      .intersect(gctx, truthParticle.position(), truthParticle.direction())
      .closest();
  if (intersection.isValid()) {
    auto lpResult = pSurface.globalToLocal(gctx, intersection.position(),
                                           truthParticle.direction());
    assert(lpResult.ok());

    truthParameter[Acts::BoundIndices::eBoundLoc0] =
      lpResult.value()[Acts::BoundIndices::eBoundLoc0];
    truthParameter[Acts::BoundIndices::eBoundLoc1] =
      lpResult.value()[Acts::BoundIndices::eBoundLoc1];
  } else {
    ACTS_ERROR("Cannot get the truth perigee parameter");
  }
  truthParameter[Acts::BoundIndices::eBoundPhi] =
    phi(truthParticle.direction());
  truthParameter[Acts::BoundIndices::eBoundTheta] =
    theta(truthParticle.direction());
  truthParameter[Acts::BoundIndices::eBoundQOverP] = truthParticle.qOverP();
  truthParameter[Acts::BoundIndices::eBoundTime] = truthParticle.time();

  // get the truth eta and pT
  const auto truthEta = eta(truthParticle.direction());
  const auto truthPt = truthParticle.transverseMomentum();
  const auto truthPhi = truthParameter[Acts::BoundIndices::eBoundPhi];

  // fill the histograms for residual and pull
  for (unsigned int parID = 0; parID < Acts::eBoundSize; parID++) {
    std::string parName = m_cfg.paramNames.at(parID);
    double residual = trackParameter[parID] - truthParameter[parID];
    m_res.at(parName).fill({residual});
    m_resVsEta.at(parName).fill({truthEta, residual});
    m_resVsPhi.at(parName).fill({truthPhi, residual});
    m_resVsPt.at(parName).fill({truthPt, residual});

    if (!fittedParamters.covariance().has_value()) {
      ACTS_WARNING("Fitted track parameter :" << parName
                                              << " has no covariance");
      continue;
    }

    auto covariance = *fittedParamters.covariance();
    if (covariance(parID, parID) <= 0.0) {
      ACTS_WARNING("Fitted track parameter :"
                   << parName << " has non-positive covariance = "
                   << covariance(parID, parID));
      continue;
    }

    double pull = residual / std::sqrt(covariance(parID, parID));
    m_pull.at(parName).fill({pull});
    m_pullVsEta.at(parName).fill({truthEta, pull});
    m_pullVsPhi.at(parName).fill({truthPhi, pull});
    m_pullVsPt.at(parName).fill({truthPt, pull});
  }

  // relative pT resolution sigma(pT)/pT = p*sigma(q/p)
  double qop_residual = trackParameter[Acts::eBoundQOverP] - truthParameter[Acts::eBoundQOverP];
  double rel_pt_residual = std::abs(1. / trackParameter[Acts::eBoundQOverP]) * qop_residual;

  // for (const auto& it : m_res){
  //   std::cout<<"PF::KEY="<<it.first<<std::endl;
  // }

  m_res.at("res_pt_o_pt").fill({rel_pt_residual});
  m_resVsEta.at("res_pt_o_pt_vs_eta").fill({truthEta, rel_pt_residual});
  m_resVsPhi.at("res_pt_o_pt_vs_phi").fill({truthPhi, rel_pt_residual});
  m_resVsEtaVsPhi.at("res_pt_o_pt_vs_eta_vs_phi").fill({truthEta, truthPhi, rel_pt_residual});
}

} // namespace AliceActsTrk
