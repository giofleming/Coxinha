// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Histogram.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace ActsExamples;

namespace AliceActsTrk
{

/// Tools to make hists to show residual, i.e. smoothed_parameter -
/// truth_parameter, and pull, i.e. (smoothed_parameter -
/// truth_parameter)/smoothed_paramter_error, of track parameters at perigee
/// surface
class ResPlotTool
{
 public:
  using AxisVariant = Acts::Experimental::AxisVariant;
  using BoostRegularAxis = Acts::Experimental::BoostRegularAxis;
  using Histogram1 = Acts::Experimental::Histogram1;
  using Histogram2 = Acts::Experimental::Histogram2;
  using Histogram3 = Acts::Experimental::Histogram<3>;

  /// @brief Nested configuration struct
  struct Config {
    /// parameter sets to do plots
    std::vector<std::string> paramNames = {"d0", "z0", "phi",
                                           "theta", "qop", "t"};

    /// Binning info for variables
    std::map<std::string, AxisVariant> varBinning = {
      {"Eta", BoostRegularAxis(40, -4, 4, "#eta")},
      {"Phi", BoostRegularAxis(70, -3.5, 3.5, "#phi")},
      {"Pt", BoostRegularAxis(40, 0, 100, "pT [GeV/c]")},
      {"Pull", BoostRegularAxis(100, -5, 5, "pull")},
      {"Residual_d0", BoostRegularAxis(100, -0.5, 0.5, "r_{d0} [mm]")},
      {"Residual_z0", BoostRegularAxis(100, -0.5, 0.5, "r_{z0} [mm]")},
      {"Residual_phi", BoostRegularAxis(100, -0.01, 0.01, "r_{#phi} [rad]")},
      {"Residual_theta",
       BoostRegularAxis(100, -0.01, 0.01, "r_{#theta} [rad]")},
      {"Residual_pt_o_pt", BoostRegularAxis(100, -0.01, 0.01, "r_{p_T} / p_T")},
      {"Residual_qop", BoostRegularAxis(100, -0.1, 0.1, "r_{q/p} [c/GeV]")},
      {"Residual_t", BoostRegularAxis(100, -1000, 1000, "r_{t} [s]")}};

    std::vector<float> etaBins{0., 0.1, 0.5, 1.5, 2.0, 2.5, 4.0};
    std::vector<float> phiBins{0., 0.1, 0.5, 1.5, 2.0, 2.5, 4.0};
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param level Message level declaration
  ResPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief fill the histograms
  ///
  /// @param gctx the geometry context
  /// @param truthParticle the truth particle
  /// @param fittedParamters the fitted parameters at perigee surface
  void fill(const Acts::GeometryContext& gctx,
            const SimParticleState& truthParticle,
            const Acts::BoundTrackParameters& fittedParamters);

  /// @brief Accessors for histograms (const reference)
  const std::map<std::string, Histogram1>& res() const { return m_res; }
  const std::map<std::string, Histogram2>& resVsEta() const
  {
    return m_resVsEta;
  }
  const std::map<std::string, Histogram2>& resVsPhi() const
  {
    return m_resVsPhi;
  }
  const std::map<std::string, Histogram3>& resVsEtaVsPhi() const
  {
    return m_resVsEtaVsPhi;
  }
  const std::map<std::string, Histogram2>& resVsPt() const { return m_resVsPt; }
  const std::map<std::string, Histogram1>& pull() const { return m_pull; }
  const std::map<std::string, Histogram2>& pullVsEta() const
  {
    return m_pullVsEta;
  }
  const std::map<std::string, Histogram2>& pullVsPhi() const
  {
    return m_pullVsPhi;
  }
  const std::map<std::string, Histogram2>& pullVsPt() const
  {
    return m_pullVsPt;
  }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Residual distribution
  std::map<std::string, Histogram1> m_res;
  /// Residual vs eta scatter plot
  std::map<std::string, Histogram2> m_resVsEta;
  /// Residual vs phi scatter plot
  std::map<std::string, Histogram2> m_resVsPhi;
  /// Residual vs pT scatter plot
  std::map<std::string, Histogram2> m_resVsPt;

  /// Residual vs eta vs phi scatter plot
  std::map<std::string, Histogram3> m_resVsEtaVsPhi;

  /// Pull distribution
  std::map<std::string, Histogram1> m_pull;
  /// Pull vs eta scatter plot
  std::map<std::string, Histogram2> m_pullVsEta;
  /// Pull vs eta scatter plot
  std::map<std::string, Histogram2> m_pullVsPhi;
  /// Pull vs pT scatter plot
  std::map<std::string, Histogram2> m_pullVsPt;
};

} // namespace AliceActsTrk
