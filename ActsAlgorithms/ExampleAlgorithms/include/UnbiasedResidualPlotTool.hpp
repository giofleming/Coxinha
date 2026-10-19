// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Histogram.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <map>
#include <memory>
#include <string>

//using namespace ActsExamples;

namespace AliceActsTrk {

/// Tool to make hists to show the unbiased hit-on-track residual, i.e.
/// calibrated measurement - unbiased track parameter (the fitted parameter
/// recomputed with this particular measurement removed from the fit), and
/// the associated pull.
///
/// Histograms are booked lazily, once per distinct geometry category
/// encountered, plus one inclusive "all" category collecting every hit.
class UnbiasedResidualPlotTool {
 public:
  using AxisVariant = Acts::Experimental::AxisVariant;
  using BoostRegularAxis = Acts::Experimental::BoostRegularAxis;
  using Histogram1 = Acts::Experimental::Histogram1;

  /// @brief Nested configuration struct
  struct Config {
    /// Binning info for variables
    std::map<std::string, AxisVariant> varBinning = {
        {"ResidualLoc0", BoostRegularAxis(100, -0.5, 0.5, "res loc0 [mm]")},
        {"ResidualLoc1", BoostRegularAxis(100, -0.5, 0.5, "res loc1 [mm]")},
        {"Pull", BoostRegularAxis(100, -5, 5, "pull")}};
  };

  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  UnbiasedResidualPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// Fill the residual/pull histograms for one measurement.
  ///
  /// @param geoId geometry identifier of the surface the hit is on
  /// @param resLoc0 unbiased residual in the first local coordinate
  /// @param pullLoc0 unbiased pull in the first local coordinate
  /// @param hasLoc1 whether the measurement also constrains the second
  ///   local coordinate
  /// @param resLoc1 unbiased residual in the second local coordinate
  ///   (only used if @p hasLoc1 is true)
  /// @param pullLoc1 unbiased pull in the second local coordinate
  ///   (only used if @p hasLoc1 is true)
  void fill(Acts::GeometryIdentifier geoId, double resLoc0, double pullLoc0,
            bool hasLoc1, double resLoc1, double pullLoc1);

  /// Category key ("all", or "vol{V}_lay{L}") to booked histogram.
  const std::map<std::string, Histogram1>& resLoc0() const {
    return m_resLoc0;
  }
  const std::map<std::string, Histogram1>& pullLoc0() const {
    return m_pullLoc0;
  }
  const std::map<std::string, Histogram1>& resLoc1() const {
    return m_resLoc1;
  }
  const std::map<std::string, Histogram1>& pullLoc1() const {
    return m_pullLoc1;
  }

 private:
  Config m_cfg;

  std::unique_ptr<const Acts::Logger> m_logger;

  std::map<std::string, Histogram1> m_resLoc0;
  std::map<std::string, Histogram1> m_pullLoc0;
  std::map<std::string, Histogram1> m_resLoc1;
  std::map<std::string, Histogram1> m_pullLoc1;

  const Acts::Logger& logger() const { return *m_logger; }

  /// Geometry category key for a geometry identifier: "vol{V}_lay{L}_sens{S}"
  static std::string category(Acts::GeometryIdentifier geoId);

  /// Return the histogram for `category`, booking it on first use.
  Histogram1& bookedHistogram(std::map<std::string, Histogram1>& hists,
                              const std::string& category,
                              const std::string& namePrefix,
                              const std::string& titlePrefix,
                              const AxisVariant& axis);

  void fillOne(std::map<std::string, Histogram1>& hists,
              const std::string& category, const std::string& namePrefix,
              const std::string& titlePrefix, const AxisVariant& axis,
              double value);
  };

}  // namespace AliceActsTrk
