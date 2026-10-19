// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "UnbiasedResidualPlotTool.hpp"

#include <format>

//using namespace ActsExamples;

namespace AliceActsTrk {

UnbiasedResidualPlotTool::UnbiasedResidualPlotTool(
    const UnbiasedResidualPlotTool::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("UnbiasedResidualPlotTool",
                                                  lvl)) {
  ACTS_DEBUG(
      "Booking the inclusive histograms for unbiased residual and pull "
      "plots; per-geometry histograms are booked lazily.");

  // Book the inclusive "all" category up front so it always exists even for
  // events without any hits.
  bookedHistogram(m_resLoc0, "all", "res_loc0", "Unbiased residual of loc0",
                  m_cfg.varBinning.at("ResidualLoc0"));
  bookedHistogram(m_pullLoc0, "all", "pull_loc0", "Unbiased pull of loc0",
                  m_cfg.varBinning.at("Pull"));
  bookedHistogram(m_resLoc1, "all", "res_loc1", "Unbiased residual of loc1",
                  m_cfg.varBinning.at("ResidualLoc1"));
  bookedHistogram(m_pullLoc1, "all", "pull_loc1", "Unbiased pull of loc1",
                  m_cfg.varBinning.at("Pull"));
}

std::string UnbiasedResidualPlotTool::category(Acts::GeometryIdentifier geoId) {
  return std::format("vol{}_lay{}", geoId.volume(), geoId.layer());
}

UnbiasedResidualPlotTool::Histogram1&
UnbiasedResidualPlotTool::bookedHistogram(
    std::map<std::string, Histogram1>& hists, const std::string& category,
    const std::string& namePrefix, const std::string& titlePrefix,
    const AxisVariant& axis) {
  auto [it, inserted] = hists.try_emplace(
      category,
      Histogram1(std::format("{}_{}", namePrefix, category),
                std::format("{} ({})", titlePrefix, category),
                std::array{axis}));
  if (inserted) {
    ACTS_VERBOSE("Booked histogram '" << it->second.name() << "'");
  }
  return it->second;
}

void UnbiasedResidualPlotTool::fillOne(
    std::map<std::string, Histogram1>& hists, const std::string& category,
    const std::string& namePrefix, const std::string& titlePrefix,
    const AxisVariant& axis, double value) {
  bookedHistogram(hists, "all", namePrefix, titlePrefix, axis)
      .fill({value});
  bookedHistogram(hists, category, namePrefix, titlePrefix, axis)
      .fill({value});
}

void UnbiasedResidualPlotTool::fill(Acts::GeometryIdentifier geoId,
                                    double resLoc0, double pullLoc0,
                                    bool hasLoc1, double resLoc1,
                                    double pullLoc1) {
  const std::string cat = category(geoId);

  fillOne(m_resLoc0, cat, "res_loc0", "Unbiased residual of loc0",
         m_cfg.varBinning.at("ResidualLoc0"), resLoc0);
  fillOne(m_pullLoc0, cat, "pull_loc0", "Unbiased pull of loc0",
         m_cfg.varBinning.at("Pull"), pullLoc0);

  if (hasLoc1) {
    fillOne(m_resLoc1, cat, "res_loc1", "Unbiased residual of loc1",
           m_cfg.varBinning.at("ResidualLoc1"), resLoc1);
    fillOne(m_pullLoc1, cat, "pull_loc1", "Unbiased pull of loc1",
           m_cfg.varBinning.at("Pull"), pullLoc1);
  }
}

}  // namespace ActsExamples
