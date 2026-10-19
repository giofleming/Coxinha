// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "UnbiasedResidualPlotTool.hpp"

#include <mutex>
#include <string>

class TFile;

using namespace ActsExamples;
namespace AliceActsTrk {

/// Write out histograms of the unbiased hit-on-track residual and pull,
/// i.e. calibrated measurement - unbiased track parameter, where the
/// unbiased parameter is the fitted track state recomputed with this
/// particular measurement removed from the fit.
///
/// Unlike RootTrackStatesWriter, this writer needs no truth information
/// (no simulated hits, no particles, no truth matching) since the unbiased
/// residual is purely a measurement-vs-fit comparison. It is therefore
/// usable directly on measurements without an associated truth record, and
/// its output (a handful of histograms) stays small regardless of the
/// number of processed events, unlike a per-hit ntuple.
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootUnbiasedResidualWriter final : public WriterT<ConstTrackContainer> {
 public:
  struct Config {
    /// Input (fitted) tracks collection.
    std::string inputTracks;
    /// Output filename.
    std::string filePath = "unbiased_residuals.root";
    /// File access mode.
    std::string fileMode = "RECREATE";
    /// Plot tool configuration.
    UnbiasedResidualPlotTool::Config plotToolConfig;
  };

  /// Construct from configuration and log level.
  /// @param config The configuration
  /// @param level The logger level
  RootUnbiasedResidualWriter(const Config& config, Acts::Logging::Level level);

  ~RootUnbiasedResidualWriter() override;

  /// End-of-run hook: write the booked histograms to file.
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const ConstTrackContainer& tracks) override;

  Config m_cfg;

  /// Mutex used to protect multi-threaded filling.
  std::mutex m_writeMutex;

  TFile* m_outputFile{nullptr};

  UnbiasedResidualPlotTool m_plotTool;
};

}  // namespace ActsExamples
