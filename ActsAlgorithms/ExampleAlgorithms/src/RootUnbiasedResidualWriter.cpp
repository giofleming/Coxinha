// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RootUnbiasedResidualWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsPlugins/Root/HistogramConverter.hpp"

#include <cmath>
#include <ios>
#include <limits>
#include <stdexcept>
#include <utility>

#include <TFile.h>
#include <TH1.h>

using namespace ActsExamples;
namespace AliceActsTrk {

RootUnbiasedResidualWriter::RootUnbiasedResidualWriter(
    const RootUnbiasedResidualWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT(config.inputTracks, "RootUnbiasedResidualWriter", level),
      m_cfg(config),
      m_plotTool(m_cfg.plotToolConfig, level) {
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }
}

RootUnbiasedResidualWriter::~RootUnbiasedResidualWriter() {
  delete m_outputFile;
}

ProcessCode RootUnbiasedResidualWriter::finalize() {
  if (m_outputFile == nullptr) {
    return ProcessCode::SUCCESS;
  }

  m_outputFile->cd();

  using ActsPlugins::toRoot;
  for (const auto& [name, hist] : m_plotTool.resLoc0()) {
    toRoot(hist)->Write();
  }
  for (const auto& [name, hist] : m_plotTool.pullLoc0()) {
    toRoot(hist)->Write();
  }
  for (const auto& [name, hist] : m_plotTool.resLoc1()) {
    toRoot(hist)->Write();
  }
  for (const auto& [name, hist] : m_plotTool.pullLoc1()) {
    toRoot(hist)->Write();
  }

  ACTS_INFO("Wrote unbiased residual histograms to '"
           << m_outputFile->GetPath() << "'");

  m_outputFile->Close();
  return ProcessCode::SUCCESS;
}

ProcessCode RootUnbiasedResidualWriter::writeT(
    const AlgorithmContext& /*ctx*/, const ConstTrackContainer& tracks) {
  constexpr double nan = std::numeric_limits<double>::quiet_NaN();

  // Exclusive access to the histograms while filling
  std::lock_guard<std::mutex> lock(m_writeMutex);

  for (const auto& track : tracks) {
    for (const auto& state : track.trackStatesReversed()) {
      const auto flag = state.typeFlags();
      if (!flag.isMeasurement() && !flag.isOutlier()) {
        continue;
      }
      if (!state.hasSmoothed() || !state.hasProjector()) {
        continue;
      }

      // Measurements are part of the fit, so the smoothed parameters are
      // biased by them: recompute the unbiased parameters with this
      // measurement removed. Outliers are already excluded from the fit, so
      // the smoothed parameters are unbiased already.
      const bool evaluateUnbiased = flag.isMeasurement();

      std::pair<Acts::BoundVector, Acts::BoundMatrix> paramCov;
      if (evaluateUnbiased) {
        paramCov = Acts::calculateUnbiasedParametersCovariance(
            Acts::AnyConstTrackStateProxy{state});
      } else {
        paramCov = {state.smoothed(), state.smoothedCovariance()};
      }
      const auto& [parameters, covariance] = paramCov;

      const Acts::DynamicMatrix H =
          state.projectorSubspaceHelper().fullProjector().topLeftCorner(
              state.calibratedSize(), Acts::eBoundSize);
      const Acts::DynamicMatrix V = state.effectiveCalibratedCovariance();
      const Acts::DynamicMatrix resCov = V + H * covariance * H.transpose();
      const Acts::DynamicVector res =
          state.effectiveCalibrated() - H * parameters;

      const double resLoc0 = res[Acts::eBoundLoc0] / Acts::UnitConstants::um;
      const double pullLoc0 =
          resCov(Acts::eBoundLoc0, Acts::eBoundLoc0) > 0
              ? res[Acts::eBoundLoc0] /
                    std::sqrt(resCov(Acts::eBoundLoc0, Acts::eBoundLoc0))
              : nan;

      const bool hasLoc1 = state.calibratedSize() >= 2;
      double resLoc1 = 0.;
      double pullLoc1 = 0.;
      if (hasLoc1) {
        resLoc1 = res[Acts::eBoundLoc1] / Acts::UnitConstants::um;
        pullLoc1 =
            resCov(Acts::eBoundLoc1, Acts::eBoundLoc1) > 0
                ? res[Acts::eBoundLoc1] /
                      std::sqrt(resCov(Acts::eBoundLoc1, Acts::eBoundLoc1))
                : nan;
      }

      m_plotTool.fill(state.referenceSurface().geometryId(), resLoc0,
                      pullLoc0, hasLoc1, resLoc1, pullLoc1);
    }
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
