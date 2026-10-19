// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "FlatLutEntry.h"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsToO2Track.hpp"
#include "ResidualAccumulator.hpp"

#include <Eigen/Core>
#include <map>
#include <string>
#include <array>
#include <cstdint>
#include <vector>
#include <cmath>
#include <cstddef>
#include <optional>
#include <TTree.h>
class TH3F;
class TH1D;

using namespace ActsExamples;

namespace AliceActsTrk
{

struct LutContainer {
  o2::delphes::lutHeader_t header;
  std::vector<o2::delphes::lutEntry_t> entries;
};

/// Matches tracks to truth particles and vice versa
class LutMaker final : public IAlgorithm
{
 public:
  struct Config {
    std::string inputTracks;                /// Input (fitted) tracks collection
    std::string inputParticles;             /// Input particles collection.
    std::string inputTrackParticleMatching; /// Input track-particle matching.
    std::string inputParticleTrackMatching; /// Input particle-track matching.
    std::string lutTag;                     /// Tag to identify the LUT (e.g. "lut_v1")
    float magneticField = 2.0f;             /// Magnetic field value to store in LUT header (in Tesla)
    bool storeCovarianceSigma = false;      /// If true store per-bin sigma of covariance-element distributions instead of means
    bool produceOutputTree = false;         /// If true produce an output tree of particles and tracks
    bool useResiduals = true;               /// If true use residuals to compute covariance, otherwise use direct covariance sums
    std::string outputDirectory = "luts";   /// Output directory for LUTs and tree
  };

  LutMaker(const Config& config, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ~LutMaker();

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<MeasurementParticlesMap> m_inputMeasurementParticlesMap{this, "InputMeasurementParticlesMap"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{this, "InputTrackParticleMatching"};
  ReadDataHandle<ParticleTrackMatching> m_inputParticleTrackMatching{this, "InputParticleTrackMatching"};

  static std::map<int, LutContainer> mLutMap;                          ///< Map of LUTs indexed by PDG code
  static std::map<int, std::vector<std::uint64_t>> mPtTruthCounts;     ///< Truth particle counts per PDG and pT bin
  static std::map<int, std::vector<std::uint64_t>> mPtMatchedCounts;   ///< Matched track counts per PDG and pT bin
  static std::map<int, std::vector<std::array<double, 15>>> mCovSums;  ///< Running sum of covariance lower-triangle elements per PDG and entry index
  static std::map<int, std::vector<std::array<double, 15>>> mCov2Sums; ///< Running sum of squared covariance lower-triangle elements per PDG and entry index
  static std::map<int, std::vector<ResidualAccumulator>> mResiduals;   ///< Running residual accumulator per PDG and LUT entry index
  static std::map<int, std::vector<std::uint64_t>> mEntryCounts;       ///< Number of accumulated tracks per PDG and entry index

  static std::map<int, TH3F*> res_acts_d0_vs_pT; ///< ACTS-native d0 residual at the track reference surface vs pT and eta
  static std::map<int, TH3F*> res_o2_y_vs_pT;    ///< Current O2 local-Y residual used by the residual-based LUT vs pT and eta
  static std::map<int, TH3F*> res_qopt_vs_pT;    ///< Resolution graphs for q/pT vs pT per PDG
  static std::map<int, TH3F*> res_pt_vs_pT;      ///< Resolution graphs for pT vs pT per PDG
  static TH1D* h_stat;                           ///< Statistics

  static TTree* mTree; ///< Output tree of particles and tracks
  struct TreeEntry {
    int pdg = 0;
    float px = 0.f;
    float pxTrue = 0.f;
    float py = 0.f;
    float pyTrue = 0.f;
    float pz = 0.f;
    float pzTrue = 0.f;
    void import(const SimParticle* p, const Acts::BoundTrackParameters* t)
    {
      pdg = p->pdg();
      px = t->momentum().x();
      pxTrue = p->momentum().x();
      py = t->momentum().y();
      pyTrue = p->momentum().y();
      pz = t->momentum().z();
      pzTrue = p->momentum().z();
      mTree->Fill();
    }
  };
  static TreeEntry mTreeContainer;
};

} // namespace AliceActsTrk
