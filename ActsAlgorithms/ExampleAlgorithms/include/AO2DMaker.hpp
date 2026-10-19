// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <cstdint>
#include <memory>
#include <mutex>
#include <string>

class TDirectory;
class TFile;
class TTree;

using namespace ActsExamples;

namespace AliceActsTrk
{

/// Minimal standalone AO2D maker for ACTS simulation/reconstruction output.
///
/// The maker has no dependency on AliceO2. It reproduces the stored ROOT
/// layout of the relevant AO2D tables directly:
///
///   DF_0/O2collision_001
///   DF_0/O2mccollision_001
///   DF_0/O2mcparticle_001
///   DF_0/O2track
///
/// One reconstructed collision and one MC collision are written for every
/// ACTS event. All input SimParticles are written as MC particles and all
/// fitted tracks with a valid reference surface are written as tracks.
class AO2DMaker final : public IAlgorithm
{
 public:
  struct Config {
    std::string inputTracks;    ///< Input fitted-track collection.
    std::string inputParticles; ///< Input simulated-particle collection.

    std::string outputFile = "AO2D.root"; ///< Output AO2D ROOT file.

    /// Conversion applied to ACTS positions before storing them in AO2D.
    /// ACTS examples normally use mm, while O2 track/vertex positions are cm.
    float positionScale = 0.1f;

    /// HepMC status stored for all SimParticles. A plain value of 1 represents
    /// a final-state particle in the backward-compatible AO2D encoding.
    int defaultHepMCStatus = 1;

    /// ROOT compression setting passed to TFile::SetCompressionSettings().
    int compressionSettings = 101;
  };

  AO2DMaker(const Config& config,
            std::unique_ptr<const Acts::Logger> logger = nullptr);

  ~AO2DMaker() override;

  AO2DMaker(const AO2DMaker&) = delete;
  AO2DMaker& operator=(const AO2DMaker&) = delete;

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

 private:
  // ------------------------------------------------------------------
  // Local copies of the stored AO2D table rows.
  // ------------------------------------------------------------------

  struct CollisionRow {
    int32_t fIndexBCs = -1;

    float fPosX = 0.f;
    float fPosY = 0.f;
    float fPosZ = 0.f;

    float fCovXX = 0.f;
    float fCovXY = 0.f;
    float fCovYY = 0.f;
    float fCovXZ = 0.f;
    float fCovYZ = 0.f;
    float fCovZZ = 0.f;

    uint16_t fFlags = 0;
    float fChi2 = 0.f;
    uint16_t fNumContrib = 0;

    float fCollisionTime = 0.f;
    float fCollisionTimeRes = 0.f;
  };

  struct MCCollisionRow {
    int32_t fIndexBCs = -1;
    int16_t fGeneratorsID = 0;

    float fPosX = 0.f;
    float fPosY = 0.f;
    float fPosZ = 0.f;

    float fT = 0.f;
    float fWeight = 1.f;
    float fImpactParameter = 0.f;
    float fEventPlaneAngle = 0.f;
  };

  /// The current AO2D MC-particle table stores a variable-size mother array.
  /// We do not have generator genealogy in the ACTS collection used here, so
  /// this maker keeps the mother list empty and the daughter slice invalid.
  /// The fixed storage below is only a ROOT branch backing buffer.
  static constexpr int kMaxMothers = 8;

  struct MCParticleRow {
    int32_t fIndexMcCollisions = -1;

    int32_t fPdgCode = 0;
    int32_t fStatusCode = 0;
    uint8_t fFlags = 0;

    int32_t fIndexArray_Mothers_size = 0;
    int32_t fIndexArray_Mothers[kMaxMothers]{};
    int32_t fIndexSlice_Daughters[2] = {-1, -1};

    float fWeight = 1.f;

    float fPx = 0.f;
    float fPy = 0.f;
    float fPz = 0.f;
    float fE = 0.f;

    float fVx = 0.f;
    float fVy = 0.f;
    float fVz = 0.f;
    float fVt = 0.f;
  };

  struct TrackRow {
    int32_t fIndexCollisions = -1;
    uint8_t fTrackType = 0;

    float fX = 0.f;
    float fAlpha = 0.f;
    float fY = 0.f;
    float fZ = 0.f;
    float fSnp = 0.f;
    float fTgl = 0.f;
    float fSigned1Pt = 0.f;
  };

  void createTrees();

  /// Convert the fitted ACTS state at its reference surface into the stored
  /// O2 five-parameter representation. Returns false for unusable states.
  bool fillTrackRow(const Acts::BoundTrackParameters& track,
                    const Acts::GeometryContext& simGeoContext,
                    int32_t collisionId,
                    TrackRow& row) const;

  static double normalizeAngle(double angle) noexcept;
  static int angleToSector(double phi) noexcept;
  static double sectorToAlpha(int sector) noexcept;
  static double angleToSectorAlpha(double phi) noexcept;

  Config m_cfg;

  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};

  std::unique_ptr<TFile> mFile;
  TDirectory* mDataFrameDirectory = nullptr;

  TTree* mCollisionTree = nullptr;
  TTree* mMCCollisionTree = nullptr;
  TTree* mMCParticleTree = nullptr;
  TTree* mTrackTree = nullptr;

  mutable CollisionRow mCollisionRow;
  mutable MCCollisionRow mMCCollisionRow;
  mutable MCParticleRow mMCParticleRow;
  mutable TrackRow mTrackRow;

  mutable int32_t mCollisionCounter = 0;
  mutable int32_t mMCCollisionCounter = 0;

  /// ROOT TTree::Fill is not thread-safe. ACTS may execute algorithms for
  /// several events concurrently, therefore serialize access to the file.
  mutable std::mutex mOutputMutex;
};

} // namespace AliceActsTrk
