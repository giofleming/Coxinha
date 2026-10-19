// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "AO2DMaker.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace AliceActsTrk
{

namespace
{
constexpr int kNSectors = 18;
constexpr double kTwoPi = 2.0 * std::numbers::pi;
constexpr double kSectorSpan = kTwoPi / static_cast<double>(kNSectors);
constexpr double kSnpSafety = 1.e-6;
} // namespace

AO2DMaker::AO2DMaker(const Config& config,
                       std::unique_ptr<const Acts::Logger> logger)
  : IAlgorithm("AO2DMaker", std::move(logger)), m_cfg(config)
{
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("AO2DMaker: missing input tracks");
  }

  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("AO2DMaker: missing input particles");
  }

  if (m_cfg.outputFile.empty()) {
    throw std::invalid_argument("AO2DMaker: empty output file name");
  }

  if (!std::isfinite(m_cfg.positionScale) || m_cfg.positionScale <= 0.f) {
    throw std::invalid_argument("AO2DMaker: positionScale must be finite and positive");
  }

  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputParticles.initialize(m_cfg.inputParticles);

  mFile.reset(TFile::Open(m_cfg.outputFile.c_str(), "RECREATE", "ACTS AO2D"));
  if (!mFile || mFile->IsZombie()) {
    throw std::runtime_error("AO2DMaker: cannot create output file " + m_cfg.outputFile);
  }

  mFile->SetCompressionSettings(m_cfg.compressionSettings);

  // AO2D tables are stored inside dataframe directories. For a standalone
  // ACTS output file we keep all events in one dataframe and use global table
  // indices within that dataframe.
  mDataFrameDirectory = mFile->mkdir("DF_0");
  if (mDataFrameDirectory == nullptr) {
    throw std::runtime_error("AO2DMaker: cannot create DF_0 directory");
  }

  createTrees();
}

AO2DMaker::~AO2DMaker()
{
  std::scoped_lock lock(mOutputMutex);

  if (!mFile || mFile->IsZombie()) {
    return;
  }

  if (mDataFrameDirectory != nullptr) {
    mDataFrameDirectory->cd();

    if (mCollisionTree != nullptr) {
      mCollisionTree->Write();
    }
    if (mMCCollisionTree != nullptr) {
      mMCCollisionTree->Write();
    }
    if (mMCParticleTree != nullptr) {
      mMCParticleTree->Write();
    }
    if (mTrackTree != nullptr) {
      mTrackTree->Write();
    }
  }

  mFile->Close();
}

void AO2DMaker::createTrees()
{
  mDataFrameDirectory->cd();

  // ------------------------------------------------------------------
  // aod::Collisions_001
  // ------------------------------------------------------------------
  mCollisionTree = new TTree("O2collision_001", "Collision tree");

  mCollisionTree->Branch("fIndexBCs", &mCollisionRow.fIndexBCs, "fIndexBCs/I");

  mCollisionTree->Branch("fPosX", &mCollisionRow.fPosX, "fPosX/F");
  mCollisionTree->Branch("fPosY", &mCollisionRow.fPosY, "fPosY/F");
  mCollisionTree->Branch("fPosZ", &mCollisionRow.fPosZ, "fPosZ/F");

  mCollisionTree->Branch("fCovXX", &mCollisionRow.fCovXX, "fCovXX/F");
  mCollisionTree->Branch("fCovXY", &mCollisionRow.fCovXY, "fCovXY/F");
  mCollisionTree->Branch("fCovYY", &mCollisionRow.fCovYY, "fCovYY/F");
  mCollisionTree->Branch("fCovXZ", &mCollisionRow.fCovXZ, "fCovXZ/F");
  mCollisionTree->Branch("fCovYZ", &mCollisionRow.fCovYZ, "fCovYZ/F");
  mCollisionTree->Branch("fCovZZ", &mCollisionRow.fCovZZ, "fCovZZ/F");

  mCollisionTree->Branch("fFlags", &mCollisionRow.fFlags, "fFlags/s");
  mCollisionTree->Branch("fChi2", &mCollisionRow.fChi2, "fChi2/F");
  mCollisionTree->Branch("fNumContrib", &mCollisionRow.fNumContrib, "fNumContrib/s");
  mCollisionTree->Branch("fCollisionTime", &mCollisionRow.fCollisionTime, "fCollisionTime/F");
  mCollisionTree->Branch("fCollisionTimeRes", &mCollisionRow.fCollisionTimeRes, "fCollisionTimeRes/F");

  // ------------------------------------------------------------------
  // aod::McCollisions_001
  // ------------------------------------------------------------------
  mMCCollisionTree = new TTree("O2mccollision_001", "MC collisions");

  mMCCollisionTree->Branch("fIndexBCs", &mMCCollisionRow.fIndexBCs, "fIndexBCs/I");
  mMCCollisionTree->Branch("fGeneratorsID", &mMCCollisionRow.fGeneratorsID, "fGeneratorsID/S");

  mMCCollisionTree->Branch("fPosX", &mMCCollisionRow.fPosX, "fPosX/F");
  mMCCollisionTree->Branch("fPosY", &mMCCollisionRow.fPosY, "fPosY/F");
  mMCCollisionTree->Branch("fPosZ", &mMCCollisionRow.fPosZ, "fPosZ/F");

  mMCCollisionTree->Branch("fT", &mMCCollisionRow.fT, "fT/F");
  mMCCollisionTree->Branch("fWeight", &mMCCollisionRow.fWeight, "fWeight/F");
  mMCCollisionTree->Branch("fImpactParameter", &mMCCollisionRow.fImpactParameter, "fImpactParameter/F");
  mMCCollisionTree->Branch("fEventPlaneAngle", &mMCCollisionRow.fEventPlaneAngle, "fEventPlaneAngle/F");

  // ------------------------------------------------------------------
  // aod::StoredMcParticles_001
  // ------------------------------------------------------------------
  mMCParticleTree = new TTree("O2mcparticle_001", "Kinematics");

  mMCParticleTree->Branch("fIndexMcCollisions", &mMCParticleRow.fIndexMcCollisions, "fIndexMcCollisions/I");
  mMCParticleTree->Branch("fPdgCode", &mMCParticleRow.fPdgCode, "fPdgCode/I");
  mMCParticleTree->Branch("fStatusCode", &mMCParticleRow.fStatusCode, "fStatusCode/I");
  mMCParticleTree->Branch("fFlags", &mMCParticleRow.fFlags, "fFlags/b");

  mMCParticleTree->Branch("fIndexArray_Mothers_size",
                          &mMCParticleRow.fIndexArray_Mothers_size,
                          "fIndexArray_Mothers_size/I");
  mMCParticleTree->Branch("fIndexArray_Mothers",
                          mMCParticleRow.fIndexArray_Mothers,
                          "fIndexArray_Mothers[fIndexArray_Mothers_size]/I");
  mMCParticleTree->Branch("fIndexSlice_Daughters",
                          mMCParticleRow.fIndexSlice_Daughters,
                          "fIndexSlice_Daughters[2]/I");

  mMCParticleTree->Branch("fWeight", &mMCParticleRow.fWeight, "fWeight/F");

  mMCParticleTree->Branch("fPx", &mMCParticleRow.fPx, "fPx/F");
  mMCParticleTree->Branch("fPy", &mMCParticleRow.fPy, "fPy/F");
  mMCParticleTree->Branch("fPz", &mMCParticleRow.fPz, "fPz/F");
  mMCParticleTree->Branch("fE", &mMCParticleRow.fE, "fE/F");

  mMCParticleTree->Branch("fVx", &mMCParticleRow.fVx, "fVx/F");
  mMCParticleTree->Branch("fVy", &mMCParticleRow.fVy, "fVy/F");
  mMCParticleTree->Branch("fVz", &mMCParticleRow.fVz, "fVz/F");
  mMCParticleTree->Branch("fVt", &mMCParticleRow.fVt, "fVt/F");

  // ------------------------------------------------------------------
  // aod::StoredTracks
  // ------------------------------------------------------------------
  mTrackTree = new TTree("O2track", "Barrel tracks Parameters");

  mTrackTree->Branch("fIndexCollisions", &mTrackRow.fIndexCollisions, "fIndexCollisions/I");
  mTrackTree->Branch("fTrackType", &mTrackRow.fTrackType, "fTrackType/b");

  mTrackTree->Branch("fX", &mTrackRow.fX, "fX/F");
  mTrackTree->Branch("fAlpha", &mTrackRow.fAlpha, "fAlpha/F");
  mTrackTree->Branch("fY", &mTrackRow.fY, "fY/F");
  mTrackTree->Branch("fZ", &mTrackRow.fZ, "fZ/F");
  mTrackTree->Branch("fSnp", &mTrackRow.fSnp, "fSnp/F");
  mTrackTree->Branch("fTgl", &mTrackRow.fTgl, "fTgl/F");
  mTrackTree->Branch("fSigned1Pt", &mTrackRow.fSigned1Pt, "fSigned1Pt/F");
}

ProcessCode AO2DMaker::execute(const AlgorithmContext& ctx) const
{
  std::scoped_lock lock(mOutputMutex);

  const ConstTrackContainer& tracks = m_inputTracks(ctx);
  const SimParticleContainer& particles = m_inputParticles(ctx);

  const int32_t collisionId = mCollisionCounter;
  const int32_t mcCollisionId = mMCCollisionCounter;

  // ------------------------------------------------------------------
  // Event / collision position.
  //
  // The ACTS generator input used by this chain does not provide a separate
  // reconstructed-primary-vertex collection here. Use the generated vertex
  // of the first SimParticle as the event vertex. If there are no particles,
  // keep the origin.
  // ------------------------------------------------------------------
  Acts::Vector3 eventVertex = Acts::Vector3::Zero();
  if (!particles.empty()) {
    eventVertex = particles.begin()->position();
  }

  const float vx = static_cast<float>(m_cfg.positionScale * eventVertex.x());
  const float vy = static_cast<float>(m_cfg.positionScale * eventVertex.y());
  const float vz = static_cast<float>(m_cfg.positionScale * eventVertex.z());

  // ------------------------------------------------------------------
  // MC collision.
  // ------------------------------------------------------------------
  mMCCollisionRow = {};
  mMCCollisionRow.fIndexBCs = -1;
  mMCCollisionRow.fGeneratorsID = 0;
  mMCCollisionRow.fPosX = vx;
  mMCCollisionRow.fPosY = vy;
  mMCCollisionRow.fPosZ = vz;
  mMCCollisionRow.fT = 0.f;
  mMCCollisionRow.fWeight = 1.f;
  mMCCollisionRow.fImpactParameter = 0.f;
  mMCCollisionRow.fEventPlaneAngle = 0.f;
  mMCCollisionTree->Fill();

  // ------------------------------------------------------------------
  // MC particles.
  // ------------------------------------------------------------------
  for (const auto& particle : particles) {
    mMCParticleRow = {};

    mMCParticleRow.fIndexMcCollisions = mcCollisionId;
    mMCParticleRow.fPdgCode = particle.pdg();
    mMCParticleRow.fStatusCode = m_cfg.defaultHepMCStatus;

    // Keep ALICE-specific MC flags unset. In particular, the transport bit is
    // not set, so these rows behave as generator particles downstream.
    mMCParticleRow.fFlags = 0;

    // SimParticleContainer in this chain does not provide the full generator
    // genealogy needed by the AO2D mother-array / daughter-slice convention.
    mMCParticleRow.fIndexArray_Mothers_size = 0;
    mMCParticleRow.fIndexSlice_Daughters[0] = -1;
    mMCParticleRow.fIndexSlice_Daughters[1] = -1;

    mMCParticleRow.fWeight = 1.f;

    const auto& momentum = particle.momentum();
    const auto& position = particle.position();

    mMCParticleRow.fPx = static_cast<float>(momentum.x());
    mMCParticleRow.fPy = static_cast<float>(momentum.y());
    mMCParticleRow.fPz = static_cast<float>(momentum.z());

    // SimParticle::mass() is already available in the ACTS input used by the
    // LUT maker. Compute the on-shell energy instead of depending on a
    // version-specific SimParticle::energy() accessor.
    const double p2 = momentum.squaredNorm();
    const double mass = particle.mass();
    mMCParticleRow.fE = static_cast<float>(std::sqrt(std::max(0.0, p2 + mass * mass)));

    mMCParticleRow.fVx = static_cast<float>(m_cfg.positionScale * position.x());
    mMCParticleRow.fVy = static_cast<float>(m_cfg.positionScale * position.y());
    mMCParticleRow.fVz = static_cast<float>(m_cfg.positionScale * position.z());

    // Keep time at zero to avoid imposing a version-specific ACTS time-unit
    // convention. It can be replaced by particle.time() if needed.
    mMCParticleRow.fVt = 0.f;

    mMCParticleTree->Fill();
  }

  // ------------------------------------------------------------------
  // Reconstructed tracks.
  // ------------------------------------------------------------------
  std::size_t nWrittenTracks = 0;

  for (const auto& trackProxy : tracks) {
    if (!trackProxy.hasReferenceSurface()) {
      continue;
    }

    const Acts::BoundTrackParameters track = trackProxy.createParametersAtReference();

    mTrackRow = {};
    if (!fillTrackRow(track, ctx.simGeoContext, collisionId, mTrackRow)) {
      continue;
    }

    mTrackTree->Fill();
    ++nWrittenTracks;
  }

  // ------------------------------------------------------------------
  // Reconstructed collision.
  //
  // No standalone vertex fit is available in the supplied ACTS chain, so the
  // reconstructed collision position is initialized from the MC event vertex
  // and its covariance/chi2/time resolution are left at zero.
  // ------------------------------------------------------------------
  mCollisionRow = {};
  mCollisionRow.fIndexBCs = -1;
  mCollisionRow.fPosX = vx;
  mCollisionRow.fPosY = vy;
  mCollisionRow.fPosZ = vz;

  mCollisionRow.fFlags = 0;
  mCollisionRow.fChi2 = 0.f;
  mCollisionRow.fNumContrib = static_cast<uint16_t>(
    std::min<std::size_t>(nWrittenTracks, std::numeric_limits<uint16_t>::max()));
  mCollisionRow.fCollisionTime = 0.f;
  mCollisionRow.fCollisionTimeRes = 0.f;

  mCollisionTree->Fill();

  ++mCollisionCounter;
  ++mMCCollisionCounter;

  return ProcessCode::SUCCESS;
}

bool AO2DMaker::fillTrackRow(const Acts::BoundTrackParameters& track,
                              const Acts::GeometryContext& simGeoContext,
                              int32_t collisionId,
                              TrackRow& row) const
{
  const auto& parameters = track.parameters();
  const auto& surface = track.referenceSurface();

  const Acts::Vector2 localPosition{
    parameters[Acts::eBoundLoc0],
    parameters[Acts::eBoundLoc1]};

  const Acts::Vector3 direction = track.direction();
  const Acts::Vector3 globalPosition =
    surface.localToGlobal(simGeoContext, localPosition, direction);
  const Acts::Vector3 momentum = track.momentum();

  if (!globalPosition.allFinite() || !momentum.allFinite()) {
    return false;
  }

  const double px = momentum.x();
  const double py = momentum.y();
  const double pz = momentum.z();
  const double pt = std::hypot(px, py);

  if (!(pt > std::numeric_limits<double>::epsilon()) || !std::isfinite(pt)) {
    return false;
  }

  // Use the same O2-like 18-sector alpha convention used by the LUT code.
  const double phi = std::atan2(py, px);
  const double alpha = angleToSectorAlpha(phi);

  const double sinAlpha = std::sin(alpha);
  const double cosAlpha = std::cos(alpha);

  const double xGlobal = m_cfg.positionScale * globalPosition.x();
  const double yGlobal = m_cfg.positionScale * globalPosition.y();
  const double zGlobal = m_cfg.positionScale * globalPosition.z();

  // Global -> O2 local frame.
  const double xLocal = xGlobal * cosAlpha + yGlobal * sinAlpha;
  const double yLocal = -xGlobal * sinAlpha + yGlobal * cosAlpha;

  const double pxLocal = px * cosAlpha + py * sinAlpha;
  const double pyLocal = -px * sinAlpha + py * cosAlpha;

  // The O2 local representation assumes propagation along +local x.
  if (!(pxLocal > 0.0)) {
    return false;
  }

  double snp = pyLocal / pt;
  snp = std::clamp(snp, -1.0 + kSnpSafety, 1.0 - kSnpSafety);

  if (!std::isfinite(xLocal) || !std::isfinite(yLocal) ||
      !std::isfinite(zGlobal) || !std::isfinite(alpha) ||
      !std::isfinite(snp) || !std::isfinite(pz / pt)) {
    return false;
  }

  // BoundTrackParameters::qOverP() carries the charge sign. AO2D stores
  // sign(q)/pT rather than q/pT with the charge magnitude.
  const double qOverP = track.qOverP();
  if (!std::isfinite(qOverP)) {
    return false;
  }
  const double chargeSign = std::signbit(qOverP) ? -1.0 : 1.0;

  row.fIndexCollisions = collisionId;

  // No detector-specific O2 track type exists in this standalone ACTS chain.
  row.fTrackType = 0;

  row.fX = static_cast<float>(xLocal);
  row.fAlpha = static_cast<float>(alpha);
  row.fY = static_cast<float>(yLocal);
  row.fZ = static_cast<float>(zGlobal);
  row.fSnp = static_cast<float>(snp);
  row.fTgl = static_cast<float>(pz / pt);
  row.fSigned1Pt = static_cast<float>(chargeSign / pt);

  return true;
}

double AO2DMaker::normalizeAngle(double angle) noexcept
{
  angle = std::remainder(angle, kTwoPi);
  if (angle <= -std::numbers::pi) {
    angle += kTwoPi;
  } else if (angle > std::numbers::pi) {
    angle -= kTwoPi;
  }
  return angle;
}

int AO2DMaker::angleToSector(double phi) noexcept
{
  phi = normalizeAngle(phi);

  constexpr double radiansToDegrees = 180.0 / std::numbers::pi;
  constexpr double sectorSpanDegrees = 360.0 / static_cast<double>(kNSectors);

  int sector = static_cast<int>(phi * radiansToDegrees / sectorSpanDegrees);
  if (phi < 0.0) {
    sector += kNSectors - 1;
  }

  sector %= kNSectors;
  if (sector < 0) {
    sector += kNSectors;
  }

  return sector;
}

double AO2DMaker::sectorToAlpha(int sector) noexcept
{
  sector %= kNSectors;
  if (sector < 0) {
    sector += kNSectors;
  }

  return normalizeAngle(kSectorSpan * (static_cast<double>(sector) + 0.5));
}

double AO2DMaker::angleToSectorAlpha(double phi) noexcept
{
  return sectorToAlpha(angleToSector(phi));
}

} // namespace AliceActsTrk
