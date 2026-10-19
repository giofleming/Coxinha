// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "AliceVertexTruthMatcher.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Utilities/VertexTruthUtility.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <stdexcept>
#include <vector>

namespace AliceActsTrk {

namespace {

// sqrt of a covariance diagonal entry, guarded against (small) negatives.
double safeSigma(double variance) {
  return variance > 0. ? std::sqrt(variance) : 0.;
}

}  // namespace

AliceVertexTruthMatcher::AliceVertexTruthMatcher(const Config& config,
                                       std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("AliceVertexTruthMatcher", std::move(logger)), m_cfg(config) {
  if (m_cfg.inputVertices.empty()) {
    throw std::invalid_argument("Collection with vertices missing");
  }
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Collection with tracks missing");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Collection with particles missing");
  }
  if (m_cfg.inputTrackParticleMatching.empty()) {
    throw std::invalid_argument(
        "Collection with track-particle matching missing");
  }
  if (m_cfg.outputVertexTruthMatching.empty()) {
    throw std::invalid_argument("Missing output vertex-truth matching");
  }
  if (m_cfg.outputTruthVertexMatching.empty()) {
    throw std::invalid_argument("Missing output truth-vertex matching");
  }

  m_inputVertices.initialize(m_cfg.inputVertices);
  m_inputTracks.maybeInitialize(m_cfg.inputTracks);
  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputTrackParticleMatching.maybeInitialize(
      m_cfg.inputTrackParticleMatching);
  // Optional: only needed for the truth-vertex printout.
  m_inputTruthVertices.maybeInitialize(m_cfg.inputTruthVertices);
  m_outputVertexTruthMatching.initialize(m_cfg.outputVertexTruthMatching);
  m_outputTruthVertexMatching.initialize(m_cfg.outputTruthVertexMatching);
}

ProcessCode AliceVertexTruthMatcher::execute(const AlgorithmContext& ctx) const {
  // Read input reconstructed vertices
  const auto& vertices = m_inputVertices(ctx);
  // Read input tracks
  const auto& tracks = m_inputTracks(ctx);
  // Read truth particle input collection
  const SimParticleContainer& particles = m_inputParticles(ctx);
  // Read track-particle matching
  const TrackParticleMatching& trackParticleMatching =
      m_inputTrackParticleMatching(ctx);

  const auto& gctx = ctx.recoGeoContext;

  // Map a particle's barcode to the truth-vertex key used for matching. With
  // collapseSecondaryVertices, decay daughters are folded onto their primary
  // pp interaction (key = vertexPrimary only); otherwise the full vertexId
  // (primary+secondary+generation) is used, so each decay vertex is distinct.
  auto vtxKey = [&](const SimBarcode& bc) -> SimVertexBarcode {
    if (m_cfg.collapseSecondaryVertices) {
      return SimVertexBarcode{SimBarcode{}.withVertexPrimary(bc.vertexPrimary())};
    }
    return SimVertexBarcode{bc};
  };

  // Lookup of truth-vertex time (native units) by vertex key, used for the
  // optional time-consistency cut. Only available if inputTruthVertices is set.
  // When collapsing, key by primary and store the PRIMARY vertex's time.
  std::map<SimVertexBarcode, double> truthVtxTime;
  const bool haveTruthVertices = m_inputTruthVertices.isInitialized();
  if (haveTruthVertices) {
    for (const auto& sv : m_inputTruthVertices(ctx)) {
      if (m_cfg.collapseSecondaryVertices) {
        if (sv.vertexId().vertexSecondary() == 0) {  // the primary pp vertex
          truthVtxTime[vtxKey(sv.vertexId().barcode())] = sv.time();
        }
      } else {
        truthVtxTime[sv.vertexId()] = sv.time();
      }
    }
  }
  const bool applyTimeCut = m_cfg.useTimeWindow && haveTruthVertices;

  // ----------------------------------------------------------------------
  //  Per-track diagnostics collected during the matching loop, so we can
  //  print them after the (split-vertex) classification has been finalised.
  // ----------------------------------------------------------------------
  struct TrackDiag {
    Acts::Vector4 pos4 = Acts::Vector4::Zero();  // global x,y,z,t at perigee
    double sigD0 = 0.;                            // perigee loc0 uncertainty
    double sigZ0 = 0.;                            // perigee loc1 uncertainty
    double sigT = 0.;                             // perigee time uncertainty
    Acts::Vector4 oPos4 = Acts::Vector4::Zero();  // original (pre-fit) perigee x,y,z,t
    double oSigD0 = 0.;
    double oSigZ0 = 0.;
    double oSigT = 0.;
    double pt = 0.;
    double weight = 0.;
    bool hasParticle = false;
    std::optional<SimVertexBarcode> truthVtx;
    SimBarcode particleId{};
    int pdg = 0;
    bool timeCompatible = true;  // passed the optional time-window cut
    bool haveDt = false;         // truth-vertex time known -> dt meaningful
    double dt = 0.;              // t_track(orig) - t_truthVtx (native units)
  };
  struct VtxDiag {
    Acts::Vector4 pos4 = Acts::Vector4::Zero();
    Acts::Vector4 sigma4 = Acts::Vector4::Zero();
    std::vector<TrackDiag> tracks;
    int nSelected = 0;       // tracks above minTrkWeight
    int nWithParticle = 0;   // of which matched to a truth particle
    int nMatchedVtx = 0;     // of which from the matched (majority) truth vertex
    int nContamination = 0;  // of which from a different truth vertex
    int nTimeRejected = 0;   // of which rejected by the time-window cut
  };
  std::vector<VtxDiag> diags;
  diags.reserve(vertices.size());

  std::vector<VertexToTruthMatching> recoToTruthMatching;
  std::map<SimVertexBarcode, VertexToRecoMatching> truthToRecoMatching;

  // Do truth matching for each reconstructed vertex
  for (const auto& [vtxIndex, vtx] : Acts::enumerate(vertices)) {
    // Reconstructed tracks that contribute to the reconstructed vertex
    const std::vector<Acts::TrackAtVertex>& tracksAtVtx = vtx.tracks();

    // Containers for storing truth particles and truth vertices that
    // contribute to the reconstructed vertex
    std::vector<std::pair<SimVertexBarcode, double>> contributingTruthVertices;

    VtxDiag diag;
    diag.pos4 = vtx.fullPosition();
    const Acts::SquareMatrix4& cov = vtx.fullCovariance();
    diag.sigma4 = Acts::Vector4(safeSigma(cov(0, 0)), safeSigma(cov(1, 1)),
                                safeSigma(cov(2, 2)), safeSigma(cov(3, 3)));

    double totalTrackWeight = 0;
    for (const Acts::TrackAtVertex& trk : tracksAtVtx) {
      if (trk.trackWeight < m_cfg.minTrkWeight) {
        continue;
      }

      totalTrackWeight += trk.trackWeight;
      ++diag.nSelected;

      // ---- collect this track's diagnostics ----
      TrackDiag td;
      td.weight = trk.trackWeight;
      td.pos4 = trk.fittedParams.fourPosition(gctx);
      td.pt = trk.fittedParams.transverseMomentum();
      if (trk.fittedParams.covariance().has_value()) {
        const auto& tcov = *trk.fittedParams.covariance();
        td.sigD0 = safeSigma(tcov(Acts::eBoundLoc0, Acts::eBoundLoc0));
        td.sigZ0 = safeSigma(tcov(Acts::eBoundLoc1, Acts::eBoundLoc1));
        td.sigT = safeSigma(tcov(Acts::eBoundTime, Acts::eBoundTime));
      }

      // Original (pre-fit) input track parameters, at their own reference
      // surface (beamline perigee). This is the track's *own* d0/z0/t0, as
      // opposed to fittedParams which are expressed at the vertex.
      const auto origParams =
          Acts::InputTrack::extractParameters(trk.originalParams);
      td.oPos4 = origParams.fourPosition(gctx);
      if (origParams.covariance().has_value()) {
        const auto& ocov = *origParams.covariance();
        td.oSigD0 = safeSigma(ocov(Acts::eBoundLoc0, Acts::eBoundLoc0));
        td.oSigZ0 = safeSigma(ocov(Acts::eBoundLoc1, Acts::eBoundLoc1));
        td.oSigT = safeSigma(ocov(Acts::eBoundTime, Acts::eBoundTime));
      }

      std::optional<ConstTrackProxy> trackOpt = findTrack(tracks, trk);
      if (!trackOpt.has_value()) {
        ACTS_DEBUG("Track has no matching input track.");
        diag.tracks.push_back(td);
        continue;
      }
      const ConstTrackProxy& inputTrk = *trackOpt;
      const SimParticle* particle =
          findParticle(particles, trackParticleMatching, inputTrk, logger());
      if (particle == nullptr) {
        ACTS_VERBOSE("Track has no matching truth particle.");
      } else {
        td.hasParticle = true;
        td.particleId = particle->particleId();
        td.pdg = static_cast<int>(particle->pdg());
        const SimVertexBarcode pVtx = vtxKey(particle->particleId());
        td.truthVtx = pVtx;
        ++diag.nWithParticle;

        // Optional time-consistency cut: compare the track's own (pre-fit)
        // time to the time of the truth vertex its particle belongs to.
        if (applyTimeCut) {
          if (auto tit = truthVtxTime.find(pVtx); tit != truthVtxTime.end()) {
            td.haveDt = true;
            td.dt = td.oPos4[3] - tit->second;
            td.timeCompatible = std::abs(td.dt) <= m_cfg.timeWindow;
          }
          // if the truth-vertex time is unknown, keep timeCompatible = true
        }

        if (td.timeCompatible) {
          contributingTruthVertices.emplace_back(pVtx, trk.trackWeight);
        } else {
          ++diag.nTimeRejected;
        }
      }
      diag.tracks.push_back(td);
    }

    // Find true vertex that contributes most to the reconstructed vertex
    std::map<SimVertexBarcode, std::pair<int, double>> fmap;
    for (const auto& [vtxId, weight] : contributingTruthVertices) {
      ++fmap[vtxId].first;
      fmap[vtxId].second += weight;
    }
    double truthMajorityVertexTrackWeights = 0;
    std::optional<SimVertexBarcode> truthMajorityVertexId = std::nullopt;
    for (const auto& [vtxId, counter] : fmap) {
      if (counter.second > truthMajorityVertexTrackWeights) {
        truthMajorityVertexId = vtxId;
        truthMajorityVertexTrackWeights = counter.second;
      }
    }

    // Track-count based fractions (1)/(2)/(3) for the printout. Tracks rejected
    // by the time-window cut are excluded from matched/contamination (they are
    // counted separately as nTimeRejected).
    if (truthMajorityVertexId.has_value()) {
      for (const auto& td : diag.tracks) {
        if (td.hasParticle && td.truthVtx.has_value() && td.timeCompatible) {
          if (td.truthVtx.value() == truthMajorityVertexId.value()) {
            ++diag.nMatchedVtx;
          } else {
            ++diag.nContamination;
          }
        }
      }
    }

    double sumPt2 = calcSumPt2(vtx, m_cfg.minTrkWeight);

    double vertexMatchFraction =
        (totalTrackWeight > 0
             ? truthMajorityVertexTrackWeights / totalTrackWeight
             : 0.);
    RecoVertexClassification recoVertexClassification =
        RecoVertexClassification::Unknown;

    if (vertexMatchFraction >= m_cfg.vertexMatchThreshold) {
      recoVertexClassification = RecoVertexClassification::Clean;
    } else {
      recoVertexClassification = RecoVertexClassification::Merged;
    }

    recoToTruthMatching.push_back({truthMajorityVertexId, totalTrackWeight,
                                   truthMajorityVertexTrackWeights,
                                   vertexMatchFraction,
                                   recoVertexClassification});

    auto& recoToTruth = recoToTruthMatching.back();

    // We have to decide if this reco vertex is a split vertex.
    if (!truthMajorityVertexId.has_value()) {
      // No truth vertex matched to this reconstructed vertex
      ACTS_DEBUG("No truth vertex matched to this reconstructed vertex.");
      diags.push_back(std::move(diag));
      continue;
    }

    if (auto it = truthToRecoMatching.find(truthMajorityVertexId.value());
        it != truthToRecoMatching.end()) {
      // This truth vertex is already matched to a reconstructed vertex so we
      // are dealing with a split vertex.

      // We have to decide which of the two reconstructed vertices is the
      // split vertex.
      if (sumPt2 <= it->second.recoSumPt2) {
        // Since the sumPt2 is smaller we can simply call this a split vertex

        recoToTruth.classification = RecoVertexClassification::Split;

        // Keep the existing truth to reco matching
      } else {
        // The sumPt2 is larger, so we call the other vertex a split vertex.

        auto& otherRecoToTruth = recoToTruthMatching.at(it->second.recoIndex);
        // Swap the classification
        recoToTruth.classification = otherRecoToTruth.classification;
        otherRecoToTruth.classification = RecoVertexClassification::Split;

        // Overwrite the truth to reco matching
        it->second = {vtxIndex, sumPt2};
      }
    } else {
      truthToRecoMatching[truthMajorityVertexId.value()] = {vtxIndex, sumPt2};
    }

    diags.push_back(std::move(diag));
  }

  // ----------------------------------------------------------------------
  //  PRINTOUTS
  // ----------------------------------------------------------------------
  if (m_cfg.printReco || m_cfg.printTruth) {
    std::cout << "\n";
    std::cout
        << "==================================================================="
           "============\n";
    std::cout << "[AliceVertexTruthMatcher] event " << ctx.eventNumber
              << " : " << vertices.size() << " reconstructed vertices";
    if (m_inputTruthVertices.isInitialized()) {
      std::cout << ", " << m_inputTruthVertices(ctx).size()
                << " truth vertices";
    }
    std::cout << "\n";
    std::cout << "  modes: collapseSecondaryVertices="
              << (m_cfg.collapseSecondaryVertices ? "ON" : "OFF")
              << ", timeWindowCut="
              << (applyTimeCut
                      ? ("ON (|dt|<=" + std::to_string(m_cfg.timeWindow) + ")")
                      : "OFF")
              << "\n";
    std::cout
        << "==================================================================="
           "============\n";
  }

  auto classificationName = [](RecoVertexClassification c) -> const char* {
    switch (c) {
      case RecoVertexClassification::Clean:
        return "Clean";
      case RecoVertexClassification::Merged:
        return "Merged";
      case RecoVertexClassification::Split:
        return "Split";
      default:
        return "Unknown";
    }
  };

  if (m_cfg.printReco) {
    std::cout << std::fixed;
    for (const auto& [iv, vtx] : Acts::enumerate(vertices)) {
      const VtxDiag& d = diags[iv];
      const VertexToTruthMatching& m = recoToTruthMatching[iv];

      std::cout << "\n----- RECO vertex #" << iv << " -----\n";
      std::cout << std::setprecision(4);
      std::cout << "  position [mm,mm,mm,ns-equiv]: "
                << "x = " << d.pos4[0] << " +/- " << d.sigma4[0] << ",  "
                << "y = " << d.pos4[1] << " +/- " << d.sigma4[1] << ",  "
                << "z = " << d.pos4[2] << " +/- " << d.sigma4[2] << ",  "
                << "t = " << d.pos4[3] << " +/- " << d.sigma4[3] << "\n";

      // matched truth vertex + fractions
      std::cout << "  classification: " << classificationName(m.classification)
                << "\n";
      if (m.vertexId.has_value()) {
        std::cout << "  matched TRUTH vertex: " << m.vertexId.value() << "\n";
      } else {
        std::cout << "  matched TRUTH vertex: (none)\n";
      }
      const double nSel = static_cast<double>(d.nSelected);
      auto frac = [&](int n) {
        return nSel > 0 ? static_cast<double>(n) / nSel : 0.;
      };
      std::cout << std::setprecision(3);
      std::cout << "  contributing tracks (weight > " << m_cfg.minTrkWeight
                << "): " << d.nSelected << "\n";
      std::cout << "    (1) frac. matched to a truth particle      : "
                << frac(d.nWithParticle) << "  (" << d.nWithParticle << "/"
                << d.nSelected << ")\n";
      std::cout << "    (2) frac. from the matched truth vertex    : "
                << frac(d.nMatchedVtx) << "  (" << d.nMatchedVtx << "/"
                << d.nSelected << ")\n";
      std::cout << "    (3) frac. from OTHER vertices (contam.)    : "
                << frac(d.nContamination) << "  (" << d.nContamination << "/"
                << d.nSelected << ")\n";
      if (applyTimeCut) {
        std::cout << "    (4) frac. time-incompatible (|dt|>"
                  << m_cfg.timeWindow << "): " << frac(d.nTimeRejected) << "  ("
                  << d.nTimeRejected << "/" << d.nSelected
                  << ")  [excluded from (2)/(3) and from the majority vote]\n";
      }
      std::cout << "    weight-based match fraction (output)       : "
                << m.matchFraction << "\n";

      // list contributing tracks
      std::cout << "    tracks: per track two lines -- 'fit' = params at the "
                   "vertex, 'orig' = own beamline perigee\n";
      std::cout << "            fit:  (x,y,z,t) | sigma(d0,z0,t) | pT | weight "
                   "| truth\n";
      int printed = 0;
      for (const auto& td : d.tracks) {
        if (m_cfg.maxListPerVertex >= 0 && printed >= m_cfg.maxListPerVertex) {
          std::cout << "      ... (" << (d.tracks.size() - printed)
                    << " more)\n";
          break;
        }
        std::cout << std::setprecision(4);
        std::cout << "      [" << std::setw(3) << printed << "] fit:  "
                  << "(" << td.pos4[0] << ", " << td.pos4[1] << ", "
                  << td.pos4[2] << ", " << td.pos4[3] << ") | "
                  << "(" << td.sigD0 << ", " << td.sigZ0 << ", " << td.sigT
                  << ") | pT=" << std::setprecision(3) << td.pt
                  << " | w=" << td.weight << " | ";
        if (td.hasParticle) {
          std::cout << "particle " << td.particleId << " pdg=" << td.pdg
                    << " vtx=" << td.truthVtx.value();
          if (applyTimeCut && td.haveDt) {
            std::cout << " dt=" << std::setprecision(3) << td.dt
                      << (td.timeCompatible ? " [t-OK]" : " [t-REJECT]");
          }
        } else {
          std::cout << "no truth particle (fake)";
        }
        std::cout << "\n";
        std::cout << std::setprecision(4);
        std::cout << "            orig: "
                  << "(" << td.oPos4[0] << ", " << td.oPos4[1] << ", "
                  << td.oPos4[2] << ", " << td.oPos4[3] << ") | "
                  << "(" << td.oSigD0 << ", " << td.oSigZ0 << ", " << td.oSigT
                  << ")\n";
        ++printed;
      }
    }
  }

  // ----------------------------------------------------------------------
  //  TRUTH vertices
  // ----------------------------------------------------------------------
  if (m_cfg.printTruth) {
    // group particles by their truth vertex
    std::map<SimVertexBarcode, std::vector<const SimParticle*>> partByVtx;
    for (const auto& p : particles) {
      SimVertexBarcode key{SimBarcode{p.particleId()}};
      partByVtx[key].push_back(&p);
    }

    std::cout << "\n----- TRUTH vertices -----\n";
    std::cout << std::setprecision(4);

    auto printTruthVertex = [&](SimVertexBarcode id,
                                const Acts::Vector4& pos4) {
      std::cout << "  TRUTH vertex " << id << " : "
                << "x = " << pos4[0] << ", y = " << pos4[1]
                << ", z = " << pos4[2] << ", t = " << pos4[3]
                << "   (truth: exact, no uncertainty)";
      const SimVertexBarcode lookupId =
          m_cfg.collapseSecondaryVertices ? vtxKey(id.barcode()) : id;
      if (auto it = truthToRecoMatching.find(lookupId);
          it != truthToRecoMatching.end()) {
        std::cout << "  -> reco vertex #" << it->second.recoIndex;
      } else {
        std::cout << "  -> not reconstructed";
      }
      std::cout << "\n";
      // contributing particles
      auto pit = partByVtx.find(id);
      if (pit != partByVtx.end()) {
        std::cout << "    contributing particles (" << pit->second.size()
                  << "): id | pdg | pT | eta\n";
        int printed = 0;
        for (const SimParticle* p : pit->second) {
          if (m_cfg.maxListPerVertex >= 0 &&
              printed >= m_cfg.maxListPerVertex) {
            std::cout << "      ... (" << (pit->second.size() - printed)
                      << " more)\n";
            break;
          }
          const double theta = p->theta();
          const double eta =
              -std::log(std::tan(std::max(theta, 1e-6) / 2.));
          std::cout << "      " << p->particleId()
                    << " | pdg=" << static_cast<int>(p->pdg())
                    << " | pT=" << std::setprecision(3)
                    << p->transverseMomentum() << " | eta="
                    << std::setprecision(3) << eta << std::setprecision(4)
                    << "\n";
          ++printed;
        }
      }
    };

    if (m_inputTruthVertices.isInitialized()) {
      const SimVertexContainer& truthVertices = m_inputTruthVertices(ctx);
      for (const auto& sv : truthVertices) {
        printTruthVertex(sv.vertexId(), sv.position4);
      }
    } else {
      // No truth-vertex collection: derive the position from the production
      // point of the first particle of each vertex.
      for (const auto& [id, plist] : partByVtx) {
        Acts::Vector4 pos4 = Acts::Vector4::Zero();
        if (!plist.empty()) {
          pos4 = plist.front()->fourPosition();
        }
        printTruthVertex(id, pos4);
      }
    }
    std::cout << std::defaultfloat;
  }

  m_outputVertexTruthMatching(ctx, std::move(recoToTruthMatching));
  m_outputTruthVertexMatching(ctx, std::move(truthToRecoMatching));

  return ProcessCode::SUCCESS;
}

}  // namespace AliceActsTrk
