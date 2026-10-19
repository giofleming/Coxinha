// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// ---------------------------------------------------------------------------
//  AliceVertexTruthMatcher
//
//  Copy of ActsExamples::VertexTruthMatcher (Examples/Algorithms/TruthTracking)
//  with the SAME matching logic + outputs, but with detailed per-event
//  terminal printouts:
//    1) reconstructed vertices (x,y,z,t + uncertainties), their contributing
//       tracks (x,y,z,t + perigee uncertainties + fit weight), the matched
//       truth vertex, and the fractions of contributing reco tracks that
//         (1) correspond to a truth particle,
//         (2) correspond to a truth particle of the matched truth vertex,
//         (3) correspond to truth particles of OTHER vertices ("contamination").
//    2) truth vertices (x,y,z,t) with their contributing particles.
// ---------------------------------------------------------------------------

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/EventData/Vertex.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <string>

using namespace ActsExamples;

namespace AliceActsTrk
{

/// Matches reconstructed vertices to truth interactions and vice versa,
/// with verbose per-event terminal printouts.
class AliceVertexTruthMatcher final : public IAlgorithm
{
 public:
  struct Config {
    /// Input vertex collection.
    std::string inputVertices;
    /// Tracks object from track finding.
    std::string inputTracks;
    /// All input truth particle collection.
    std::string inputParticles;
    /// Input track-particle matching.
    std::string inputTrackParticleMatching;
    /// Input truth vertex collection (optional, e.g. "vertices_truth").
    /// Used only for the truth-vertex printout. If empty, truth-vertex
    /// positions are derived from the particles' production points.
    std::string inputTruthVertices;
    /// Output vertex-truth matching.
    std::string outputVertexTruthMatching;
    /// Output truth-vertex matching.
    std::string outputTruthVertexMatching;

    /// Minimum fraction of track weight matched between truth
    /// and reco vertices to consider as truth matched.
    double vertexMatchThreshold = 0.7;
    /// Minimum track weight for track to be considered as part of the fit
    double minTrkWeight = 0.1;

    /// Collapse secondary (decay) vertices onto their primary pp interaction
    /// for the matching. When true, a track whose truth particle comes from a
    /// heavy-flavour (or any) decay vertex is counted toward the PRIMARY
    /// vertex of the same pp collision (keying by vertexPrimary() only),
    /// instead of being treated as a separate truth vertex. This is what you
    /// want for primary-vertex purity; leave false to study secondary/HF
    /// vertexing separately. Default OFF -> identical to VertexTruthMatcher.
    bool collapseSecondaryVertices = false;

    /// Optional time-consistency cut for the truth matching. When enabled, a
    /// contributing track is only counted toward its truth particle's vertex
    /// if the track's own (pre-fit) time is within +/- timeWindow of that
    /// truth vertex's time. This removes physically-impossible associations
    /// (e.g. a track measured at 50 ns labelled to a 475 ns vertex) so the
    /// majority truth vertex and the match/contamination fractions become
    /// time-aware. Requires inputTruthVertices to be set (otherwise the cut
    /// is skipped for vertices whose truth time is unknown).
    /// Default OFF -> identical behaviour to ActsExamples::VertexTruthMatcher.
    bool useTimeWindow = false;
    /// Half-width of the time window, in native time units (e.g. 1.0 * u.ns).
    double timeWindow = 0.;

    /// Enable the per-event terminal printouts.
    bool printReco = true;
    bool printTruth = true;
    /// Maximum number of contributing tracks/particles to list per vertex
    /// (-1 = list all).
    int maxListPerVertex = -1;
  };

  explicit AliceVertexTruthMatcher(
      const Config& config,
      std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<VertexContainer> m_inputVertices{this, "InputVertices"};
  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{
      this, "InputTrackParticleMatching"};
  ReadDataHandle<SimVertexContainer> m_inputTruthVertices{
      this, "InputTruthVertices"};

  WriteDataHandle<VertexTruthMatching> m_outputVertexTruthMatching{
      this, "OutputVertexTruthMatching"};
  WriteDataHandle<TruthVertexMatching> m_outputTruthVertexMatching{
      this, "OutputTruthVertexMatching"};
};

} // namespace AliceActsTrk
