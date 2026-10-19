#pragma once

// ----------------------------------------------------------------------------
//  SimVertexBuilder
//
//  Builds a SimVertexContainer from a SimParticleContainer, capturing BOTH
//  generator-level and Geant4-level vertices (e.g. weak-decay vertices), which
//  the Geant4 simulation itself does not emit as a vertex collection.
//
//  - Generator particles (generation() == 0) are grouped by their generator
//    vertex id (primary + secondary), which the generator numbers uniquely.
//  - Geant4 transport secondaries (generation() > 0, e.g. K0s/Lambda decay
//    products) are grouped by their PARENT particle, i.e. one vertex per
//    physical decay / interaction point, and are assigned freshly-numbered
//    secondary-vertex ids so that each is a distinct SimVertex (the ACTS
//    barcode scheme numbers Geant4 secondaries by generation, which is not
//    unique per physical vertex).
//
//  The real particle barcodes are kept in each vertex' incoming/outgoing lists,
//  so the link to the particles is preserved. The output can be written with
//  the standard RootVertexWriter.
// ----------------------------------------------------------------------------

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>
#include <string>

using namespace ActsExamples;

namespace AliceActsTrk {

class SimVertexBuilder final : public IAlgorithm {
 public:
  struct Config {
    /// Input simulated-particle collection (must include the Geant4
    /// secondaries), e.g. "particles_simulated".
    std::string inputParticles;
    /// Output vertex collection.
    std::string outputVertices;
  };

  SimVertexBuilder(const Config& cfg,
                   std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  WriteDataHandle<SimVertexContainer> m_outputVertices{this, "OutputVertices"};
};

}  // namespace AliceActsTrk
