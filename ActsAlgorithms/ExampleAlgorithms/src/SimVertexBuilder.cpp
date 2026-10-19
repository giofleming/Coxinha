#include "SimVertexBuilder.hpp"

#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <vector>

namespace AliceActsTrk {

SimVertexBuilder::SimVertexBuilder(const Config& cfg,
                                   std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("SimVertexBuilder", std::move(logger)), m_cfg(cfg) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particle collection");
  }
  if (m_cfg.outputVertices.empty()) {
    throw std::invalid_argument("Missing output vertex collection");
  }
  m_inputParticles.initialize(m_cfg.inputParticles);
  m_outputVertices.initialize(m_cfg.outputVertices);
}

ProcessCode SimVertexBuilder::execute(const AlgorithmContext& ctx) const {
  const auto& particles = m_inputParticles(ctx);

  using PrimaryId = SimBarcode::PrimaryVertexId;
  using SecondaryId = SimBarcode::SecondaryVertexId;

  // Generator vertices (generation 0), keyed by their (already unique) vertex
  // id = (primary, secondary).
  std::map<SimVertexBarcode, SimVertex> genVerts;
  // Highest generator secondary-vertex number seen per primary, so the fresh
  // numbers assigned to Geant4 vertices do not collide with them.
  std::map<PrimaryId, SecondaryId> maxSec;
  // Geant4 transport secondaries (generation > 0), grouped by parent particle
  // = one bucket per physical decay / interaction point.
  std::map<SimBarcode, std::vector<const SimParticle*>> secByParent;

  for (const auto& p : particles) {
    const SimBarcode pid = p.particleId();
    if (pid.generation() == 0) {
      const SimVertexBarcode id{pid};  // = vertexId (primary, secondary, gen 0)
      auto it = genVerts.find(id);
      if (it == genVerts.end()) {
        it = genVerts
                 .emplace(id, SimVertex(id, p.fourPosition(), p.process()))
                 .first;
      }
      it->second.outgoing.insert(pid);
      it->second.incoming.insert(p.parentParticleId());
      auto& mx = maxSec[pid.vertexPrimary()];
      mx = std::max<SecondaryId>(mx, pid.vertexSecondary());
    } else {
      secByParent[p.parentParticleId()].push_back(&p);
    }
  }

  SimVertexContainer vertices;
  for (auto& [id, v] : genVerts) {
    vertices.insert(std::move(v));
  }

  // Assign each physical Geant4 secondary vertex a fresh secondary number
  // (above the generator ones for the same primary) so it is a distinct entry.
  std::map<PrimaryId, SecondaryId> nextSec;
  for (auto& [parent, daughters] : secByParent) {
    const PrimaryId prim = parent.vertexPrimary();
    SecondaryId& n = nextSec[prim];
    if (n == 0) {
      n = maxSec[prim];  // start just above the generator secondaries
    }
    ++n;
    const SimVertexBarcode id{
        SimBarcode{}.withVertexPrimary(prim).withVertexSecondary(n)};

    const SimParticle& first = *daughters.front();
    SimVertex v(id, first.fourPosition(), first.process());
    for (const auto* d : daughters) {
      v.outgoing.insert(d->particleId());
    }
    v.incoming.insert(parent);
    vertices.insert(std::move(v));
  }

  ACTS_DEBUG("Built " << vertices.size() << " sim vertices ("
                      << genVerts.size() << " generator, " << secByParent.size()
                      << " Geant4) from " << particles.size() << " particles");

  m_outputVertices(ctx, std::move(vertices));
  return ProcessCode::SUCCESS;
}

}  // namespace AliceActsTrk
