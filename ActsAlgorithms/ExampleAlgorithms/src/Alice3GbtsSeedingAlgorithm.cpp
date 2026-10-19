#include "Alice3GbtsSeedingAlgorithm.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Seeding/GbtsNodeStorage.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace AliceActsTrk
{

namespace
{

/// The (r, z) envelope of one detector layer, plus the running sum of the
/// per-surface reference coordinate used for the layer's refCoord.
struct LayerExtent {
  std::uint32_t nSurfaces{};
  double rMin = std::numeric_limits<double>::max();
  double rMax = std::numeric_limits<double>::lowest();
  double zMin = std::numeric_limits<double>::max();
  double zMax = std::numeric_limits<double>::lowest();
  double rSum{};
  double zSum{};

  void add(double rLow, double rHigh, double zLow, double zHigh, double rRef,
           double zRef)
  {
    ++nSurfaces;
    rMin = std::min(rMin, rLow);
    rMax = std::max(rMax, rHigh);
    zMin = std::min(zMin, zLow);
    zMax = std::max(zMax, zHigh);
    rSum += rRef;
    zSum += zRef;
  }
};

/// Measure one sensitive surface, dispatching on its type. The three ALICE 3
/// vertex layers are single cylinders, the FT3 discs may carry a radial disc
/// and everything else is tiled from rectangular chips - the bounds values()
/// mean something different in each case (see docs/GBTS.md section 3).
void measureSurface(const Acts::GeometryContext& gctx,
                    const Acts::Surface& surface, LayerExtent& extent)
{
  const std::vector<double> values = surface.bounds().values();
  const Acts::Vector3 centre = surface.center(gctx);

  if (surface.type() == Acts::Surface::SurfaceType::Cylinder &&
      values.size() >= 2) {
    // CylinderBounds::values() = (R, halfZ, halfPhi, avgPhi, bevelMin, bevelMax)
    const double radius = values[0];
    const double halfZ = values[1];
    extent.add(radius, radius, centre.z() - halfZ, centre.z() + halfZ, radius,
               centre.z());
    return;
  }

  if (surface.type() == Acts::Surface::SurfaceType::Disc && values.size() >= 2) {
    // RadialBounds::values() = (rMin, rMax, halfPhi, avgPhi)
    extent.add(values[0], values[1], centre.z(), centre.z(),
               0.5 * (values[0] + values[1]), centre.z());
    return;
  }

  if (values.size() < 4) {
    // Nothing sensible to measure; fall back to the centre.
    const double r = std::hypot(centre.x(), centre.y());
    extent.add(r, r, centre.z(), centre.z(), r, centre.z());
    return;
  }

  // Planar sensor. RectangleBounds::values() = (-hlX, -hlY, hlX, hlY); use the
  // true corner positions, since the centre alone undersizes a large module.
  const Acts::Vector3 fakeMomentum(1., 1., 1.);
  const std::array<Acts::Vector2, 4> corners{
      Acts::Vector2(values[0], values[1]), Acts::Vector2(values[2], values[1]),
      Acts::Vector2(values[2], values[3]), Acts::Vector2(values[0], values[3])};

  double rLow = std::numeric_limits<double>::max();
  double rHigh = std::numeric_limits<double>::lowest();
  double zLow = std::numeric_limits<double>::max();
  double zHigh = std::numeric_limits<double>::lowest();
  for (const auto& corner : corners) {
    const Acts::Vector3 global = surface.localToGlobal(gctx, corner, fakeMomentum);
    const double r = std::hypot(global.x(), global.y());
    rLow = std::min(rLow, r);
    rHigh = std::max(rHigh, r);
    zLow = std::min(zLow, global.z());
    zHigh = std::max(zHigh, global.z());
  }
  extent.add(rLow, rHigh, zLow, zHigh, std::hypot(centre.x(), centre.y()),
             centre.z());
}

}  // namespace

Alice3GbtsSeedingAlgorithm::Alice3GbtsSeedingAlgorithm(
    Config cfg, std::unique_ptr<const Acts::Logger> inputLogger)
    : IAlgorithm("Alice3GbtsSeedingAlgorithm", std::move(inputLogger)),
      m_cfg(std::move(cfg))
{
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing input space point collection");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing output seed collection");
  }
  if (m_cfg.trackingGeometry == nullptr) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  if (m_cfg.layers.empty()) {
    throw std::invalid_argument("No GBTS layers configured");
  }

  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_outputSeeds.initialize(m_cfg.outputSeeds);

  const auto descriptions = buildLayerDescriptions();
  const auto connections = buildLayerConnections(descriptions);

  if (connections.empty()) {
    throw std::runtime_error(
        "GBTS layer connection table is empty - no pair of the configured "
        "layers is reachable by a track from the configured z0 range and ROI");
  }

  const Acts::Experimental::GbtsZ0Range z0Range{.min = m_cfg.binTableZ0Min,
                                                .max = m_cfg.binTableZ0Max};

  m_geometry = std::make_shared<Acts::Experimental::GbtsGeometry>(
      descriptions, connections, m_cfg.etaBinWidth, z0Range, this->logger());

  // Resolve the dense layer index of every configured layer once, so the
  // per-space-point lookup is a single map probe.
  for (std::size_t position = 0; position < m_cfg.layers.size(); ++position) {
    const auto& spec = m_cfg.layers[position];
    const auto index = m_geometry->layerIndex(gbtsLayerId(position));
    if (!index.has_value()) {
      ACTS_WARNING("GBTS geometry has no layer for volume "
                   << spec.volume << " layer " << spec.layer
                   << "; its space points will be dropped");
      continue;
    }
    m_layerIndices[layerKey(spec.volume, spec.layer)] = *index;
  }

  m_roi.emplace(m_cfg.roiEtaMin, m_cfg.roiEtaMax, m_cfg.roiZMin, m_cfg.roiZMax);
  // The seeder cuts doublets on z0; the ROI is what defines it.
  m_cfg.seedFinderConfig.minZ0 = static_cast<float>(m_roi->zMin());
  m_cfg.seedFinderConfig.maxZ0 = static_cast<float>(m_roi->zMax());

  // Cluster widths and local positions are not available here, so the
  // machine-learning tau narrowing cannot run.
  if (m_cfg.seedFinderConfig.useClusterWidthCuts) {
    ACTS_WARNING(
        "useClusterWidthCuts requested but no cluster widths are provided; "
        "disabling it");
    m_cfg.seedFinderConfig.useClusterWidthCuts = false;
  }

  m_finder.emplace(Acts::Experimental::GraphBasedTrackSeeder::DerivedConfig(
                       m_cfg.seedFinderConfig),
                   m_geometry, this->logger().cloneWithSuffix("Finder"));

  m_filter.emplace(m_cfg.trackingFilterConfig, m_geometry,
                   this->logger().cloneWithSuffix("Filter"));

  ACTS_INFO("GBTS geometry: " << m_geometry->numLayers() << " layers, "
                              << m_geometry->numBins() << " eta bins, "
                              << connections.size() << " layer connections, "
                              << m_geometry->binGroups().size()
                              << " bin groups");
  ACTS_INFO("GBTS cuts: minPt = " << m_cfg.seedFinderConfig.minPt
                                  << ", minSeedLevel = "
                                  << m_cfg.seedFinderConfig.minSeedLevel
                                  << ", z0 in [" << m_cfg.seedFinderConfig.minZ0
                                  << ", " << m_cfg.seedFinderConfig.maxZ0
                                  << "]");
}

std::vector<Acts::Experimental::GbtsLayerDescription>
Alice3GbtsSeedingAlgorithm::buildLayerDescriptions() const
{
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();

  std::map<std::uint64_t, LayerExtent> extents;
  for (const auto& spec : m_cfg.layers) {
    extents[layerKey(spec.volume, spec.layer)] = LayerExtent{};
  }

  m_cfg.trackingGeometry->visitSurfaces(
      [&](const Acts::Surface* surface) {
        if (surface == nullptr || !surface->isSensitive()) {
          return;
        }
        const Acts::GeometryIdentifier geoId = surface->geometryId();
        const auto it = extents.find(
            layerKey(static_cast<std::uint32_t>(geoId.volume()),
                     static_cast<std::uint32_t>(geoId.layer())));
        if (it == extents.end()) {
          return;
        }
        measureSurface(gctx, *surface, it->second);
      },
      /*restrictToSensitives=*/false);

  std::vector<Acts::Experimental::GbtsLayerDescription> descriptions;
  descriptions.reserve(m_cfg.layers.size());

  for (std::size_t position = 0; position < m_cfg.layers.size(); ++position) {
    const auto& spec = m_cfg.layers[position];
    const LayerExtent& extent = extents.at(layerKey(spec.volume, spec.layer));

    if (extent.nSurfaces == 0) {
      throw std::runtime_error(
          "No sensitive surface found for GBTS layer with volume " +
          std::to_string(spec.volume) + " layer " + std::to_string(spec.layer));
    }

    bool barrel{};
    switch (spec.kind) {
      case LayerKind::Barrel:
        barrel = true;
        break;
      case LayerKind::Endcap:
        barrel = false;
        break;
      case LayerKind::Auto:
      default:
        barrel = (extent.zMax - extent.zMin) > (extent.rMax - extent.rMin);
        break;
    }

    const double meanR = extent.rSum / extent.nSurfaces;
    const double meanZ = extent.zSum / extent.nSurfaces;

    Acts::Experimental::GbtsLayerDescription description{};
    description.id = gbtsLayerId(position);
    description.type = barrel ? Acts::Experimental::GbtsLayerType::Barrel
                              : Acts::Experimental::GbtsLayerType::Endcap;
    description.technology = Acts::Experimental::GbtsLayerTechnology::Pixel;
    description.refCoord = static_cast<float>(barrel ? meanR : meanZ);
    description.minBound = static_cast<float>(barrel ? extent.zMin : extent.rMin);
    description.maxBound = static_cast<float>(barrel ? extent.zMax : extent.rMax);
    // barrelOrder stays at -1: GbtsGeometry derives it by sorting the pixel
    // barrel layers on refCoord, which is exactly what we would do here.

    ACTS_DEBUG("GBTS layer " << description.id << " (volume " << spec.volume
                             << " layer " << spec.layer << "): "
                             << (barrel ? "barrel" : "endcap") << ", refCoord "
                             << description.refCoord << ", bounds ["
                             << description.minBound << ", "
                             << description.maxBound << "], " << extent.nSurfaces
                             << " sensitive surfaces");

    descriptions.push_back(description);
  }

  return descriptions;
}

std::vector<Acts::Experimental::GbtsLayerConnection>
Alice3GbtsSeedingAlgorithm::buildLayerConnections(
    const std::vector<Acts::Experimental::GbtsLayerDescription>& descriptions)
    const
{
  std::vector<Acts::Experimental::GbtsLayerConnection> connections;

  // An explicit table wins, for tuning.
  if (!m_cfg.layerConnections.empty()) {
    connections.reserve(m_cfg.layerConnections.size());
    for (const auto& [src, dst] : m_cfg.layerConnections) {
      if (src >= descriptions.size() || dst >= descriptions.size()) {
        throw std::invalid_argument(
            "layerConnections refers to a layer index outside `layers`");
      }
      connections.push_back(Acts::Experimental::GbtsLayerConnection{
          .src = descriptions[src].id, .dst = descriptions[dst].id});
    }
    return connections;
  }

  // Otherwise scan straight lines out of the luminous region: a line from
  // (r = 0, z = z0) with tau = cot(theta) has z(r) = z0 + tau * r. Its
  // crossings, ordered by r, are the layer sequence a track along it would
  // hit; consecutive crossings (and, with maxLayerSkip > 1, crossings a few
  // apart, so a missing layer is tolerated) become connections.
  const std::uint32_t nZ0 = std::max<std::uint32_t>(m_cfg.connectionScanZ0Steps, 1);
  const std::uint32_t nEta =
      std::max<std::uint32_t>(m_cfg.connectionScanEtaSteps, 2);
  const std::uint32_t maxSkip = std::max<std::uint32_t>(m_cfg.maxLayerSkip, 1);

  // Counted rather than just collected: if two layers are seen in both radial
  // orders the table would carry a cycle, and GbtsGeometry's stage-finding
  // ablation loop never terminates on one. The count decides which direction
  // survives; see the guard below.
  std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> found;

  std::vector<std::pair<double, std::uint32_t>> crossings;
  crossings.reserve(descriptions.size());

  for (std::uint32_t iz = 0; iz < nZ0; ++iz) {
    const double z0 =
        nZ0 == 1 ? 0.5 * (m_cfg.binTableZ0Min + m_cfg.binTableZ0Max)
                 : m_cfg.binTableZ0Min + (m_cfg.binTableZ0Max - m_cfg.binTableZ0Min) *
                                             iz / (nZ0 - 1.0);

    for (std::uint32_t ie = 0; ie < nEta; ++ie) {
      const double eta = m_cfg.roiEtaMin + (m_cfg.roiEtaMax - m_cfg.roiEtaMin) *
                                               ie / (nEta - 1.0);
      const double tau = std::sinh(eta);  // cot(theta)

      crossings.clear();
      for (std::uint32_t layer = 0; layer < descriptions.size(); ++layer) {
        const auto& description = descriptions[layer];
        if (description.type == Acts::Experimental::GbtsLayerType::Barrel) {
          const double r = description.refCoord;
          const double z = z0 + tau * r;
          if (z >= description.minBound && z <= description.maxBound) {
            crossings.emplace_back(r, layer);
          }
        } else {
          // endcap: solve z0 + tau * r == refCoord
          if (std::abs(tau) < 1e-9) {
            continue;
          }
          const double r = (description.refCoord - z0) / tau;
          if (r > 0 && r >= description.minBound && r <= description.maxBound) {
            crossings.emplace_back(r, layer);
          }
        }
      }

      std::ranges::sort(crossings);

      for (std::size_t inner = 0; inner + 1 < crossings.size(); ++inner) {
        const std::size_t last =
            std::min(crossings.size(), inner + 1 + std::size_t{maxSkip});
        for (std::size_t outer = inner + 1; outer < last; ++outer) {
          // GBTS builds edges outside in: src is the outer end
          ++found[{crossings[outer].second, crossings[inner].second}];
        }
      }
    }
  }

  // Break any two-layer cycle by keeping the direction the scan saw more
  // often. On the ALICE 3 layer sets this never triggers - the acceptance of
  // a barrel layer and a disc never overlaps in both radial orders - but a
  // future geometry or a wider ROI could change that, and a cycle would hang
  // the seeder rather than fail.
  std::vector<std::pair<std::uint32_t, std::uint32_t>> ordered;
  ordered.reserve(found.size());
  for (const auto& [pair, count] : found) {
    const auto reverse = found.find({pair.second, pair.first});
    if (reverse != found.end()) {
      // keep the more frequent direction; ties go to the lower index pair so
      // the outcome does not depend on map iteration order
      const bool keep = count > reverse->second ||
                        (count == reverse->second && pair.first < pair.second);
      ACTS_WARNING("GBTS layers "
                   << descriptions[pair.first].id << " and "
                   << descriptions[pair.second].id
                   << " are crossed in both radial orders (" << count << " vs "
                   << reverse->second << " samples); keeping only "
                   << (keep ? "this" : "the other") << " direction");
      if (!keep) {
        continue;
      }
    }
    ordered.push_back(pair);
  }

  // Longer cycles cannot come out of a radial ordering, but check rather than
  // risk the non-terminating ablation loop in GbtsGeometry.
  {
    std::map<std::uint32_t, std::vector<std::uint32_t>> outgoing;
    std::map<std::uint32_t, std::uint32_t> inDegree;
    for (const auto& [src, dst] : ordered) {
      outgoing[src].push_back(dst);
      ++inDegree[dst];
      inDegree.try_emplace(src, 0);
    }
    std::vector<std::uint32_t> ready;
    for (const auto& [layer, degree] : inDegree) {
      if (degree == 0) {
        ready.push_back(layer);
      }
    }
    std::size_t visited = 0;
    while (!ready.empty()) {
      const std::uint32_t layer = ready.back();
      ready.pop_back();
      ++visited;
      for (const std::uint32_t next : outgoing[layer]) {
        if (--inDegree[next] == 0) {
          ready.push_back(next);
        }
      }
    }
    if (visited != inDegree.size()) {
      throw std::runtime_error(
          "Generated GBTS layer connection table contains a cycle; GbtsGeometry "
          "would not terminate. Pass an explicit `layerConnections` list.");
    }
  }

  connections.reserve(ordered.size());
  for (const auto& [src, dst] : ordered) {
    ACTS_DEBUG("GBTS connection " << descriptions[src].id << " -> "
                                  << descriptions[dst].id);
    connections.push_back(Acts::Experimental::GbtsLayerConnection{
        .src = descriptions[src].id, .dst = descriptions[dst].id});
  }

  return connections;
}

ProcessCode Alice3GbtsSeedingAlgorithm::execute(
    const AlgorithmContext& ctx) const
{
  const SpacePointContainer& spacePoints = m_inputSpacePoints(ctx);

  Acts::Experimental::GbtsNodeStorage nodeStorage = m_finder->makeNodeStorage();

  std::uint32_t nUnmapped = 0;
  for (const auto& spacePoint : spacePoints) {
    const auto& sourceLinks = spacePoint.sourceLinks();
    if (sourceLinks.empty()) {
      ++nUnmapped;
      continue;
    }
    const auto& sourceLink = sourceLinks.front().get<IndexSourceLink>();
    const Acts::GeometryIdentifier geoId = sourceLink.geometryId();

    const auto it = m_layerIndices.find(
        layerKey(static_cast<std::uint32_t>(geoId.volume()),
                 static_cast<std::uint32_t>(geoId.layer())));
    if (it == m_layerIndices.end()) {
      ++nUnmapped;
      continue;
    }

    nodeStorage.insert(spacePoint.index(), spacePoint.x(), spacePoint.y(),
                       spacePoint.z(), it->second);
  }

  nodeStorage.finalize();

  ACTS_DEBUG("Loaded " << nodeStorage.numberOfNodes() << " graph nodes from "
                       << spacePoints.size() << " space points, " << nUnmapped
                       << " not on a GBTS layer");

  Acts::SeedContainer seeds;
  seeds.assignSpacePointContainer(spacePoints);

  const Acts::Experimental::GraphBasedTrackSeeder::Options options{
      .bFieldInZ = m_cfg.bFieldInZ};

  m_finder->createSeeds(nodeStorage, m_roi.value(), *m_filter, options, seeds);

  ACTS_DEBUG("Created " << seeds.size() << " seeds");

  m_outputSeeds(ctx, std::move(seeds));

  return ProcessCode::SUCCESS;
}

}  // namespace AliceActsTrk
