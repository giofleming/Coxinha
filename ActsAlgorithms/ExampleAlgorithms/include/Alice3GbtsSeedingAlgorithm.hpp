#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsLayerConnection.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"
#include "Acts/Seeding/GbtsRoiDescriptor.hpp"
#include "Acts/Seeding/GbtsTrackingFilter.hpp"
#include "Acts/Seeding/GraphBasedTrackSeeder.hpp"

#include "ActsExamples/EventData/Seed.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

using namespace ActsExamples;

namespace AliceActsTrk
{

/// GBTS (graph-based track seeding) driver for ALICE 3.
///
/// Wraps Acts::Experimental::GraphBasedTrackSeeder from ACTS Core. Unlike
/// ActsExamples::GraphBasedSeedingAlgorithm, which decodes ATLAS/ITk layer ids
/// out of a CSV and a trained connector file, this builds the GBTS geometry
/// straight from the ALICE 3 tracking geometry and generates the layer
/// connection table geometrically. See docs/GBTS.md.
class Alice3GbtsSeedingAlgorithm final : public IAlgorithm
{
 public:
  /// How a layer's GBTS type is decided.
  enum class LayerKind : std::uint8_t {
    /// from the (r, z) envelope: longer in z than deep in r is a barrel
    Auto = 0,
    Barrel = 1,
    Endcap = 2,
  };

  /// One ALICE 3 detector layer GBTS should see, named by its ACTS ids.
  struct LayerSpec {
    /// Acts::GeometryIdentifier volume of the layer
    std::uint32_t volume{};
    /// Acts::GeometryIdentifier layer within that volume
    std::uint32_t layer{};
    /// Barrel/endcap override; Auto derives it from the geometry
    LayerKind kind{LayerKind::Auto};
  };

  struct Config {
    /// input space point container
    std::string inputSpacePoints;
    /// output seed container
    std::string outputSeeds;

    /// ALICE 3 tracking geometry, used to measure the layers
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

    /// the detector layers GBTS builds its graph over, innermost first. The
    /// GBTS layer ids are assigned in this order. Should match the geo
    /// selection the input space points were made with.
    std::vector<LayerSpec> layers;

    /// explicit layer connections as pairs of indices into `layers`
    /// (first = outer/src, second = inner/dst). Empty means: generate them
    /// geometrically, see `maxLayerSkip` and `connectionScan*`.
    std::vector<std::pair<std::uint32_t, std::uint32_t>> layerConnections;

    /// how many layers a generated connection may skip. 1 connects only
    /// consecutive layers along a trajectory, 2 tolerates one missing layer.
    std::uint32_t maxLayerSkip = 2;
    /// number of z0 samples used when generating the connections
    std::uint32_t connectionScanZ0Steps = 21;
    /// number of eta samples used when generating the connections
    std::uint32_t connectionScanEtaSteps = 2001;

    /// eta bin width each GBTS layer is split into
    float etaBinWidth = 0.2f;
    /// z0 range the GBTS eta bin table is built against. Separate from the
    /// seeder's own minZ0/maxZ0 cut, which is taken from the ROI below.
    float binTableZ0Min = -60.0f;
    float binTableZ0Max = 60.0f;

    /// region of interest the seeder searches. There is no phi range; GBTS
    /// always covers the full azimuth.
    float roiEtaMin = -4.0f;
    float roiEtaMax = 4.0f;
    float roiZMin = -60.0f;
    float roiZMax = 60.0f;

    /// magnetic field along z, in ACTS units
    float bFieldInZ = 2.0f * Acts::UnitConstants::T;

    /// graph building and seed extraction
    Acts::Experimental::GraphBasedTrackSeeder::Config seedFinderConfig;
    /// chain-following Kalman filter used during seed extraction
    Acts::Experimental::GbtsTrackingFilter::Config trackingFilterConfig;
  };

  explicit Alice3GbtsSeedingAlgorithm(
      Config cfg,
      std::unique_ptr<const Acts::Logger> inputLogger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

 private:
  /// key into the layer lookup: volume * 1000 + layer
  static std::uint64_t layerKey(std::uint32_t volume, std::uint32_t layer)
  {
    return static_cast<std::uint64_t>(volume) * 1000 + layer;
  }

  /// GBTS layer id we hand out for the n-th entry of Config::layers. Kept
  /// away from zero so an unset id is recognisable in a log line.
  static Acts::Experimental::GbtsExperimentLayerId gbtsLayerId(
      std::size_t position)
  {
    return static_cast<Acts::Experimental::GbtsExperimentLayerId>(1000 +
                                                                  position);
  }

  /// Measure every configured layer against the tracking geometry.
  std::vector<Acts::Experimental::GbtsLayerDescription> buildLayerDescriptions()
      const;

  /// Generate the layer connection table by scanning straight lines from the
  /// luminous region across the configured ROI.
  std::vector<Acts::Experimental::GbtsLayerConnection> buildLayerConnections(
      const std::vector<Acts::Experimental::GbtsLayerDescription>& descriptions)
      const;

  Config m_cfg;

  /// the seeder and the filter, both built once in the constructor
  std::shared_ptr<Acts::Experimental::GbtsGeometry> m_geometry;
  std::optional<Acts::Experimental::GraphBasedTrackSeeder> m_finder;
  std::optional<Acts::Experimental::GbtsTrackingFilter> m_filter;
  std::optional<Acts::Experimental::GbtsRoiDescriptor> m_roi;

  /// (volume, layer) -> dense GBTS layer index, for the per-space-point lookup
  std::map<std::uint64_t, Acts::Experimental::GbtsLayerIndex> m_layerIndices;

  ReadDataHandle<SpacePointContainer> m_inputSpacePoints{this,
                                                        "InputSpacePoints"};
  WriteDataHandle<SeedContainer> m_outputSeeds{this, "OutputSeeds"};
};

}  // namespace AliceActsTrk
