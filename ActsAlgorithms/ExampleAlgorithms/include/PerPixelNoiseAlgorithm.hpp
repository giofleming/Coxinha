#pragma once

// ----------------------------------------------------------------------------
//  PerPixelNoiseAlgorithm
//
//  Injects fake "in-pixel noise" measurements on selected pixel layers, AFTER
//  digitization and BEFORE space-point making / seeding.
//
//  For every sensitive pixel surface in the configured (volume, layer) set, the
//  number of noise pixels is drawn as Poisson( nPixels * noiseProbPerPixel ),
//  where nPixels is taken from the same digitization-config segmentation used by
//  the digitizer. Each noise pixel is a single activated cell at a random
//  (col,row).
//
//  To let noise pixels participate in clustering (option 1 / re-clustering):
//  on every module that receives noise, the existing measurements' cells (read
//  from the input ClusterContainer) and the new noise cells are fed into the
//  ACTS `ModuleClusters` engine, which runs the *same* connected-component
//  clustering + squashing as the digitizer. A noise cell adjacent to a real
//  cluster is therefore absorbed into it (enlarging it / shifting its centroid),
//  and adjacent noise cells form a single multi-pixel noise cluster. Only the
//  measurements connected (through their cells) to a noise cell are re-built;
//  the rest of the module, and modules that receive no noise this event, are
//  passed through unchanged. This matters on the VD, where a layer is a single
//  surface carrying O(10^4) cells per event and would otherwise be re-clustered
//  in full for a handful of noise pixels.
//
//  Truth provenance is preserved: a re-built cluster inherits the contributing
//  particles / sim hits of the real cells that fell into it; pure-noise clusters
//  get no truth entry (so they behave as genuine fakes downstream).
//
//  All output collections are produced under (by default) the canonical
//  digitization names, so the intended wiring is: rename the digitizer outputs
//  (e.g. with a "_raw" suffix) and let this algorithm produce the canonical ones
//  (no downstream changes needed).
// ----------------------------------------------------------------------------

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/IMultiAxis.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"

#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"

#include "ActsExamples/Digitization/DigitizationConfig.hpp"

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace ActsExamples;

namespace AliceActsTrk
{

class PerPixelNoiseAlgorithm final : public IAlgorithm
{
 public:
  struct Config {
    /// Tracking geometry (to enumerate the sensitive pixel surfaces).
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
    /// Same digitization config as the digitizer (segmentation per surface).
    DigiConfigContainer digitizationConfigs;
    /// Random number service.
    std::shared_ptr<const RandomNumbers> randomNumbers = nullptr;

    /// Fake-hit probability per pixel per event.
    double noiseProbPerPixel = 1e-7;
    /// Activation assigned to a noise cell (units of the cluster activation).
    double noiseActivation = 1.0;

    /// Which (volume, layer) get noise. Key = volume id; value = list of layer
    /// ids (empty list => all layers of that volume). Empty map => no-op.
    std::map<std::uint32_t, std::vector<std::uint32_t>> layersByVolume;

    /// Clustering parameters (mirror the digitizer's settings).
    bool mergeCommonCorner = true;
    double mergeNsigma = 1.0;

    /// Full width of the uniform time window assigned to *pure-noise* pixels, in
    /// ACTS native time units (c*t, i.e. mm; 1 ns = 299.792 mm). Only has any
    /// effect on surfaces whose digitization smears eBoundTime (index 5); on all
    /// other surfaces (no time dimension) this is unused and the algorithm is
    /// identical to before. A pure-noise cluster gets t ~ U(0, W) with the
    /// layer's own time variance (sigma^2 read from the smearing config); a noise
    /// cell absorbed into a real cluster keeps the real hit's time. Noise time is
    /// drawn as U(0, W). W = 0 (default) places all pure-noise hits at t = 0.
    double noiseTimeWindow = 0.0;

    // ---- inputs (rename the digitizer outputs to these) ----
    std::string inputMeasurements = "measurements";
    std::string inputClusters = "clusters";
    std::string inputMeasurementParticlesMap = "measurement_particles_map";
    std::string inputMeasurementSimHitsMap = "measurement_simhits_map";

    // ---- outputs (canonical names by default) ----
    std::string outputMeasurements = "measurements";
    std::string outputMeasurementSubset = "measurement_subset";
    std::string outputClusters = "clusters";
    std::string outputMeasurementParticlesMap = "measurement_particles_map";
    std::string outputMeasurementSimHitsMap = "measurement_simhits_map";
    std::string outputParticleMeasurementsMap = "particle_measurements_map";
    std::string outputSimHitMeasurementsMap = "simhit_measurements_map";
  };

  PerPixelNoiseAlgorithm(Config cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  /// One pixel surface that may receive noise (precomputed at construction).
  struct NoiseSurface {
    Acts::GeometryIdentifier geoId;
    const Acts::Surface* surface = nullptr;
    std::shared_ptr<const Acts::IMultiAxis> segmentation;
    std::vector<Acts::BoundIndices> indices;  // geometric indices (loc0, loc1)
    std::size_t nBins0 = 0;
    std::size_t nBins1 = 0;
    bool hasTime = false;   // true if this surface's digi config smears eBoundTime
    double timeVar = 0.0;   // time measurement variance (sigma^2) for noise, native units
  };
  std::vector<NoiseSurface> m_noiseSurfaces;
  /// Lookup into m_noiseSurfaces (filled once, after the vector is complete).
  std::map<Acts::GeometryIdentifier, const NoiseSurface*> m_noiseSurfaceByGeo;

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this, "InputMeasurements"};
  ReadDataHandle<ClusterContainer> m_inputClusters{this, "InputClusters"};
  ReadDataHandle<MeasurementParticlesMap> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMap"};
  ReadDataHandle<MeasurementSimHitsMap> m_inputMeasurementSimHitsMap{
      this, "InputMeasurementSimHitsMap"};

  WriteDataHandle<MeasurementContainer> m_outputMeasurements{this, "OutputMeasurements"};
  WriteDataHandle<MeasurementSubset> m_outputMeasurementSubset{this, "OutputMeasurementSubset"};
  WriteDataHandle<ClusterContainer> m_outputClusters{this, "OutputClusters"};
  WriteDataHandle<MeasurementParticlesMap> m_outputMeasurementParticlesMap{
      this, "OutputMeasurementParticlesMap"};
  WriteDataHandle<MeasurementSimHitsMap> m_outputMeasurementSimHitsMap{
      this, "OutputMeasurementSimHitsMap"};
  WriteDataHandle<ParticleMeasurementsMap> m_outputParticleMeasurementsMap{
      this, "OutputParticleMeasurementsMap"};
  WriteDataHandle<SimHitMeasurementsMap> m_outputSimHitMeasurementsMap{
      this, "OutputSimHitMeasurementsMap"};
};

}  // namespace AliceActsTrk
