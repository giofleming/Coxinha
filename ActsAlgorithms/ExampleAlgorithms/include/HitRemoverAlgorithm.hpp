#pragma once

#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"

#include "ActsExamples/EventData/SimParticle.hpp"

#include <memory>
#include <string>
#include <vector>
#include <unordered_set>

using namespace ActsExamples;

namespace AliceActsTrk
{

class HitRemoverAlgorithm final : public IAlgorithm
{

 public:
  struct Config {
    std::string inputMeasurements = "";
    std::string inputTracks = "";
    std::string inputMeasurementParticlesMap = "";
    bool sortByOldIndex = true;

    std::string usedIndices = "";
    std::string outputMeasurements = "";
    /// Output collection to map measured hits to contributing particles.
    std::string outputMeasurementParticlesMap = "";
    /// Output collection to map particles to measurements.
    std::string outputParticleMeasurementsMap = "";
    /// Mapping between old and new measurements
    /// It is a vector. The position is the new measurement, the value stored is the old index
    /// Needed when merging tracks collections to restore proper indexing
    std::string outputIndexingMap;
  };

  HitRemoverAlgorithm(Config cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  void computeUsedHits(const ConstTrackContainer& tracks,
                       const MeasurementContainer& measurements,
                       const IndexMultimap<SimBarcode>& inputMeasurementParticlesMap,
                       std::unordered_set<size_t>& usedIndices,
                       MeasurementContainer& filteredMeasurements,
                       IndexMultimap<SimBarcode>& measurementParticlesMap,
                       std::vector<size_t>& outputIndexingMap) const;

  void createReindexedMultimap(const IndexMultimap<SimBarcode>& original,
                               const std::unordered_set<size_t>& usedIndices,
                               IndexMultimap<SimBarcode>& reindexed,
                               bool sortByOldIndex = true) const;

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};
  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  ReadDataHandle<IndexMultimap<SimBarcode>> m_inputMeasurementParticlesMap{this, "measurements_particles_map"};

  WriteDataHandle<std::unordered_set<size_t>> m_usedIndices{this, "UsedIndices"};
  WriteDataHandle<MeasurementContainer> m_outputMeasurements{this, "OutputMeasurments"};

  WriteDataHandle<IndexMultimap<SimBarcode>> m_outputMeasurementParticlesMap{
    this, "OutputMeasurementParticlesMap"};

  WriteDataHandle<InverseMultimap<SimBarcode>> m_outputParticleMeasurementsMap{
    this, "OutputParticleMeasurementsMap"};

  WriteDataHandle<std::vector<size_t>> m_outputIndexingMap{
    this, "OutputIndexingMap"};

}; // HitRemoverAlgorithm

} // namespace AliceActsTrk
