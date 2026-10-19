#include "HitRemoverAlgorithm.hpp"
#include "Acts/Utilities/Logger.hpp"

#include "ActsExamples/EventData/Index.hpp"

using namespace ActsExamples;

namespace AliceActsTrk
{

HitRemoverAlgorithm::HitRemoverAlgorithm(HitRemoverAlgorithm::Config cfg,
                                         std::unique_ptr<const Acts::Logger> logger)
  : IAlgorithm("HitRemoverAlgorithm", std::move(logger)), m_cfg(std::move(cfg))
{

  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing input hit collection");
  }

  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing input tracks");
  }

  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing output hit collection");
  }

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);

  // Outputs
  m_usedIndices.initialize(m_cfg.usedIndices);
  m_outputMeasurements.initialize(m_cfg.outputMeasurements);

  // Copies of the maps. They should go if we use indices
  m_outputMeasurementParticlesMap.initialize(
    m_cfg.outputMeasurementParticlesMap);
  m_outputParticleMeasurementsMap.initialize(
    m_cfg.outputParticleMeasurementsMap);

  m_outputIndexingMap.initialize(m_cfg.outputIndexingMap);
}

ProcessCode HitRemoverAlgorithm::execute(const AlgorithmContext& ctx) const
{

  // Read input data
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& tracks = m_inputTracks(ctx);
  const auto& inputMeasurementParticlesMap = m_inputMeasurementParticlesMap(ctx);

  // The trimmed measurement<->particles map
  // using IndexMultimap = boost::container::flat_multimap<Index, value_t>;
  IndexMultimap<SimBarcode> measurementParticlesMap;
  measurementParticlesMap.reserve(measurements.size());

  ACTS_DEBUG("Measurements size:: " << measurements.size());
  ACTS_DEBUG("Tracks size:: " << tracks.size());

  // O(1) lookup and shared hits accounted for
  std::unordered_set<size_t> usedIndices;
  MeasurementContainer filteredMeasurements;
  std::vector<size_t> outputIndexingMap;

  computeUsedHits(tracks,
                  measurements,
                  inputMeasurementParticlesMap,
                  usedIndices,
                  filteredMeasurements,
                  measurementParticlesMap,
                  outputIndexingMap);

  ACTS_DEBUG("Used Measurements size:: " << usedIndices.size());
  ACTS_DEBUG("Filtered Measurements size:: " << filteredMeasurements.size());

  m_usedIndices(ctx, std::move(usedIndices));
  m_outputMeasurements(ctx, std::move(filteredMeasurements));

  ACTS_DEBUG("measurementParticlesMap size:: " << measurementParticlesMap.size());

  // invert them before they are moved
  m_outputParticleMeasurementsMap(
    ctx, invertIndexMultimap(measurementParticlesMap));

  m_outputMeasurementParticlesMap(ctx, std::move(measurementParticlesMap));

  m_outputIndexingMap(ctx, std::move(outputIndexingMap));

  return ProcessCode::SUCCESS;
}

void HitRemoverAlgorithm::computeUsedHits(const ConstTrackContainer& tracks,
                                          const MeasurementContainer& measurements,
                                          const IndexMultimap<SimBarcode>& inputMeasurementParticlesMap,
                                          std::unordered_set<size_t>& usedIndices,
                                          MeasurementContainer& filteredMeasurements,
                                          IndexMultimap<SimBarcode>& measurementParticlesMap,
                                          std::vector<size_t>& outputIndexingMap) const
{

  // Compute used hits from all the reconstruced tracks

  for (auto track : tracks) {
    for (auto state : track.trackStatesReversed()) {
      if (!state.typeFlags().isMeasurement()) {
        continue;
      }

      std::size_t hitIndex = state.getUncalibratedSourceLink()
                               .template get<IndexSourceLink>()
                               .index();

      usedIndices.insert(hitIndex);
    }
  } // loop on tracks

  if (measurements.size() > usedIndices.size()) {
    filteredMeasurements.reserve(measurements.size() - usedIndices.size());
    outputIndexingMap.reserve(measurements.size() - usedIndices.size());
  }

  for (size_t i = 0; i < measurements.size(); ++i) {

    if (usedIndices.find(i) == usedIndices.end()) {

      auto meas = measurements.getMeasurement(i);
      filteredMeasurements.copyMeasurement(meas);
      outputIndexingMap.push_back(i);

      ACTS_DEBUG("Found unused measurement with position " << i << " and index" << meas.index());

      /*
      // Now add the measurement to the multi map.
      auto range = inputMeasurementParticlesMap.equal_range(meas.index());
      //ACTS_DEBUG("inputMeasurementParticlesMap ["<< meas.index()<<"] = range("<<*(range.first)<<","<<*(range.second)<<")");

      ACTS_DEBUG("Barcodes for measurement " << meas.index() << ": ");
      for (auto it = range.first; it != range.second; ++it) {
        auto [index, barcode] = *it;
        ACTS_DEBUG("barcode "<< barcode << " ");
      }
      ACTS_DEBUG("\n");

      measurementParticlesMap.insert(range.first,range.second);
      */
    }
  } // loop on measurements

  // Create the reindexed map
  createReindexedMultimap(inputMeasurementParticlesMap,
                          usedIndices,
                          measurementParticlesMap,
                          m_cfg.sortByOldIndex);
}

void HitRemoverAlgorithm::createReindexedMultimap(
  const IndexMultimap<SimBarcode>& original,
  const std::unordered_set<size_t>& usedIndices,
  IndexMultimap<SimBarcode>& reindexed,
  bool sortByOldIndex) const
{

  // Collect all unique indices from original that are NOT in indicesToSkip
  std::unordered_set<Index> selectedIndices;
  for (const auto& [idx, barcode] : original) {
    if (usedIndices.count(idx) == 0) { // Not in skip list
      selectedIndices.insert(idx);
    }
  }

  if (sortByOldIndex) {
    // Convert unordered set to vector

    std::vector<size_t> orderedIndices(selectedIndices.begin(), selectedIndices.end());
    std::sort(orderedIndices.begin(), orderedIndices.end());

    for (size_t newIdx = 0; newIdx < orderedIndices.size(); ++newIdx) {
      size_t oldIdx = orderedIndices[newIdx];

      auto range = original.equal_range(oldIdx);
      for (auto it = range.first; it != range.second; ++it) {
        reindexed.insert({newIdx, it->second});
      }

      ACTS_DEBUG("Barcodes for measurement " << oldIdx << ": " << " new Idx = " << newIdx);
      for (auto it = range.first; it != range.second; ++it) {
        auto [index, barcode] = *it;
        ACTS_DEBUG("barcode " << barcode << " ");
      }
      ACTS_DEBUG("\n");

    } // index loop
  } // sort by index
  else {
    // Direct iteration over set without sorting
    Index newIdx = 0;
    for (const auto& oldIdx : selectedIndices) {
      auto range = original.equal_range(oldIdx);
      for (auto it = range.first; it != range.second; ++it) {
        reindexed.insert({newIdx, it->second});
      }
      ACTS_DEBUG("Barcodes for measurement " << oldIdx << ": " << " new Idx = " << newIdx);
      for (auto it = range.first; it != range.second; ++it) {
        auto [index, barcode] = *it;
        ACTS_DEBUG("barcode " << barcode << " ");
      }
      ACTS_DEBUG("\n");

      newIdx++;
    }
  } // direct insertion
}
} // namespace AliceActsTrk
