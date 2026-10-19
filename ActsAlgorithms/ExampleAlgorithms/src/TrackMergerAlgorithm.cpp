#include "TrackMergerAlgorithm.hpp"
#include "Acts/Utilities/Logger.hpp"

using namespace ActsExamples;

namespace AliceActsTrk
{

TrackMergerAlgorithm::TrackMergerAlgorithm(TrackMergerAlgorithm::Config cfg,
                                           std::unique_ptr<const Acts::Logger> _logger)
  : IAlgorithm("TrackMergerAlgorithm", std::move(_logger)), m_cfg(std::move(cfg))
{

  m_inputTrackCollections.reserve(m_cfg.inputTrackCollections.size());

  for (size_t itc = 0; itc < m_cfg.inputTrackCollections.size(); itc++) {

    auto tc = m_cfg.inputTrackCollections[itc];

    ACTS_DEBUG("Adding trackCollection " << tc << " to the merged output");

    m_inputTrackCollections.emplace_back(this, tc);
    m_inputTrackCollections[itc].initialize(tc);
  }

  m_outputTrackCollection.initialize(m_cfg.outputTrackCollection);
}

ProcessCode TrackMergerAlgorithm::execute(const AlgorithmContext& ctx) const
{

  // Read input data
  std::vector<ConstTrackContainer> trackCollections;
  trackCollections.reserve(m_cfg.inputTrackCollections.size());

  for (size_t itc = 0; itc < m_cfg.inputTrackCollections.size(); itc++) {
    trackCollections.push_back(m_inputTrackCollections[itc](ctx));
  }

  ACTS_DEBUG("Retrieved and merging " << trackCollections.size() << " track collections");

  if (trackCollections.empty()) {
    ACTS_WARNING("No input track collections were provided.");
    return ProcessCode::SUCCESS;
  }

  // Mutable Output Collection
  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();

  TrackContainer actsTracksContainer(trackContainer, trackStateContainer);

  // Make sure the containers have the same dynamic content. I assume the track Collection size is not empty.
  // And all tracks have the same dynamic columns

  for (const auto& tracks : trackCollections) {
    actsTracksContainer.ensureDynamicColumns(tracks);
  }

  for (size_t itc = 0; itc < trackCollections.size(); ++itc) {
    ACTS_DEBUG("Track collection " << itc
                                  << " has size "
                                  << trackCollections[itc].size());

    for (const auto& track : trackCollections[itc]) {

      auto actsDestProxy = actsTracksContainer.makeTrack();

      // ACTS v46.4.0 multi-pass tracking:
      //
      // With MeasurementSubset, SourceLink indices already remain in
      // the original MeasurementContainer index space.
      // The track-level properties are copied first, then the track states
      // are copied explicitly without remapping the SourceLinks.
      actsDestProxy.copyFromWithoutStates(track);

      for (const auto& srcTrackState : track.trackStatesReversed()) {
        auto destTrackState =
            actsDestProxy.appendTrackState(srcTrackState.getMask());

        destTrackState.copyFrom(srcTrackState,
                                srcTrackState.getMask(),
                                true);
      }

      actsDestProxy.reverseTrackStates();
    }
  }

  ACTS_DEBUG("Merged outputTrackContainer size " << actsTracksContainer.size());

  auto constTrackStateContainer =
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*trackStateContainer));

  auto constTrackContainer = std::make_shared<Acts::ConstVectorTrackContainer>(
      std::move(*trackContainer));

  ConstTrackContainer constTracks{constTrackContainer,
                                  constTrackStateContainer};

  m_outputTrackCollection(ctx, std::move(constTracks));

  return ProcessCode::SUCCESS;
}

} // namespace AliceActsTrk
