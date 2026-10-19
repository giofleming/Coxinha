#include "TrackStateCreatorFactory.hpp"

namespace AliceActsTrk
{

namespace
{

/// Considers all measurements on a surface as candidates (the ACTS default).
class DefaultEventTrackStateCreator final
  : public TrackStateCreatorFactory::EventTrackStateCreator
{
 public:
  explicit DefaultEventTrackStateCreator(
    const TrackStateCreatorFactory::EventInput& input)
  {
    m_slAccessor.container = &input.measurements.orderedIndices();

    m_creator.sourceLinkAccessor
      .template connect<&ActsExamples::IndexSourceLinkAccessor::range>(
        &m_slAccessor);
    m_creator.calibrator = input.calibrator;
    m_creator.measurementSelector = input.measurementSelector;
  }

  // the source link accessor points into this object, so it must not be
  // copied or moved
  DefaultEventTrackStateCreator(const DefaultEventTrackStateCreator&) = delete;
  DefaultEventTrackStateCreator& operator=(
    const DefaultEventTrackStateCreator&) = delete;

  DefaultTrackStateCreator::TrackStatesResult createTrackStates(
    const Acts::GeometryContext& gctx,
    const Acts::CalibrationContext& calibrationContext,
    const Acts::Surface& surface,
    const DefaultTrackStateCreator::BoundState& boundState,
    Acts::TrackIndexType prevTip,
    std::vector<DefaultTrackStateCreator::TrackStateProxy>&
      trackStateCandidates,
    DefaultTrackStateCreator::TrackStateContainerBackend& trajectory,
    const Acts::Logger& logger) const override
  {
    return m_creator.createTrackStates(gctx, calibrationContext, surface,
                                       boundState, prevTip,
                                       trackStateCandidates, trajectory, logger);
  }

 private:
  ActsExamples::IndexSourceLinkAccessor m_slAccessor;
  DefaultTrackStateCreator m_creator;
};

class DefaultTrackStateCreatorFactory final : public TrackStateCreatorFactory
{
 public:
  std::unique_ptr<EventTrackStateCreator> makeEventTrackStateCreator(
    const EventInput& input) const override
  {
    return std::make_unique<DefaultEventTrackStateCreator>(input);
  }
};

} // namespace

std::shared_ptr<TrackStateCreatorFactory> makeDefaultTrackStateCreatorFactory()
{
  return std::make_shared<DefaultTrackStateCreatorFactory>();
}

} // namespace AliceActsTrk
