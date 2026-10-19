#pragma once

// CKF track state creator interface of AliceTrackFindingAlgorithm.
//
// Same design as the swappable creator of the ACTS examples
// TrackFindingAlgorithm on the acts tdr-branch (c8c80da4d), but owned by
// actso2, so that the track finding does not need that ACTS change.

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/TrackFinding/TrackStateCreator.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <memory>
#include <vector>

namespace AliceActsTrk
{

/// Track state creator of the default factory. Its calibrator and measurement
/// selector delegate types do not depend on the source link iterator, so they
/// are shared by all creators.
using DefaultTrackStateCreator =
  Acts::TrackStateCreator<ActsExamples::IndexSourceLinkAccessor::Iterator,
                          ActsExamples::TrackContainer>;

/// Factory for the track state creator, i.e. the CKF extension that turns the
/// measurements on a surface into track state candidates for a branch.
///
/// AliceTrackFindingAlgorithm calls it once per event and keeps the returned
/// object alive for the whole track finding of that event. If the returned
/// creator also implements SeedAwareTrackStateCreator, it is given the
/// measurements of every seed before the CKF runs on it.
class TrackStateCreatorFactory
{
 public:
  /// Event data and algorithm pieces a creator is built from
  struct EventInput {
    /// The context of the event being processed
    const ActsExamples::AlgorithmContext& ctx;
    /// The measurements of the event
    const ActsExamples::MeasurementSubset& measurements;
    /// The CKF measurement selector config, e.g. for the chi2 cuts
    const Acts::MeasurementSelector::Config& measurementSelectorCfg;
    /// Calibrator to be used for the measurements
    DefaultTrackStateCreator::Calibrator calibrator;
    /// Measurement selector to be used on the calibrated candidates
    DefaultTrackStateCreator::MeasurementSelector measurementSelector;
  };

  /// Per-event track state creator
  class EventTrackStateCreator
  {
   public:
    virtual ~EventTrackStateCreator() = default;

    /// Create the track states for the measurements on a surface, see
    /// `Acts::CombinatorialKalmanFilterExtensions::createTrackStates`
    virtual DefaultTrackStateCreator::TrackStatesResult createTrackStates(
      const Acts::GeometryContext& gctx,
      const Acts::CalibrationContext& calibrationContext,
      const Acts::Surface& surface,
      const DefaultTrackStateCreator::BoundState& boundState,
      Acts::TrackIndexType prevTip,
      std::vector<DefaultTrackStateCreator::TrackStateProxy>&
        trackStateCandidates,
      DefaultTrackStateCreator::TrackStateContainerBackend& trajectory,
      const Acts::Logger& logger) const = 0;
  };

  virtual ~TrackStateCreatorFactory() = default;

  /// Create the creator for one event
  virtual std::unique_ptr<EventTrackStateCreator> makeEventTrackStateCreator(
    const EventInput& input) const = 0;
};

/// The default factory: every measurement on a surface is a candidate, as in
/// the ACTS examples TrackFindingAlgorithm.
std::shared_ptr<TrackStateCreatorFactory> makeDefaultTrackStateCreatorFactory();

} // namespace AliceActsTrk
