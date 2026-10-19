#pragma once

#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <span>

namespace AliceActsTrk
{

/// Optional interface of a per-event CKF track state creator that wants to
/// know the seed the CKF is following.
///
/// AliceTrackFindingAlgorithm checks whether the creator of the event
/// implements it and, if so, passes the measurements of every seed before
/// running the CKF on it. The ACTS TrackFindingAlgorithm does not, so a creator
/// relying on it must only be used with AliceTrackFindingAlgorithm.
class SeedAwareTrackStateCreator
{
 public:
  virtual ~SeedAwareTrackStateCreator() = default;

  /// Set the measurements of the seed about to be followed, empty if there is
  /// no seed. The span is only valid during the call.
  virtual void setSeedMeasurements(
    std::span<const ActsExamples::IndexSourceLink> seedMeasurements) = 0;
};

} // namespace AliceActsTrk
