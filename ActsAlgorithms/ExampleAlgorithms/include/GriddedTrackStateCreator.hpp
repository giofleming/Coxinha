#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "TrackStateCreatorFactory.hpp"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace AliceActsTrk
{

/// CKF track state creator factory that bins the measurements of crowded
/// surfaces in a (loc0, loc1) grid.
///
/// The default creator calibrates and chi2-tests every measurement on a surface.
/// That is fine for small sensors, but the ALICE 3 vertex detector layers are
/// single full cylinders with O(1000) measurements each in pp with pileup.
/// This creator only hands the measurements inside a window around the
/// prediction to the (unchanged) calibrator and measurement selector.
///
/// The window does not change the selection. For a measurement with residual r
/// and residual covariance S = V + H P H^T,
///     chi2 = r^T S^-1 r >= r_i^2 / S_ii    for every component i.
/// A measurement with r_i^2 > chi2Max * S_ii in loc0 or in loc1 therefore fails
/// both the measurement and the outlier cut, where chi2Max is the largest chi2
/// cut in the measurement selector config, and dropping it changes nothing.
/// The candidates inside the window keep the order of the default creator.
///
/// Assumptions, both true for the ACTS examples TrackFindingAlgorithm:
///  - the calibrator passes the measurement through unchanged
///    (PassThroughCalibrator), since the window uses the measurement
///    parameters and covariance;
///  - stayOnSeed is off. With stayOnSeed the selector keeps the seed
///    measurements if any is a candidate, so removing a far away seed
///    measurement can change which candidates are kept.
///
/// Optionally (`forceSeedMeasurements`) the seed measurements are forced, as
/// ATLAS does with forceTrackOnSeed: on a surface where the current seed has a
/// measurement, that measurement is the only candidate and it is accepted as a
/// measurement without chi2 cut. This changes the tracks. It needs the seed,
/// so it only works with AliceTrackFindingAlgorithm (SeedAwareTrackStateCreator).
class GriddedTrackStateCreatorFactory final
  : public TrackStateCreatorFactory
{
 public:
  struct Config {
    /// Surfaces with at least this many measurements in the event get a grid.
    /// The others use the default creator.
    std::size_t minMeasurements = 64;
    /// If not empty, only surfaces in these volumes get a grid
    std::vector<std::uint32_t> volumeIds{};
    /// Target average number of measurements per grid cell
    double measurementsPerCell = 1.;
    /// Relative margin added to the window, against rounding in the chi2
    double windowTolerance = 1e-6;
    /// Force the seed measurements on the seed surfaces, see above
    bool forceSeedMeasurements = false;
  };

  /// Counters summed over all events and algorithm instances using this
  /// factory
  struct Statistics {
    /// Number of createTrackStates calls on gridded surfaces
    std::size_t nCalls = 0;
    /// Measurements the default creator would have calibrated on those calls
    std::size_t nMeasurementsOnSurface = 0;
    /// Measurements actually passed on to calibration and selection
    std::size_t nCandidates = 0;
    /// Number of (event, surface) grids built
    std::size_t nGrids = 0;
    /// Number of seeds whose measurements were set
    std::size_t nSeeds = 0;
    /// Number of createTrackStates calls on a surface with a forced seed
    /// measurement
    std::size_t nForced = 0;
    /// Forced seed measurements whose predicted chi2 is above the measurement
    /// cut, i.e. that the selector would have made outliers (or dropped)
    std::size_t nForcedAboveChi2Measurement = 0;
    /// Forced seed measurements whose predicted chi2 is above the outlier cut,
    /// i.e. that the selector would have dropped
    std::size_t nForcedAboveChi2Outlier = 0;
  };

  GriddedTrackStateCreatorFactory(
    const Config& cfg, Acts::Logging::Level level = Acts::Logging::INFO);
  ~GriddedTrackStateCreatorFactory() override;

  std::unique_ptr<EventTrackStateCreator> makeEventTrackStateCreator(
    const EventInput& input) const override;

  const Config& config() const { return m_cfg; }

  Statistics statistics() const;

  /// Called by the per event creators when they are destroyed
  void addStatistics(const Statistics& stats) const;

 private:
  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  mutable std::atomic<std::size_t> m_nCalls{0};
  mutable std::atomic<std::size_t> m_nMeasurementsOnSurface{0};
  mutable std::atomic<std::size_t> m_nCandidates{0};
  mutable std::atomic<std::size_t> m_nGrids{0};
  mutable std::atomic<std::size_t> m_nSeeds{0};
  mutable std::atomic<std::size_t> m_nForced{0};
  mutable std::atomic<std::size_t> m_nForcedAboveChi2Measurement{0};
  mutable std::atomic<std::size_t> m_nForcedAboveChi2Outlier{0};
  mutable std::atomic<std::size_t> m_nEvents{0};

  const Acts::Logger& logger() const { return *m_logger; }
};

} // namespace AliceActsTrk
