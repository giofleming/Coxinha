#pragma once

// Copy of ActsExamples::TrackFindingAlgorithm (acts tdr-branch, commit
// 3cf41def8) with ALICE 3 specific options for the CKF branch stopper. With
// the default options it gives the same tracks as the ACTS algorithm.
//
// It reuses the ACTS track finder function types, so
// `acts.examples.TrackFindingAlgorithm.makeTrackFinderFunction` plugs in
// unchanged. The track state creator interface (TrackStateCreatorFactory.hpp)
// is owned by actso2, so no ACTS change is needed: it builds against upstream
// ACTS main.

#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/TrackFinding/TrackSelector.hpp"
#include "Acts/TrackFinding/TrackStateCreator.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Seed.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "TrackStateCreatorFactory.hpp"

#include <atomic>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <tbb/combinable.h>
#pragma GCC diagnostic pop

namespace AliceActsTrk
{

class AliceTrackFindingAlgorithm final : public ActsExamples::IAlgorithm
{
 public:
  using ActsAlgorithm = ActsExamples::TrackFindingAlgorithm;

  using TrackFinderOptions = ActsAlgorithm::TrackFinderOptions;
  using TrackFinderResult = ActsAlgorithm::TrackFinderResult;
  using TrackFinderFunction = ActsAlgorithm::TrackFinderFunction;
  /// CKF extensions the track state creator is connected to
  using Extensions =
      Acts::CombinatorialKalmanFilterExtensions<ActsExamples::TrackContainer>;
  using DefaultTrackStateCreator = AliceActsTrk::DefaultTrackStateCreator;
  using TrackStateCreatorFactory = AliceActsTrk::TrackStateCreatorFactory;

  struct Config {
    /// Input measurements collection.
    std::string inputMeasurements;
    /// Input initial track parameter estimates for for each proto track.
    std::string inputInitialTrackParameters;
    /// Input seeds. These are optional and allow for seed deduplication.
    /// The seeds must match the initial track parameters.
    std::string inputSeeds;
    /// Output find trajectories collection.
    std::string outputTracks;

    /// The tracking geometry that should be used.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// The magnetic field that should be used.
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;

    /// Type erased track finder function.
    std::shared_ptr<TrackFinderFunction> findTracks;
    /// Type erased track finder with brem recovery function.
    std::shared_ptr<TrackFinderFunction> findTracksBrem;
    /// Type erased track state creator factory. ALICE 3: if the per event
    /// creator also implements SeedAwareTrackStateCreator, it is given the
    /// measurements of every seed before the CKF runs on it.
    std::shared_ptr<TrackStateCreatorFactory> trackStateCreatorFactory =
        AliceActsTrk::makeDefaultTrackStateCreatorFactory();
    /// CKF measurement selector config
    Acts::MeasurementSelector::Config measurementSelectorCfg;
    /// Track selector config
    std::optional<std::variant<Acts::TrackSelector::Config,
                               Acts::TrackSelector::EtaBinnedConfig>>
        trackSelectorCfg = std::nullopt;

    /// Maximum number of propagation steps
    unsigned int maxSteps = 100000;
    /// Extrapolation strategy
    Acts::TrackExtrapolationStrategy extrapolationStrategy =
        Acts::TrackExtrapolationStrategy::firstOrLast;
    /// Run finding in two directions
    bool twoWay = true;
    /// Whether to run the finding in seed parameter direction or reverse
    /// direction
    bool reverseSearch = false;
    /// Whether to use seed deduplication
    /// This is only available if `inputSeeds` is set.
    bool seedDeduplication = false;
    /// Whether to stick on the seed measurements during track finding.
    /// This is only available if `inputSeeds` is set.
    bool stayOnSeed = false;
    /// Compute shared hit information
    bool computeSharedHits = false;
    /// Whether to trim the tracks
    bool trimTracks = true;
    /// Whether to record track states on material-only surfaces. Disabling this
    /// is faster but leaves those surfaces out of the track, which a refit with
    /// the `DirectNavigator` needs. Not applied to the bremsstrahlung finder,
    /// which does not support it.
    bool recordMaterialStates = true;

    /// Whether to use the Joseph formulation for the Kalman filter update. This
    /// is typically more stable but also more computationally expensive.
    bool useJosephFormulation = false;

    // Pixel and strip volume ids to be used for maxPixel/StripHoles cuts
    std::vector<std::uint32_t> pixelVolumeIds;
    std::vector<std::uint32_t> stripVolumeIds;

    // additional track selector settings
    std::size_t maxPixelHoles = std::numeric_limits<std::size_t>::max();
    std::size_t maxStripHoles = std::numeric_limits<std::size_t>::max();

    /// The volume ids to constrain the track finding to
    std::vector<std::uint32_t> constrainToVolumeIds;
    /// The volume ids to stop the track finding at
    std::vector<std::uint32_t> endOfWorldVolumeIds;

    // ---- ALICE 3 additions ------------------------------------------------

    /// The branch stopper cuts on the hole and outlier counts of the branch,
    /// i.e. those of the track found so far. In the two-way finding the
    /// second (inward) pass starts from the first pass track, whose counts
    /// include the trailing holes and outliers after its last measurement.
    /// Those are trimmed from the final track, and the final selection does
    /// not count them, but in the ACTS algorithm they use up the hole budget
    /// of the second pass. A first pass stopped by too many holes therefore
    /// starts the second pass already over budget, and the second pass stops
    /// at its first surface.
    ///
    /// If true, the second pass starts from the counts of the first pass track
    /// without its trailing holes and outliers, i.e. what the final track
    /// selection will see. Changes the tracks, hence false by default.
    bool secondPassTrimmedCounts = false;
  };

  /// Constructor of the track finding algorithm
  ///
  /// @param config is the config struct to configure the algorithm
  /// @param level is the logging level
  explicit AliceTrackFindingAlgorithm(
      Config config, std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Framework execute method of the track finding algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algorithm flow
  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  void computeSharedHits(ActsExamples::TrackContainer& tracks,
                         const ActsExamples::MeasurementSubset& measurements) const;

  ActsExamples::ProcessCode finalize() override;

 private:
  Config m_cfg;
  std::optional<Acts::TrackSelector> m_trackSelector;

  ActsExamples::ReadDataHandle<ActsExamples::MeasurementSubset>
      m_inputMeasurements{this, "InputMeasurements"};
  ActsExamples::ReadDataHandle<ActsExamples::TrackParametersContainer>
      m_inputInitialTrackParameters{this, "InputInitialTrackParameters"};
  ActsExamples::ReadDataHandle<ActsExamples::SeedContainer> m_inputSeeds{
      this, "InputSeeds"};

  ActsExamples::WriteDataHandle<ActsExamples::ConstTrackContainer>
      m_outputTracks{this, "OutputTracks"};

  mutable std::atomic<std::size_t> m_nTotalSeeds{0};
  mutable std::atomic<std::size_t> m_nDeduplicatedSeeds{0};
  mutable std::atomic<std::size_t> m_nFailedSeeds{0};
  mutable std::atomic<std::size_t> m_nFailedSmoothing{0};
  mutable std::atomic<std::size_t> m_nFailedExtrapolation{0};
  mutable std::atomic<std::size_t> m_nFoundTracks{0};
  mutable std::atomic<std::size_t> m_nSelectedTracks{0};
  mutable std::atomic<std::size_t> m_nStoppedBranches{0};
  mutable std::atomic<std::size_t> m_nSkippedSecondPass{0};
  /// second passes started
  mutable std::atomic<std::size_t> m_nSecondPasses{0};
  /// second passes whose inherited hole or outlier count was already above
  /// the cut, i.e. that stop at their first surface
  mutable std::atomic<std::size_t> m_nSecondPassesOverBudget{0};

  mutable tbb::combinable<Acts::VectorMultiTrajectory::Statistics>
      m_memoryStatistics{[]() {
        auto mtj = std::make_shared<Acts::VectorMultiTrajectory>();
        return mtj->statistics();
      }};
};

}  // namespace AliceActsTrk
