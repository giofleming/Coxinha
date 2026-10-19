#include "GriddedTrackStateCreator.hpp"

#include "SeedAwareTrackStateCreator.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFinding/TrackStateCreator.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace AliceActsTrk
{

namespace
{

using Factory = TrackStateCreatorFactory;
using DefaultCreator = DefaultTrackStateCreator;
using SourceLinkIterator = std::vector<Acts::SourceLink>::const_iterator;
// Same calibrator and measurement selector delegate types as the default
// creator, only the source link iterator differs
using GridCreator =
  Acts::TrackStateCreator<SourceLinkIterator, ActsExamples::TrackContainer>;

/// Largest chi2 cut (measurement or outlier) in the selector config. A
/// candidate above it can neither be selected nor become the outlier.
double maxChi2Cut(const Acts::MeasurementSelector::Config& config)
{
  double maxChi2 = 0.;
  for (std::size_t i = 0; i < config.size(); ++i) {
    const Acts::MeasurementSelectorCuts& cuts = config.valueAt(i);
    // An empty outlier cut means "any best candidate is an outlier", so no
    // measurement can be skipped
    if (cuts.chi2CutOff.empty() || cuts.chi2CutOffOutlier.empty()) {
      return std::numeric_limits<double>::infinity();
    }
    for (double c : cuts.chi2CutOff) {
      maxChi2 = std::max(maxChi2, c);
    }
    for (double c : cuts.chi2CutOffOutlier) {
      maxChi2 = std::max(maxChi2, c);
    }
  }
  return maxChi2;
}

/// Measurement and outlier chi2 cuts of the selector config for a surface and
/// direction, looked up as Acts::MeasurementSelector does. Only used for the
/// forced seed measurement diagnostics.
std::pair<double, double> chi2Cuts(
  const Acts::MeasurementSelector::Config& config,
  Acts::GeometryIdentifier geoId, double theta)
{
  constexpr double kMax = std::numeric_limits<double>::max();
  const auto it = config.find(geoId);
  if (it == config.end()) {
    return {kMax, kMax};
  }
  const Acts::MeasurementSelectorCuts& cuts = *it;
  // Acts::MeasurementSelector::getCutsByTheta: the bin edges are compared
  // with min(theta, pi - theta)
  const double constrainedTheta = std::min(theta, std::numbers::pi - theta);
  std::size_t bin = 0;
  while (bin < cuts.etaBins.size() && !(constrainedTheta < cuts.etaBins[bin])) {
    ++bin;
  }
  auto binOrBackOrMax = [bin](const std::vector<double>& v) {
    return v.empty() ? kMax : (bin < v.size() ? v[bin] : v.back());
  };
  return {binOrBackOrMax(cuts.chi2CutOff), binOrBackOrMax(cuts.chi2CutOffOutlier)};
}

/// Predicted chi2 of a calibrated candidate, as Acts::MeasurementSelector
/// computes it
double predictedChi2(const DefaultCreator::TrackStateProxy& trackState)
{
  constexpr int kMaxSize = Acts::eBoundSize;
  const auto n = static_cast<Eigen::Index>(trackState.calibratedSize());
  const auto calibrated = trackState.effectiveCalibrated();
  const auto calibratedCov = trackState.effectiveCalibratedCovariance();
  const auto predicted = trackState.predicted();
  const auto predictedCov = trackState.predictedCovariance();
  const auto indices = trackState.projectorSubspaceIndices();
  Eigen::Matrix<double, Eigen::Dynamic, 1, 0, kMaxSize, 1> residual(n);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, kMaxSize, kMaxSize>
    covariance(n, n);
  for (Eigen::Index i = 0; i < n; ++i) {
    residual(i) = calibrated(i) - predicted(indices[i]);
    for (Eigen::Index j = 0; j < n; ++j) {
      covariance(i, j) = calibratedCov(i, j) + predictedCov(indices[i], indices[j]);
    }
  }
  return residual.dot(covariance.inverse() * residual);
}

/// The measurements of one surface, binned in (loc0, loc1)
struct SurfaceGrid {
  struct Entry {
    double loc0 = 0.;
    double loc1 = 0.;
    double var0 = 0.;
    double var1 = 0.;
    /// position in the default creator's candidate list for this surface
    std::uint32_t rank = 0;
    ActsExamples::Index index = 0;
  };

  Acts::GeometryIdentifier geoId;
  std::size_t nMeasurements = 0;
  double min0 = 0.;
  double min1 = 0.;
  double max0 = 0.;
  double max1 = 0.;
  double invWidth0 = 1.;
  double invWidth1 = 1.;
  std::size_t n0 = 1;
  std::size_t n1 = 1;
  double maxVar0 = 0.;
  double maxVar1 = 0.;
  /// entries of cell c are entries[cellOffsets[c], cellOffsets[c+1])
  std::vector<std::uint32_t> cellOffsets;
  std::vector<Entry> entries;

  std::size_t bin0(double loc0) const
  {
    const double b = std::floor((loc0 - min0) * invWidth0);
    return static_cast<std::size_t>(
      std::clamp(b, 0., static_cast<double>(n0 - 1)));
  }
  std::size_t bin1(double loc1) const
  {
    const double b = std::floor((loc1 - min1) * invWidth1);
    return static_cast<std::size_t>(
      std::clamp(b, 0., static_cast<double>(n1 - 1)));
  }

  /// Fill from the unbinned entries, in default creator order
  void build(std::vector<Entry> unbinned, double measurementsPerCell)
  {
    nMeasurements = unbinned.size();
    min0 = min1 = std::numeric_limits<double>::max();
    max0 = max1 = std::numeric_limits<double>::lowest();
    for (const Entry& e : unbinned) {
      min0 = std::min(min0, e.loc0);
      max0 = std::max(max0, e.loc0);
      min1 = std::min(min1, e.loc1);
      max1 = std::max(max1, e.loc1);
      maxVar0 = std::max(maxVar0, e.var0);
      maxVar1 = std::max(maxVar1, e.var1);
    }

    // cells with roughly the aspect ratio of the occupied area
    const double ext0 = std::max(max0 - min0, 1e-6);
    const double ext1 = std::max(max1 - min1, 1e-6);
    const double nCells = std::clamp(
      static_cast<double>(nMeasurements) / measurementsPerCell, 1., 4e6);
    n0 = static_cast<std::size_t>(
      std::clamp(std::round(std::sqrt(nCells * ext0 / ext1)), 1., nCells));
    n1 = static_cast<std::size_t>(
      std::clamp(std::round(nCells / static_cast<double>(n0)), 1., nCells));
    invWidth0 = static_cast<double>(n0) / ext0;
    invWidth1 = static_cast<double>(n1) / ext1;

    // counting sort into the cells, keeping the default order inside a cell
    std::vector<std::uint32_t> cellOf(unbinned.size());
    cellOffsets.assign(n0 * n1 + 1, 0);
    for (std::size_t i = 0; i < unbinned.size(); ++i) {
      cellOf[i] = static_cast<std::uint32_t>(bin0(unbinned[i].loc0) * n1 +
                                             bin1(unbinned[i].loc1));
      ++cellOffsets[cellOf[i] + 1];
    }
    for (std::size_t c = 0; c < n0 * n1; ++c) {
      cellOffsets[c + 1] += cellOffsets[c];
    }
    entries.resize(unbinned.size());
    std::vector<std::uint32_t> fill(cellOffsets.begin(), cellOffsets.end() - 1);
    for (std::size_t i = 0; i < unbinned.size(); ++i) {
      entries[fill[cellOf[i]]++] = unbinned[i];
    }
  }
};

class GriddedEventTrackStateCreator final
  : public Factory::EventTrackStateCreator,
    public SeedAwareTrackStateCreator
{
 public:
  GriddedEventTrackStateCreator(const Factory::EventInput& input,
                                const GriddedTrackStateCreatorFactory& factory,
                                const Acts::Logger& logger)
    : m_factory(&factory), m_logger(&logger)
  {
    const GriddedTrackStateCreatorFactory::Config& cfg = factory.config();

    m_force = cfg.forceSeedMeasurements;
    m_measurementSelectorCfg = &input.measurementSelectorCfg;
    m_forcedCreator.calibrator = input.calibrator;
    m_forcedCreator.measurementSelector
      .template connect<&GriddedEventTrackStateCreator::selectForced>(this);

    m_slAccessor.container = &input.measurements.orderedIndices();
    m_defaultCreator.sourceLinkAccessor
      .template connect<&ActsExamples::IndexSourceLinkAccessor::range>(
        &m_slAccessor);
    m_defaultCreator.calibrator = input.calibrator;
    m_defaultCreator.measurementSelector = input.measurementSelector;

    m_gridCreator.calibrator = input.calibrator;
    m_gridCreator.measurementSelector = input.measurementSelector;

    const double chi2Max = maxChi2Cut(input.measurementSelectorCfg);
    if (!std::isfinite(chi2Max)) {
      ACTS_WARNING(
        "No finite outlier chi2 cut, every measurement can be selected. "
        "Falling back to the default track state creator.");
      return;
    }
    m_windowChi2 = chi2Max * (1. + cfg.windowTolerance);

    // the ordered indices are sorted by geometry id, in default creator order
    // inside a surface
    const auto& ordered = input.measurements.orderedIndices();
    auto it = ordered.begin();
    while (it != ordered.end()) {
      const Acts::GeometryIdentifier geoId = it->geometryId();
      auto end = std::find_if(it, ordered.end(), [&](const auto& sl) {
        return sl.geometryId() != geoId;
      });
      const auto nOnSurface = static_cast<std::size_t>(std::distance(it, end));
      const bool volumeSelected =
        cfg.volumeIds.empty() ||
        std::ranges::find(cfg.volumeIds, geoId.volume()) !=
          cfg.volumeIds.end();
      if (nOnSurface >= cfg.minMeasurements && volumeSelected) {
        buildGrid(input.measurements, geoId, it, end, cfg, logger);
      }
      it = end;
    }
    m_stats.nGrids = m_grids.size();
  }

  // the source link accessor points into this object
  GriddedEventTrackStateCreator(const GriddedEventTrackStateCreator&) = delete;
  GriddedEventTrackStateCreator& operator=(
    const GriddedEventTrackStateCreator&) = delete;

  ~GriddedEventTrackStateCreator() override
  {
    if (m_force && m_stats.nSeeds == 0 && m_nAllCalls > 0) {
      const Acts::Logger& logger = *m_logger;
      ACTS_ERROR(
        "forceSeedMeasurements is set but no seed was ever given: "
        "it needs AliceTrackFindingAlgorithm with input seeds");
    }
    m_factory->addStatistics(m_stats);
  }

  void setSeedMeasurements(
    std::span<const ActsExamples::IndexSourceLink> seedMeasurements) override
  {
    m_seedMeasurements.assign(seedMeasurements.begin(), seedMeasurements.end());
    if (!m_seedMeasurements.empty()) {
      ++m_stats.nSeeds;
    }
  }

  DefaultCreator::TrackStatesResult createTrackStates(
    const Acts::GeometryContext& gctx,
    const Acts::CalibrationContext& calibrationContext,
    const Acts::Surface& surface, const DefaultCreator::BoundState& boundState,
    Acts::TrackIndexType prevTip,
    std::vector<DefaultCreator::TrackStateProxy>& trackStateCandidates,
    DefaultCreator::TrackStateContainerBackend& trajectory,
    const Acts::Logger& logger) const override
  {
    ++m_nAllCalls;

    if (m_force && !m_seedMeasurements.empty()) {
      // the seed measurements on this surface, if any, are the only candidates
      m_sourceLinks.clear();
      for (const ActsExamples::IndexSourceLink& sl : m_seedMeasurements) {
        if (sl.geometryId() == surface.geometryId()) {
          m_sourceLinks.emplace_back(sl);
        }
      }
      if (!m_sourceLinks.empty()) {
        ++m_stats.nForced;
        return m_forcedCreator.createSourceLinkTrackStates(
          gctx, calibrationContext, surface, boundState, m_sourceLinks.cbegin(),
          m_sourceLinks.cend(), prevTip, trackStateCandidates, trajectory,
          logger);
      }
    }

    const SurfaceGrid* grid = findGrid(surface.geometryId());
    const auto& boundParams = std::get<0>(boundState);
    if (grid == nullptr || !boundParams.covariance().has_value()) {
      return m_defaultCreator.createTrackStates(
        gctx, calibrationContext, surface, boundState, prevTip,
        trackStateCandidates, trajectory, logger);
    }

    const Acts::BoundVector& pred = boundParams.parameters();
    const Acts::BoundMatrix& cov = *boundParams.covariance();
    const double p0 = pred[Acts::eBoundLoc0];
    const double p1 = pred[Acts::eBoundLoc1];
    const double c00 = cov(Acts::eBoundLoc0, Acts::eBoundLoc0);
    const double c11 = cov(Acts::eBoundLoc1, Acts::eBoundLoc1);
    if (!std::isfinite(p0) || !std::isfinite(p1) || !std::isfinite(c00) ||
        !std::isfinite(c11)) {
      // the window cannot decide; the default path applies the selector as is
      // (older ACTS selectors even let a NaN chi2 pass)
      return m_defaultCreator.createTrackStates(
        gctx, calibrationContext, surface, boundState, prevTip,
        trackStateCandidates, trajectory, logger);
    }

    ++m_stats.nCalls;
    m_stats.nMeasurementsOnSurface += grid->nMeasurements;

    // half widths of the window covering every measurement on the surface
    const double w0 = std::sqrt(m_windowChi2 * (c00 + grid->maxVar0));
    const double w1 = std::sqrt(m_windowChi2 * (c11 + grid->maxVar1));

    m_selected.clear();
    if (p0 + w0 >= grid->min0 && p0 - w0 <= grid->max0 &&
        p1 + w1 >= grid->min1 && p1 - w1 <= grid->max1) {
      const std::size_t lo0 = grid->bin0(p0 - w0);
      const std::size_t hi0 = grid->bin0(p0 + w0);
      const std::size_t lo1 = grid->bin1(p1 - w1);
      const std::size_t hi1 = grid->bin1(p1 + w1);
      for (std::size_t i0 = lo0; i0 <= hi0; ++i0) {
        const std::size_t row = i0 * grid->n1;
        for (std::size_t c = grid->cellOffsets[row + lo1];
             c < grid->cellOffsets[row + hi1 + 1]; ++c) {
          const SurfaceGrid::Entry& e = grid->entries[c];
          // per measurement window, with its own variance
          const double r0 = e.loc0 - p0;
          const double r1 = e.loc1 - p1;
          if (r0 * r0 <= m_windowChi2 * (c00 + e.var0) &&
              r1 * r1 <= m_windowChi2 * (c11 + e.var1)) {
            m_selected.emplace_back(e.rank, e.index);
          }
        }
      }
    }

    m_stats.nCandidates += m_selected.size();
    if (m_selected.empty()) {
      // the default creator would find all measurements incompatible and
      // return no track state, i.e. a hole
      return DefaultCreator::TrackStatesResult::success({});
    }

    // restore the default creator order, the selector depends on it for ties
    std::ranges::sort(m_selected);
    m_sourceLinks.clear();
    for (const auto& [rank, index] : m_selected) {
      m_sourceLinks.emplace_back(
        ActsExamples::IndexSourceLink(grid->geoId, index));
    }

    return m_gridCreator.createSourceLinkTrackStates(
      gctx, calibrationContext, surface, boundState, m_sourceLinks.cbegin(),
      m_sourceLinks.cend(), prevTip, trackStateCandidates, trajectory,
      logger);
  }

 private:
  const GriddedTrackStateCreatorFactory* m_factory;
  const Acts::Logger* m_logger;

  ActsExamples::IndexSourceLinkAccessor m_slAccessor;
  DefaultCreator m_defaultCreator;
  GridCreator m_gridCreator;
  /// creator for the forced seed measurements, its selector accepts all
  GridCreator m_forcedCreator;

  bool m_force = false;
  const Acts::MeasurementSelector::Config* m_measurementSelectorCfg = nullptr;
  /// measurements of the seed currently followed
  std::vector<ActsExamples::IndexSourceLink> m_seedMeasurements;
  mutable std::size_t m_nAllCalls = 0;

  /// Measurement selector of the forced seed measurements: all candidates are
  /// accepted as measurements. Records how the chi2 cuts would have treated
  /// them. The Kalman update overwrites the chi2 set here.
  Acts::Result<std::pair<std::vector<DefaultCreator::TrackStateProxy>::iterator,
                         std::vector<DefaultCreator::TrackStateProxy>::iterator>>
    selectForced(std::vector<DefaultCreator::TrackStateProxy>& candidates,
                 bool& isOutlier, const Acts::Logger& /*logger*/) const
  {
    isOutlier = false;
    for (auto& candidate : candidates) {
      const double chi2 = predictedChi2(candidate);
      candidate.chi2() = chi2;
      const auto [cutMeasurement, cutOutlier] =
        chi2Cuts(*m_measurementSelectorCfg,
                 candidate.referenceSurface().geometryId(),
                 candidate.predicted()[Acts::eBoundTheta]);
      if (!(chi2 < cutMeasurement)) {
        ++m_stats.nForcedAboveChi2Measurement;
      }
      if (!(chi2 < cutOutlier)) {
        ++m_stats.nForcedAboveChi2Outlier;
      }
    }
    return std::pair{candidates.begin(), candidates.end()};
  }

  double m_windowChi2 = 0.;
  /// sorted by geometry id
  std::vector<SurfaceGrid> m_grids;

  // per call scratch, the creator is only used by the thread of its event
  mutable std::vector<std::pair<std::uint32_t, ActsExamples::Index>> m_selected;
  mutable std::vector<Acts::SourceLink> m_sourceLinks;
  mutable GriddedTrackStateCreatorFactory::Statistics m_stats;

  const SurfaceGrid* findGrid(Acts::GeometryIdentifier geoId) const
  {
    // a handful of grids, a linear search is fine
    for (const SurfaceGrid& grid : m_grids) {
      if (grid.geoId == geoId) {
        return &grid;
      }
    }
    return nullptr;
  }

  template <typename iterator_t>
  void buildGrid(const ActsExamples::MeasurementSubset& measurements,
                 Acts::GeometryIdentifier geoId, iterator_t begin,
                 iterator_t end,
                 const GriddedTrackStateCreatorFactory::Config& cfg,
                 const Acts::Logger& logger)
  {
    std::vector<SurfaceGrid::Entry> unbinned;
    unbinned.reserve(static_cast<std::size_t>(std::distance(begin, end)));
    std::uint32_t rank = 0;
    for (auto it = begin; it != end; ++it, ++rank) {
      const ActsExamples::Index index = it->index();
      const auto meas = measurements.getMeasurement(index);
      if (!meas.contains(Acts::eBoundLoc0) || !meas.contains(Acts::eBoundLoc1)) {
        ACTS_DEBUG("Surface " << geoId
                              << " has a measurement without loc0 and loc1, "
                                 "using the default creator");
        return;
      }
      const std::size_t i0 = meas.indexOf(Acts::eBoundLoc0);
      const std::size_t i1 = meas.indexOf(Acts::eBoundLoc1);
      SurfaceGrid::Entry& e = unbinned.emplace_back();
      e.loc0 = meas.parameters()[i0];
      e.loc1 = meas.parameters()[i1];
      e.var0 = meas.covariance()(i0, i0);
      e.var1 = meas.covariance()(i1, i1);
      e.rank = rank;
      e.index = index;
    }

    SurfaceGrid& grid = m_grids.emplace_back();
    grid.geoId = geoId;
    grid.build(std::move(unbinned), cfg.measurementsPerCell);
    ACTS_DEBUG("Grid on surface " << geoId << ": " << grid.nMeasurements
                                  << " measurements in " << grid.n0 << " x "
                                  << grid.n1 << " cells");
  }
};

} // namespace

GriddedTrackStateCreatorFactory::GriddedTrackStateCreatorFactory(
  const Config& cfg, Acts::Logging::Level level)
  : m_cfg(cfg),
    m_logger(Acts::getDefaultLogger("GriddedTrackStateCreator", level))
{
  if (m_cfg.measurementsPerCell <= 0.) {
    throw std::invalid_argument("measurementsPerCell must be positive");
  }
  if (m_cfg.windowTolerance < 0.) {
    throw std::invalid_argument("windowTolerance must not be negative");
  }
}

GriddedTrackStateCreatorFactory::~GriddedTrackStateCreatorFactory()
{
  const Statistics stats = statistics();
  if (stats.nCalls > 0) {
    ACTS_INFO("Gridded track state creator: "
              << m_nEvents << " event creators, " << stats.nGrids
              << " surface grids, " << stats.nCalls
              << " calls on gridded surfaces, " << stats.nCandidates << " of "
              << stats.nMeasurementsOnSurface << " measurements calibrated ("
              << 100. * static_cast<double>(stats.nCandidates) /
                   static_cast<double>(stats.nMeasurementsOnSurface)
              << "%)");
  }
  if (stats.nForced > 0) {
    ACTS_INFO("Forced seed measurements: "
              << stats.nSeeds << " seeds, " << stats.nForced
              << " forced surfaces, " << stats.nForcedAboveChi2Measurement
              << " above the measurement chi2 cut, "
              << stats.nForcedAboveChi2Outlier << " above the outlier cut");
  }
}

std::unique_ptr<Factory::EventTrackStateCreator>
  GriddedTrackStateCreatorFactory::makeEventTrackStateCreator(
    const EventInput& input) const
{
  ++m_nEvents;
  return std::make_unique<GriddedEventTrackStateCreator>(input, *this,
                                                         logger());
}

GriddedTrackStateCreatorFactory::Statistics
  GriddedTrackStateCreatorFactory::statistics() const
{
  Statistics stats;
  stats.nCalls = m_nCalls;
  stats.nMeasurementsOnSurface = m_nMeasurementsOnSurface;
  stats.nCandidates = m_nCandidates;
  stats.nGrids = m_nGrids;
  stats.nSeeds = m_nSeeds;
  stats.nForced = m_nForced;
  stats.nForcedAboveChi2Measurement = m_nForcedAboveChi2Measurement;
  stats.nForcedAboveChi2Outlier = m_nForcedAboveChi2Outlier;
  return stats;
}

void GriddedTrackStateCreatorFactory::addStatistics(
  const Statistics& stats) const
{
  m_nCalls += stats.nCalls;
  m_nMeasurementsOnSurface += stats.nMeasurementsOnSurface;
  m_nCandidates += stats.nCandidates;
  m_nGrids += stats.nGrids;
  m_nSeeds += stats.nSeeds;
  m_nForced += stats.nForced;
  m_nForcedAboveChi2Measurement += stats.nForcedAboveChi2Measurement;
  m_nForcedAboveChi2Outlier += stats.nForcedAboveChi2Outlier;
}

} // namespace AliceActsTrk
