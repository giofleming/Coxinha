#include "PerPixelNoiseAlgorithm.hpp"

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "ActsExamples/Digitization/ModuleClusters.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"

#include <array>
#include <cstdint>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <utility>

using namespace ActsExamples;

namespace AliceActsTrk
{

PerPixelNoiseAlgorithm::PerPixelNoiseAlgorithm(PerPixelNoiseAlgorithm::Config cfg,
                                               std::unique_ptr<const Acts::Logger> logger)
  : IAlgorithm("PerPixelNoiseAlgorithm", std::move(logger)), m_cfg(std::move(cfg))
{
  if (m_cfg.inputMeasurements.empty())
    throw std::invalid_argument("Missing input measurements");
  if (m_cfg.inputClusters.empty())
    throw std::invalid_argument("Missing input clusters (needed for re-clustering)");
  if (m_cfg.outputMeasurements.empty())
    throw std::invalid_argument("Missing output measurements");
  if (m_cfg.trackingGeometry == nullptr)
    throw std::invalid_argument("Missing tracking geometry");
  if (m_cfg.randomNumbers == nullptr)
    throw std::invalid_argument("Missing random numbers service");

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputClusters.initialize(m_cfg.inputClusters);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);

  m_outputMeasurements.initialize(m_cfg.outputMeasurements);
  m_outputMeasurementSubset.initialize(m_cfg.outputMeasurementSubset);
  m_outputClusters.initialize(m_cfg.outputClusters);
  m_outputMeasurementParticlesMap.initialize(m_cfg.outputMeasurementParticlesMap);
  m_outputMeasurementSimHitsMap.initialize(m_cfg.outputMeasurementSimHitsMap);
  m_outputParticleMeasurementsMap.initialize(m_cfg.outputParticleMeasurementsMap);
  m_outputSimHitMeasurementsMap.initialize(m_cfg.outputSimHitMeasurementsMap);

  // Precompute the sensitive pixel surfaces that may receive noise: those whose
  // (volume, layer) are selected AND that have a geometric digitization config.
  for (const auto& [geoId, surface] : m_cfg.trackingGeometry->geoIdSurfaceMap()) {
    auto vit = m_cfg.layersByVolume.find(geoId.volume());
    if (vit == m_cfg.layersByVolume.end()) continue;  // volume not selected
    const auto& layers = vit->second;
    if (!layers.empty() &&
        std::find(layers.begin(), layers.end(), geoId.layer()) == layers.end())
      continue;  // layer not selected

    auto dit = m_cfg.digitizationConfigs.find(geoId);
    if (dit == m_cfg.digitizationConfigs.end()) continue;
    const auto& geo = dit->geometricDigiConfig;
    if (geo.indices.empty()) continue;  // not a geometric (pixel) digitizer

    NoiseSurface ns;
    ns.geoId = geoId;
    ns.surface = surface;
    ns.segmentation = geo.segmentation;
    ns.indices = geo.indices;

    if (!geo.segmentation || geo.segmentation->getNAxes() < 2) {
      continue;  // not a valid 2D pixel segmentation, skip this surface
    }
    
    ns.nBins0 = geo.segmentation->getAxis(0).getNBins();
    ns.nBins1 = geo.segmentation->getAxis(1).getNBins();

    // Optional time coordinate: if this surface's digitization *smears*
    // eBoundTime (index 5), noise clusters on it will also carry a time. The
    // measurement variance for a pure-noise time is the layer's own time
    // resolution (sigma^2), recovered from the configured smear function.
    for (const auto& psc : dit->smearingDigiConfig.params) {
      if (psc.index != Acts::eBoundTime) continue;
      double sigma = 0.0;
      if (auto* g = psc.smearFunction.target<const Digitization::Gauss>()) {
        sigma = g->sigma;
      } else if (auto* g = psc.smearFunction.target<const Digitization::GaussTrunc>()) {
        sigma = g->sigma;
      } else if (auto* g = psc.smearFunction.target<const Digitization::GaussClipped>()) {
        sigma = g->sigma;
      } else if (auto* g = psc.smearFunction.target<const Digitization::Exact>()) {
        sigma = g->sigma;
      }
      if (sigma > 0.0) {
        ns.hasTime = true;
        ns.timeVar = sigma * sigma;
      }
      break;
    }

    m_noiseSurfaces.push_back(std::move(ns));
  }
  for (const auto& ns : m_noiseSurfaces) m_noiseSurfaceByGeo[ns.geoId] = &ns;

  // ACTS_INFO("PerPixelNoiseAlgorithm: " << m_noiseSurfaces.size()
            // << " sensitive surfaces eligible for noise, prob/pixel = "
            // << m_cfg.noiseProbPerPixel);
}

ProcessCode PerPixelNoiseAlgorithm::execute(const AlgorithmContext& ctx) const
{
  const auto& measIn = m_inputMeasurements(ctx);
  const auto& clusIn = m_inputClusters(ctx);
  const auto& mpIn = m_inputMeasurementParticlesMap(ctx);
  const auto& msIn = m_inputMeasurementSimHitsMap(ctx);

  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);
  const SimHitIndex kNoiseSource = static_cast<SimHitIndex>(measIn.size());  // sentinel >= all real indices

  // Existing measurements grouped by module (geometry id).
  std::map<Acts::GeometryIdentifier, std::vector<Index>> measByModule;
  for (Index i = 0; i < measIn.size(); ++i)
    measByModule[measIn.getMeasurement(i).geometryId()].push_back(i);

  // Draw the noise cells for each eligible surface this event.
  std::map<Acts::GeometryIdentifier, std::vector<Cluster::Cell>> noiseByModule;
  long nNoiseTotal = 0;
  for (const auto& ns : m_noiseSurfaces) {
    const double mean =
        static_cast<double>(ns.nBins0) * static_cast<double>(ns.nBins1) * m_cfg.noiseProbPerPixel;
    if (mean <= 0.0) continue;
    std::poisson_distribution<int> pois(mean);
    const int nNoise = pois(rng);
    if (nNoise <= 0) continue;
    std::uniform_int_distribution<unsigned int> u0(0, static_cast<unsigned int>(ns.nBins0 - 1));
    std::uniform_int_distribution<unsigned int> u1(0, static_cast<unsigned int>(ns.nBins1 - 1));
    auto& cells = noiseByModule[ns.geoId];
    cells.reserve(nNoise);
    for (int k = 0; k < nNoise; ++k) {
      cells.emplace_back(Cluster::Cell{{u0(rng), u1(rng)},
                                       {Acts::Vector2::Zero(), Acts::Vector2::Zero()},
                                       m_cfg.noiseActivation});
    }
    nNoiseTotal += nNoise;
  }

  // Outputs being built.
  MeasurementContainer measOut;
  ClusterContainer clusOut;
  MeasurementParticlesMap mpOut;
  MeasurementSimHitsMap msOut;

  auto copyProvenance = [&](Index origIdx, Index newIdx) {
    auto pr = mpIn.equal_range(origIdx);
    for (auto it = pr.first; it != pr.second; ++it) mpOut.emplace(newIdx, it->second);
    auto sr = msIn.equal_range(origIdx);
    for (auto it = sr.first; it != sr.second; ++it) msOut.emplace(newIdx, it->second);
  };

  auto passThrough = [&](Index origIdx) {
    auto p = measIn.getMeasurement(origIdx);
    Index newIdx = measOut.copyMeasurement(p).index();
    clusOut.push_back(clusIn.at(origIdx));
    copyProvenance(origIdx, newIdx);
  };

  // Re-cluster a module's existing measurements together with its noise cells.
  auto reCluster = [&](const NoiseSurface& ns, const std::vector<Index>& origIdxs,
                       const std::vector<Cluster::Cell>& noiseCells) {
    ModuleClusters mc(ns.segmentation, ns.indices, /*merge=*/true, m_cfg.mergeNsigma,
                      m_cfg.mergeCommonCorner);

    // Real hits' time (index 5), kept so a re-clustered measurement preserves it
    // (only populated when this surface carries a time dimension).
    std::map<SimHitIndex, std::pair<double, double>> realTime;  // origIdx -> (t, var)

    // Add the existing measurements (their cells), tagged by their original index.
    for (Index origIdx : origIdxs) {
      auto p = measIn.getMeasurement(origIdx);
      const auto fullPar = p.fullParameters();
      const auto fullCov = p.fullCovariance();
      DigitizedParameters dp;
      dp.indices = ns.indices;
      for (auto bi : ns.indices) {
        dp.values.push_back(fullPar[bi]);
        dp.variances.push_back(fullCov(bi, bi));
      }
      dp.cluster = clusIn.at(origIdx);  // carries the real cells
      if (ns.hasTime) {
        realTime[static_cast<SimHitIndex>(origIdx)] = {
            fullPar[Acts::eBoundTime], fullCov(Acts::eBoundTime, Acts::eBoundTime)};
      }
      mc.add(std::move(dp), static_cast<SimHitIndex>(origIdx));
    }

    // Add the noise cells, tagged with the noise sentinel source.
    for (const auto& cell : noiseCells) {
      DigitizedParameters dp;
      dp.indices = ns.indices;
      dp.values.assign(ns.indices.size(), 0.0);
      dp.variances.assign(ns.indices.size(), 0.0);
      dp.cluster.channels.push_back(cell);
      mc.add(std::move(dp), kNoiseSource);
    }

    // Re-run the digitizer's clustering and emit merged measurements.
    for (auto& [dParams, sources] : mc.digitizedParameters()) {
      dParams.cluster.globalPosition =
          measurementGlobalPosition(dParams, *ns.surface, ctx.simGeoContext);
      // Append the time coordinate (index 5) when this surface carries one. A
      // cluster containing any real hit inherits that hit's time; a pure-noise
      // cluster gets a random time in the readout window with the layer's time
      // resolution. Guarded by hasTime so time-less digi is unchanged.
      if (ns.hasTime) {
        double tVal = 0.0, tVar = ns.timeVar;
        bool haveReal = false;
        for (SimHitIndex s : sources) {
          auto it = realTime.find(s);
          if (it != realTime.end()) {
            tVal = it->second.first;
            tVar = it->second.second;
            haveReal = true;
            break;
          }
        }
        if (!haveReal && m_cfg.noiseTimeWindow > 0.0) {
          std::uniform_real_distribution<double> ut(0.0, m_cfg.noiseTimeWindow);
          tVal = ut(rng);
        }
        dParams.indices.push_back(Acts::eBoundTime);
        dParams.values.push_back(tVal);
        dParams.variances.push_back(tVar);
      }
      auto proxy = createMeasurement(measOut, ns.geoId, dParams);
      Index newIdx = proxy.index();
      clusOut.push_back(dParams.cluster);
      for (SimHitIndex s : sources)
        if (s < measIn.size()) copyProvenance(static_cast<Index>(s), newIdx);
      // noise sources (s == kNoiseSource) contribute no truth -> genuine fake.
    }
  };

  // Existing measurements on a module that the noise cells can reach: flood-fill
  // from the noise cells through adjacent cells (same connectivity as the
  // clustering), hopping from a reached measurement to all of its cells. Only
  // these can change when re-clustered; the others would come out identical.
  auto binKey = [](unsigned int b0, unsigned int b1) {
    return (static_cast<std::uint64_t>(b0) << 32) | b1;
  };
  auto connectedToNoise = [&](const std::vector<Index>& idxs,
                              const std::vector<Cluster::Cell>& noiseCells) {
    std::unordered_multimap<std::uint64_t, std::size_t> cellOwner;  // bin -> position in idxs
    for (std::size_t k = 0; k < idxs.size(); ++k)
      for (const auto& cell : clusIn.at(idxs[k]).channels)
        cellOwner.emplace(binKey(cell.bin[0], cell.bin[1]), k);

    std::vector<bool> reached(idxs.size(), false);
    std::vector<std::array<unsigned int, 2>> frontier;
    for (const auto& cell : noiseCells) frontier.push_back(cell.bin);
    while (!frontier.empty()) {
      const auto bin = frontier.back();
      frontier.pop_back();
      for (int d0 = -1; d0 <= 1; ++d0) {
        for (int d1 = -1; d1 <= 1; ++d1) {
          if (!m_cfg.mergeCommonCorner && d0 != 0 && d1 != 0) continue;
          // unsigned wrap-around at the module edge just produces a key no cell has
          auto range = cellOwner.equal_range(binKey(bin[0] + d0, bin[1] + d1));
          for (auto it = range.first; it != range.second; ++it) {
            if (reached[it->second]) continue;
            reached[it->second] = true;
            for (const auto& cell : clusIn.at(idxs[it->second]).channels)
              frontier.push_back(cell.bin);
          }
        }
      }
    }
    return reached;
  };

  // 1) Modules that have existing measurements.
  for (const auto& [geoId, idxs] : measByModule) {
    auto nit = noiseByModule.find(geoId);
    if (nit == noiseByModule.end() || nit->second.empty()) {
      for (Index origIdx : idxs) passThrough(origIdx);  // untouched
      continue;
    }
    const auto reached = connectedToNoise(idxs, nit->second);
    std::vector<Index> toReCluster;
    for (std::size_t k = 0; k < idxs.size(); ++k) {
      if (reached[k]) {
        toReCluster.push_back(idxs[k]);
      } else {
        passThrough(idxs[k]);
      }
    }
    reCluster(*m_noiseSurfaceByGeo.at(geoId), toReCluster, nit->second);
  }

  // 2) Modules that received noise but had no existing measurement (pure noise).
  for (const auto& [geoId, cells] : noiseByModule) {
    if (cells.empty()) continue;
    if (measByModule.find(geoId) != measByModule.end()) continue;  // already handled
    reCluster(*m_noiseSurfaceByGeo.at(geoId), {}, cells);
  }

  ACTS_DEBUG("PerPixelNoiseAlgorithm: injected " << nNoiseTotal << " noise pixels; "
             << measIn.size() << " -> " << measOut.size() << " measurements.");

  // Publish. Build the subset from the stored container (as the digitizer does).
  const auto& stored = m_outputMeasurements(ctx, std::move(measOut));
  std::vector<MeasurementContainer::Index> allIdx(stored.size());
  std::iota(allIdx.begin(), allIdx.end(), Index{0});
  m_outputMeasurementSubset(ctx, MeasurementSubset(stored, std::move(allIdx)));

  m_outputClusters(ctx, std::move(clusOut));

  m_outputParticleMeasurementsMap(ctx, invertIndexMultimap(mpOut));
  m_outputSimHitMeasurementsMap(ctx, invertIndexMultimap(msOut));
  m_outputMeasurementParticlesMap(ctx, std::move(mpOut));
  m_outputMeasurementSimHitsMap(ctx, std::move(msOut));

  return ProcessCode::SUCCESS;
}

}  // namespace AliceActsTrk
