#pragma once

#include "Alice3GbtsSeedingAlgorithm.hpp"

#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace AliceActsPython
{

/// Bind the ACTS GBTS configuration structs and the ALICE 3 driver algorithm.
/// The upstream ACTS bindings expose only `minPt` and `nMaxPhiSlice`, which is
/// not enough to retune the seeder for ALICE 3, so everything is bound here.
inline void addAlice3GbtsSeedingAlgorithm(pybind11::module& mex)
{
  using Seeder = Acts::Experimental::GraphBasedTrackSeeder;
  using Filter = Acts::Experimental::GbtsTrackingFilter;
  using Algorithm = AliceActsTrk::Alice3GbtsSeedingAlgorithm;

  {
    // Acts::Experimental::GraphBasedTrackSeeder::Config is already registered
    // by the upstream ACTS bindings, as `acts.examples.GraphBasedSeedingConfig`
    // - but with only two of its fields. Registering the same C++ type twice
    // throws, so extend the existing class object instead and re-export it
    // here under a name that says what it configures. If a future ACTS stops
    // registering it, bind it from scratch.
    py::object actsExamples = py::module_::import("acts.examples");
    py::object existing =
        py::getattr(actsExamples, "GraphBasedSeedingConfig", py::none());

    const bool reuse = !existing.is_none();
    auto c = reuse ? py::reinterpret_borrow<py::class_<Seeder::Config>>(existing)
                   : py::class_<Seeder::Config>(mex, "GbtsSeedFinderConfig");
    if (!reuse) {
      c.def(py::init<>());
    }

    ACTS_PYTHON_STRUCT(
        c,
        // graph building
        minPt, tuningPt, tripletPtFraction, nMaxPhiSlice, nMaxEdges,
        minDeltaRadius, maxAbsTau, maxOuterRadius, doubletFilterRZ,
        useEtaBinning, minZ0, maxZ0,
        // edge linking
        tauRatioCut, tauRatioPrecut, tauRatioCorr, cutDPhiMax, cutDCurvMax,
        useAdaptiveCuts, validateTriplets, d0Max, matchBeforeCreate,
        matchBeforeCreateMaxEdges, matchBeforeCreateMaxBarrelOrder,
        z0HistogramMaxBarrelOrder, z0Resolution,
        // curvature / phi window shape
        maxCurvatureHighEta, maxCurvatureLowEta, curvatureSplitAbsTau,
        phiWindowSplitDeltaRadius, phiWindowNearOffset, phiWindowNearSlope,
        phiWindowFarOffset, phiWindowFarSlope,
        // seed extraction
        ccaMaxIterations, minSeedLevel, addTriplets, maxAbsEtaAddTriplets,
        edgeMaskMinEta, hitShareThreshold, minSplitSeedSize, maxSplitSeedSize,
        maxSeedSplitEta, maxInvRadDiff,
        // node storage / ML (the cluster width path stays off here)
        useClusterWidthCuts, maxEndcapClusterWidth, moduleHalfLengthY,
        moduleEdgeTolerance, tauLutBinWidth, phiIndexMargin, phiSortBuckets,
        beamSpotCorrection);

    if (reuse) {
      // upstream already patched the kwargs constructor on this class
      mex.attr("GbtsSeedFinderConfig") = existing;
    } else {
      ActsPython::patchKwargsConstructor(c);
    }
  }

  {
    auto c = py::class_<Filter::Config>(mex, "GbtsTrackingFilterConfig")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, sigmaMS, radLen, sigmaX, sigmaY, weightX, weightY,
                       maxDChi2X, maxDChi2Y, addHit, maxCurvature, maxZ0);
    ActsPython::patchKwargsConstructor(c);
  }

  py::enum_<Algorithm::LayerKind>(mex, "GbtsLayerKind")
      .value("Auto", Algorithm::LayerKind::Auto)
      .value("Barrel", Algorithm::LayerKind::Barrel)
      .value("Endcap", Algorithm::LayerKind::Endcap);

  py::class_<Algorithm::LayerSpec>(mex, "GbtsLayerSpec")
      .def(py::init<>())
      .def(py::init([](std::uint32_t volume, std::uint32_t layer,
                       Algorithm::LayerKind kind) {
             return Algorithm::LayerSpec{
                 .volume = volume, .layer = layer, .kind = kind};
           }),
           "volume"_a, "layer"_a, "kind"_a = Algorithm::LayerKind::Auto)
      .def_readwrite("volume", &Algorithm::LayerSpec::volume)
      .def_readwrite("layer", &Algorithm::LayerSpec::layer)
      .def_readwrite("kind", &Algorithm::LayerSpec::kind);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      Algorithm, mex, "Alice3GbtsSeedingAlgorithm", inputSpacePoints,
      outputSeeds, trackingGeometry, layers, layerConnections, maxLayerSkip,
      connectionScanZ0Steps, connectionScanEtaSteps, etaBinWidth, binTableZ0Min,
      binTableZ0Max, roiEtaMin, roiEtaMax, roiZMin, roiZMax, bFieldInZ,
      seedFinderConfig, trackingFilterConfig);
}

}  // namespace AliceActsPython
