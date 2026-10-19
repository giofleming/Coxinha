#pragma once

#include "ExampleAlgorithm.hpp"
#include "TrackTruthMatcher.hpp"
#include "HitRemoverAlgorithm.hpp"
#include "TrackMergerAlgorithm.hpp"
#include "LutMaker.hpp"
#include "PerPixelNoiseAlgorithm.hpp"
#include "AliceVertexTruthMatcher.hpp"
#include "SimVertexBuilder.hpp"
#include "AliceVertexGenerators.hpp"
#include "ActsExamples/Utilities/VertexGenerators.hpp"

#include "SeedAnalysisPlotTool.hpp"
#include "SeedAnalysisAlg.hpp"
#include "UnbiasedResidualPlotTool.hpp"
#include "RootUnbiasedResidualWriter.hpp"

#include <algorithm>
#include "ActsPython/Utilities/Macros.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace AliceActsPython
{

// Inline ensures the symbol is available in the .so
inline void addExampleAlgorithm(pybind11::module& mex)
{
  ACTS_PYTHON_DECLARE_ALGORITHM(
    AliceActsTrk::ExampleAlgorithm, mex, "ExampleAlgorithm", testString);
}

inline void addTrackTruthMatcher(pybind11::module& mex)
{
  ACTS_PYTHON_DECLARE_ALGORITHM(
    AliceActsTrk::TrackTruthMatcher, mex, "TrackTruthMatcher",
    inputTracks, inputParticles,
    inputMeasurementParticlesMap, outputTrackParticleMatching,
    outputParticleTrackMatching, matchingRatio, doubleMatching,
    looperProtection, loop_absEta, loop_maxPt, loop_maxParticleHits);
}

inline void addHitRemoverAlgorithm(pybind11::module& mex)
{
  ACTS_PYTHON_DECLARE_ALGORITHM(
    AliceActsTrk::HitRemoverAlgorithm, mex, "HitRemoverAlgorithm",
    inputMeasurements, inputTracks, inputMeasurementParticlesMap,
    sortByOldIndex,
    usedIndices, outputMeasurements,
    outputMeasurementParticlesMap, outputParticleMeasurementsMap,
    outputIndexingMap);
}

inline void addTrackMergerAlgorithm(pybind11::module& mex)
{
  ACTS_PYTHON_DECLARE_ALGORITHM(
    AliceActsTrk::TrackMergerAlgorithm, mex, "TrackMergerAlgorithm",
    inputTrackCollections, outputTrackCollection);
}

inline void addLutMakerAlgorithm(pybind11::module& mex)
{
  ACTS_PYTHON_DECLARE_ALGORITHM(
    AliceActsTrk::LutMaker, mex, "LutMakerAlgorithm",
    inputTracks,
    inputParticles,
    // inputMeasurementParticlesMap,
    inputTrackParticleMatching,
    inputParticleTrackMatching,
    lutTag,
    magneticField,
    storeCovarianceSigma,
    produceOutputTree,
    useResiduals,
    outputDirectory
    // , matchingRatio, doubleMatching,
    // looperProtection, loop_absEta, loop_maxPt, loop_maxParticleHits
  );
}

inline void addPerPixelNoiseAlgorithm(pybind11::module& mex)
{
  ACTS_PYTHON_DECLARE_ALGORITHM(
    AliceActsTrk::PerPixelNoiseAlgorithm, mex, "PerPixelNoiseAlgorithm",
    trackingGeometry, digitizationConfigs, randomNumbers,
    noiseProbPerPixel, noiseActivation, layersByVolume,
    mergeCommonCorner, mergeNsigma, noiseTimeWindow,
    inputMeasurements, inputClusters,
    inputMeasurementParticlesMap, inputMeasurementSimHitsMap,
    outputMeasurements, outputMeasurementSubset, outputClusters,
    outputMeasurementParticlesMap, outputMeasurementSimHitsMap,
    outputParticleMeasurementsMap, outputSimHitMeasurementsMap);
}

inline void addAliceVertexTruthMatcher(pybind11::module& mex)
{
  ACTS_PYTHON_DECLARE_ALGORITHM(
    AliceActsTrk::AliceVertexTruthMatcher, mex, "AliceVertexTruthMatcher",
    inputVertices, inputTracks, inputParticles,
    inputTrackParticleMatching, inputTruthVertices,
    outputVertexTruthMatching, outputTruthVertexMatching,
    vertexMatchThreshold, minTrkWeight,
    collapseSecondaryVertices, useTimeWindow, timeWindow,
    printReco, printTruth, maxListPerVertex);
}

inline void addSimVertexBuilder(pybind11::module& mex)
{
  ACTS_PYTHON_DECLARE_ALGORITHM(
    AliceActsTrk::SimVertexBuilder, mex, "SimVertexBuilder",
    inputParticles, outputVertices);
}

// Vertex position generator with a discrete bunch-crossing time distribution.
// Bound as a subclass of the (already-registered) PrimaryVertexPositionGenerator
// so it can be passed straight to addPythia8(vtxGen=...).
inline void addBunchCrossingVertexGenerator(pybind11::module& mex)
{
  py::class_<AliceActsTrk::BunchCrossingVertexPositionGenerator,
             ActsExamples::PrimaryVertexPositionGenerator,
             std::shared_ptr<AliceActsTrk::BunchCrossingVertexPositionGenerator>>(
    mex, "BunchCrossingVertexGenerator")
    .def(py::init<>())
    .def_readwrite(
      "stddev", &AliceActsTrk::BunchCrossingVertexPositionGenerator::stddev)
    .def_readwrite(
      "mean", &AliceActsTrk::BunchCrossingVertexPositionGenerator::mean)
    .def_readwrite(
      "nBunchCrossings",
      &AliceActsTrk::BunchCrossingVertexPositionGenerator::nBunchCrossings)
    .def_readwrite(
      "bunchSpacing",
      &AliceActsTrk::BunchCrossingVertexPositionGenerator::bunchSpacing);
}

// Seed Analysis Plotting tool and Algorithm
inline void addSeedAnalysisPlotTool(pybind11::module& mex)
{
  py::class_<AliceActsTrk::SeedAnalysisPlotTool::Config>(mex, "SeedAnalysisPlotToolConfig")
    .def(py::init<>())
    .def_readwrite("varBinning", &AliceActsTrk::SeedAnalysisPlotTool::Config::varBinning);
}

inline void addSeedAnalysisAlg(pybind11::module& mex)
{
  ACTS_PYTHON_DECLARE_WRITER(AliceActsTrk::SeedAnalysisAlg, mex, "SeedAnalysisAlg",
                             inputSeeds, inputTracks, inputParticles, inputTrackParticleMatching,
                             inputMeasurementParticlesMap,
                             filePath, seedAnalysisPlotToolConfig);
}

// Seed Analysis Plotting tool and Algorithm
inline void addUnbiasedResidualPlotTool(pybind11::module& mex)
{
  py::class_<AliceActsTrk::UnbiasedResidualPlotTool::Config>(mex, "UnbiasedResidualPlotToolConfig")
  .def(py::init<>())
  .def_readwrite("varBinning", &AliceActsTrk::UnbiasedResidualPlotTool::Config::varBinning);
}

inline void addRootUnbiasedResidualWriter(pybind11::module& mex)
{
  ACTS_PYTHON_DECLARE_WRITER(AliceActsTrk::RootUnbiasedResidualWriter, mex, "RootUnbiasedResidualWriter",
                             inputTracks, filePath, fileMode, plotToolConfig);
  }
} // namespace AliceActsPython

