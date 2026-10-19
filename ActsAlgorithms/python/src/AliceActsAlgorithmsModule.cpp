#include "ExampleAlgorithmBindings.hpp"
#include "GeometryBindings.hpp"
#include "Alice3GbtsBindings.hpp"
#include "TrackFindingBindings.hpp"
#include <algorithm>

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(AliceActsPythonBindings, mex)
{
  mex.doc() = "Alice Acts Examples";

  // Import acts.examples to register IAlgorithm base class
  py::module_::import("acts.examples");

  // Call all algorithm binding functions here
  AliceActsPython::addExampleAlgorithm(mex);

  AliceActsPython::addTrackTruthMatcher(mex);

  AliceActsPython::addHitRemoverAlgorithm(mex);

  AliceActsPython::addTrackMergerAlgorithm(mex);

  AliceActsPython::addLutMakerAlgorithm(mex);

  AliceActsPython::addPerPixelNoiseAlgorithm(mex);

  AliceActsPython::addAliceVertexTruthMatcher(mex);

  AliceActsPython::addSimVertexBuilder(mex);

  AliceActsPython::addBunchCrossingVertexGenerator(mex);

  // Gen3 (Blueprint) tracking geometry factory
  AliceActsPython::addALICE3Gen3Geometry(mex);

  // GBTS graph-based seeding (ALICE 3 driver, see docs/GBTS.md)
  AliceActsPython::addAlice3GbtsSeedingAlgorithm(mex);

  // Seed analysis
  AliceActsPython::addSeedAnalysisPlotTool(mex);
  AliceActsPython::addSeedAnalysisAlg(mex);

  // Unbiased Residuals - Misalignments
  AliceActsPython::addUnbiasedResidualPlotTool(mex);
  AliceActsPython::addRootUnbiasedResidualWriter(mex);

  // CKF track state creator with binned measurements
  AliceActsPython::addTrackStateCreatorFactory(mex);
  AliceActsPython::addGriddedTrackStateCreatorFactory(mex);
  AliceActsPython::addAliceTrackFindingAlgorithm(mex);

  // Validation
  // AliceActsPython::addResPlotTool(mex);
  // AliceActsPython::addRootTrackFitterPerformanceWriter(mex);
  // AliceActsPython::addRootTrackFinderPerformanceWriter(mex);
}
