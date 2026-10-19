#pragma once

#include "AliceTrackFindingAlgorithm.hpp"
#include "GriddedTrackStateCreator.hpp"
#include "TrackStateCreatorFactory.hpp"

#include "ActsPython/Utilities/Macros.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace AliceActsPython
{

// actso2 track state creator interface (base of the factories below) and the
// default factory, which tests every measurement on a surface
inline void addTrackStateCreatorFactory(pybind11::module& mex)
{
  py::class_<AliceActsTrk::TrackStateCreatorFactory,
             std::shared_ptr<AliceActsTrk::TrackStateCreatorFactory>>(
    mex, "TrackStateCreatorFactory");
  mex.def("makeDefaultTrackStateCreatorFactory",
          &AliceActsTrk::makeDefaultTrackStateCreatorFactory);
}

// CKF track state creator binning the measurements of crowded surfaces, to be
// set as AliceTrackFindingAlgorithm.Config.trackStateCreatorFactory
inline void addGriddedTrackStateCreatorFactory(pybind11::module& mex)
{
  using Factory = AliceActsTrk::GriddedTrackStateCreatorFactory;

  auto f = py::class_<Factory,
                      AliceActsTrk::TrackStateCreatorFactory,
                      std::shared_ptr<Factory>>(mex, "GriddedTrackStateCreatorFactory")
             .def(py::init<const Factory::Config&, Acts::Logging::Level>(),
                  "config"_a, "level"_a = Acts::Logging::INFO)
             .def_property_readonly("config", &Factory::config)
             .def("statistics", &Factory::statistics);

  auto c = py::class_<Factory::Config>(f, "Config").def(py::init<>());
  ACTS_PYTHON_STRUCT(c, minMeasurements, volumeIds, measurementsPerCell,
                     windowTolerance, forceSeedMeasurements);

  py::class_<Factory::Statistics>(f, "Statistics")
    .def_readonly("nCalls", &Factory::Statistics::nCalls)
    .def_readonly("nMeasurementsOnSurface",
                  &Factory::Statistics::nMeasurementsOnSurface)
    .def_readonly("nCandidates", &Factory::Statistics::nCandidates)
    .def_readonly("nGrids", &Factory::Statistics::nGrids)
    .def_readonly("nSeeds", &Factory::Statistics::nSeeds)
    .def_readonly("nForced", &Factory::Statistics::nForced)
    .def_readonly("nForcedAboveChi2Measurement",
                  &Factory::Statistics::nForcedAboveChi2Measurement)
    .def_readonly("nForcedAboveChi2Outlier",
                  &Factory::Statistics::nForcedAboveChi2Outlier);
}

// Copy of the ACTS TrackFindingAlgorithm with ALICE 3 branch stopper options
inline void addAliceTrackFindingAlgorithm(pybind11::module& mex)
{
  ACTS_PYTHON_DECLARE_ALGORITHM(
    AliceActsTrk::AliceTrackFindingAlgorithm, mex, "AliceTrackFindingAlgorithm",
    inputMeasurements, inputInitialTrackParameters, inputSeeds, outputTracks,
    trackingGeometry, magneticField, findTracks, findTracksBrem,
    trackStateCreatorFactory, measurementSelectorCfg, trackSelectorCfg,
    maxSteps, twoWay, reverseSearch, seedDeduplication, stayOnSeed,
    computeSharedHits, pixelVolumeIds, stripVolumeIds, maxPixelHoles,
    maxStripHoles, trimTracks, recordMaterialStates, useJosephFormulation,
    constrainToVolumeIds, endOfWorldVolumeIds, secondPassTrimmedCounts);
}

} // namespace AliceActsPython
