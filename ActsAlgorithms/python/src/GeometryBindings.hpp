#pragma once

// ---------------------------------------------------------------------------
//  Python bindings for the ALICE 3 tracking-geometry builders.
//
//  Kept separate from ExampleAlgorithmBindings.hpp because none of this is an
//  IAlgorithm: the entry points are plain factory functions returning a
//  TrackingGeometry, and they are bound with a normal def() rather than
//  ACTS_PYTHON_DECLARE_ALGORITHM. The implementation lives in the AliceGeometry
//  target (ActsAlgorithms/Geometry).
// ---------------------------------------------------------------------------

#include "ALICE3Gen3Geometry.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>

namespace py = pybind11;
using namespace pybind11::literals;

namespace AliceActsPython
{

/// Gen3 ("Blueprint") tracking geometry for ALICE 3.
inline void addALICE3Gen3Geometry(pybind11::module& mex)
{
  py::enum_<AliceActsTrk::Gen3NavigationPolicy>(mex, "Gen3NavigationPolicy")
      .value("TryAll", AliceActsTrk::Gen3NavigationPolicy::TryAll)
      .value("SurfaceArray", AliceActsTrk::Gen3NavigationPolicy::SurfaceArray);

  mex.def("buildALICE3GeometryGen3",
          &AliceActsTrk::buildALICE3GeometryGen3,
          "rootFile"_a,
          "level"_a = Acts::Logging::INFO,
          "graphvizFile"_a = std::string{},
          "withMaterial"_a = true,
          "configFile"_a = std::string{},
          "navigationPolicy"_a = AliceActsTrk::Gen3NavigationPolicy::TryAll,
          "Build the ALICE 3 Gen3 (Blueprint) tracking geometry from a TGeo "
          "geometry file (o2sim_geometry.root). With withMaterial=True the "
          "layer faces are marked as binned material receivers and passive "
          "cylinders (carbon supports, services, vessel wall) are added, so "
          "that extractMaterialSurfaces() returns them for material mapping. "
          "configFile is the geometry-config JSON; if empty, "
          "gen3_geometry_config.json next to rootFile is used.");
}

} // namespace AliceActsPython
