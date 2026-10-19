# ---------------------------------------------------------------------------
#  AliceGeometry
#
#  ALICE 3 tracking-geometry construction. This is deliberately NOT part of
#  ExampleAlgorithms: nothing here derives from ActsExamples::IAlgorithm - the
#  entry points are plain factory functions that return a TrackingGeometry.
#
#  Keeping it separate also confines the TGeo dependency. ALICE3Gen3Geometry is
#  the only code in this project that needs ROOT::Geom and Acts::PluginRoot;
#  while it lived in ExampleAlgorithms those were linked PUBLIC from there, so
#  every consumer of that library pulled them in for one file's benefit.
#
#  Target name is AliceGeometry rather than a bare "Geometry", which is generic
#  enough to risk a target-name collision in a larger build.
# ---------------------------------------------------------------------------

add_library(AliceGeometry SHARED)

# ----------------------------------------
# Sources
# ----------------------------------------
target_sources(AliceGeometry
    PRIVATE
    src/ALICE3Gen3Geometry.cpp
    src/Gen3GeometryConfig.cpp
    # add more geometry builders here
)

# ----------------------------------------
# Dependencies
# ----------------------------------------

# Geom is needed for the TGeo* classes used when parsing o2sim_geometry.root
find_package(ROOT REQUIRED COMPONENTS Core Geom)

# The Gen3 builder needs the ACTS Root/TGeo plugin
find_package(Acts QUIET COMPONENTS PluginRoot)

# nlohmann/json for the geometry-config loader
find_package(nlohmann_json CONFIG QUIET)

# Alignment-aware detector elements (ActsExamples::Aligned<TGeoDetectorElement>)
# so that an AlignmentDecorator can misalign gen3 surfaces, as for gen1.
# AlignedTGeoDetectorElement.hpp is installed, but the DetectorCommons headers
# it includes (Aligned.hpp, AlignmentContext.hpp) are not: take them from the
# ACTS source tree, as ExampleAlgorithms does for the Digitization headers.
# requires ACTS_SOURCE_DIR=/path/to/acts, or pass -DACTS_SOURCE_DIR=/path/to/acts
find_path(ACTS_DETECTORCOMMONS_INCLUDE ActsExamples/DetectorCommons/Aligned.hpp
    HINTS ${ACTS_SOURCE_DIR} $ENV{ACTS_SOURCE_DIR}
    PATH_SUFFIXES Examples/Detectors/Common/include REQUIRED)
find_library(ACTS_EXAMPLES_DETECTOR_TGEO ActsExamplesDetectorTGeo HINTS ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES lib lib64 REQUIRED)

message(STATUS "ACTS_EXAMPLES_DETECTOR_TGEO: ${ACTS_EXAMPLES_DETECTOR_TGEO}")

target_include_directories(AliceGeometry
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${ACTS_DETECTORCOMMONS_INCLUDE}>
)

target_link_libraries(AliceGeometry
    PUBLIC
    Acts::Core
    Acts::PluginRoot
    ${ACTS_EXAMPLES_DETECTOR_TGEO}
    ROOT::Core
    ROOT::Geom
    nlohmann_json::nlohmann_json
)

target_compile_features(AliceGeometry PUBLIC cxx_std_20)

# ----------------------------------------
# Install
# ----------------------------------------
include(GNUInstallDirs)
install(TARGETS AliceGeometry
                LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
                )
