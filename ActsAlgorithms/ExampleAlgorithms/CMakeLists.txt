add_library(ExampleAlgorithms SHARED)

# ----------------------------------------
# Sources
# ----------------------------------------
target_sources(ExampleAlgorithms
    PRIVATE
    src/ExampleAlgorithm.cpp
    src/TrackTruthMatcher.cpp
    src/HitRemoverAlgorithm.cpp
    src/TrackMergerAlgorithm.cpp
    src/LutMaker.cpp
    src/ResidualAccumulator.cpp
    src/AO2DMaker.cpp
    src/PerPixelNoiseAlgorithm.cpp
    src/AliceVertexTruthMatcher.cpp
    src/SimVertexBuilder.cpp
    src/SeedAnalysisAlg.cpp
    src/Alice3GbtsSeedingAlgorithm.cpp
    src/SeedAnalysisPlotTool.cpp
    src/UnbiasedResidualPlotTool.cpp
    src/RootUnbiasedResidualWriter.cpp
    src/TrackStateCreatorFactory.cpp
    src/GriddedTrackStateCreator.cpp
    src/AliceTrackFindingAlgorithm.cpp

        # add more .cpp files here
)

# ----------------------------------------
# Include directories
# ----------------------------------------

find_package(ROOT REQUIRED COMPONENTS Core Graf Hist RIO)

find_library(ACTS_EXAMPLES_FRAMEWORK
    NAMES ActsExamplesFramework
    HINTS ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES lib lib64
    REQUIRED
)

message(STATUS "ACTS_EXAMPLES_FRAMEWORK: ${ACTS_EXAMPLES_FRAMEWORK}")

# Needed for Digitization with Noise:
# requires ACTS_SOURCE_DIR=/path/to/acts, or pass -DACTS_SOURCE_DIR=/path/to/acts
find_path(ACTS_DIGI_INCLUDE ActsExamples/Digitization/ModuleClusters.hpp
    HINTS ${ACTS_SOURCE_DIR} $ENV{ACTS_SOURCE_DIR}
    PATH_SUFFIXES Examples/Algorithms/Digitization/include REQUIRED)
find_library(ACTS_EXAMPLES_DIGI ActsExamplesDigitization HINTS ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES lib lib64 REQUIRED)

message(STATUS "ACTS_EXAMPLES_DIGI: ${ACTS_EXAMPLES_DIGI}")

# AliceTrackFindingAlgorithm reuses the CKF finder function of the ACTS examples
# TrackFindingAlgorithm (ActsExamplesTrackFinding), whose header is not
# installed by ACTS either
find_path(ACTS_TRACKFINDING_INCLUDE ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp
    HINTS ${ACTS_SOURCE_DIR} $ENV{ACTS_SOURCE_DIR}
    PATH_SUFFIXES Examples/Algorithms/TrackFinding/include REQUIRED)
find_library(ACTS_EXAMPLES_TRACKFINDING ActsExamplesTrackFinding HINTS ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES lib lib64 REQUIRED)

message(STATUS "ACTS_EXAMPLES_TRACKFINDING: ${ACTS_EXAMPLES_TRACKFINDING}")

target_include_directories(ExampleAlgorithms 
    PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${ACTS_DIGI_INCLUDE}>
    $<BUILD_INTERFACE:${ACTS_TRACKFINDING_INCLUDE}>)

# ----------------------------------------
# Link against Acts
# ----------------------------------------
# SeedAnalysisAlg uses ActsPlugins::toRoot (HistogramConverter) from the ACTS
# Root plugin. The top-level CMakeLists does a component-less find_package(Acts),
# so the PluginRoot target may not be imported yet - request it explicitly here.
find_package(Acts QUIET COMPONENTS PluginRoot)

target_link_libraries(ExampleAlgorithms
    PUBLIC
    Acts::Core
    Acts::PluginRoot
    ${ACTS_EXAMPLES_FRAMEWORK}
    ${ACTS_EXAMPLES_DIGI}
    ${ACTS_EXAMPLES_TRACKFINDING}
    O2Interface
    ROOT::Core
    ROOT::Graf
    ROOT::Hist
    ROOT::RIO)

# Install the ExampleAlgorithms library
include(GNUInstallDirs)
install(TARGETS ExampleAlgorithms
                LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
                )
