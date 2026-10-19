add_library(Validation SHARED)

# ----------------------------------------
# Sources
# ----------------------------------------
target_sources(Validation
    PRIVATE
    src/ResPlotTool.cpp
    src/RootTrackFitterPerformanceWriter.cpp
    # add more .cpp files here
)

# ----------------------------------------
# Include directories
# ----------------------------------------
target_include_directories(Validation
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#        $<INSTALL_INTERFACE:include>
)


find_library(ACTS_EXAMPLES_FRAMEWORK
    NAMES ActsExamplesFramework
    HINTS ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES lib lib64
    REQUIRED
)

#link_directories(${CMAKE_PREFIX_PATH}/lib64)

# ----------------------------------------
# Link against Acts
# ----------------------------------------
target_link_libraries(Validation
    PUBLIC
    Acts::Core
    Acts::PluginRoot
    ${ACTS_EXAMPLES_FRAMEWORK}
  )
