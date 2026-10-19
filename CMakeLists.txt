cmake_minimum_required(VERSION 3.16)

project(Alice3ActsTracking
  VERSION 0.1.0
  LANGUAGES CXX
)

#------------------
# GLOBAL settings
#------------------

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

cmake_policy(SET CMP0167 NEW)

# ----------------------------------------
# Find external Acts installation
# ----------------------------------------
# Assumes Acts was installed with CMake config files
# e.g. ActsConfig.cmake
# ActsPluginMille links against the Mille2 target but ActsConfig.cmake does not
# pull in its dependency, so import it here before Acts is searched.
find_package(Mille QUIET)
find_package(Acts REQUIRED)


# Include the sources for the external dependencies.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")


# --- pybind11 dependency ---
set(_acts_pybind11_version 3.0.1)
include(ExternSources)

add_subdirectory(thirdparty/pybind11)

# Require Python interpreter and development headers
#find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development)
#add_subdirectory(/home/pibutti/sw/acts/thirdparty/pybind11)   #OUT OF TREE


get_property(_all_targets GLOBAL PROPERTY TARGETS)
message(STATUS "All CMake targets currently defined:")
foreach(t IN LISTS _all_targets)
  message(STATUS "  ${t}")
endforeach()

# ----------------------------------------
# Add subdirectories
# ----------------------------------------
add_subdirectory(ActsAlgorithms)
