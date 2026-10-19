# ---------------------------------------------------------------------------
#  O2Interface
#
#  This library is to interface ACTS with O2.
#  The first example is the conversion of ACTS tracks to O2 tracks.
#  Or the preparation of AO2Ds.
#
# ---------------------------------------------------------------------------

message(STATUS "")
message(STATUS "============================================================")
message(STATUS "                 COMPILING O2Interface LIBRARY               ")
message(STATUS "============================================================")
message(STATUS "")

add_library(O2Interface SHARED)

# ----------------------------------------
# Sources
# ----------------------------------------
target_sources(O2Interface
    PRIVATE
        src/ActsToO2Track.cpp
        src/LutRootData.cpp
        src/FlatLutEntry.cpp
    PUBLIC
        FILE_SET HEADERS
        BASE_DIRS include
        FILES include/ActsToO2Track.hpp
)

# ----------------------------------------
# Dependencies
# ----------------------------------------

find_package(Acts QUIET COMPONENTS PluginRoot)
find_package(nlohmann_json CONFIG QUIET)
find_package(Eigen3 REQUIRED)

target_include_directories(O2Interface
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)


target_link_libraries(O2Interface
    PUBLIC
    Acts::Core
    Acts::PluginRoot
    Eigen3::Eigen
    ROOT::Core
    ROOT::Geom
    nlohmann_json::nlohmann_json
)

target_compile_features(O2Interface PUBLIC cxx_std_20)

# ----------------------------------------
# Install
# ----------------------------------------
include(GNUInstallDirs)
install(TARGETS O2Interface
                LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
                )

install(
    TARGETS O2Interface
    FILE_SET HEADERS
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/O2Interface"
)