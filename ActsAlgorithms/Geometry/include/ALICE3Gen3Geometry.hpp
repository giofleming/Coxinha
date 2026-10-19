// ---------------------------------------------------------------------------
//  ALICE3Gen3Geometry.hpp
//
//  Gen3 ("Blueprint") tracking-geometry builder for ALICE 3.
//
//  Replaces the Gen1 `TGeoDetector`-based buildALICE3Geometry.py: instead of
//  letting TGeoDetector infer layers/volumes, we parse the TGeo tree, group the
//  sensitive surfaces into layers, assemble a Blueprint tree of container /
//  layer nodes and let ACTS run the Build / Connect / Finalize phases.
//
//  The result is handed to python as a ready-to-use TrackingGeometry, so the
//  ALICE 3 chain only has to swap the geometry-construction line.
//
//  Material is NOT applied yet - see the follow-up step (proto-material nodes
//  via MaterialDesignatorBlueprintNode::configureFace + map loading through
//  TrackingGeometry::apply).
// ---------------------------------------------------------------------------

#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>

namespace Acts {
class TrackingGeometry;
}

namespace AliceActsTrk {

/// Build the ALICE 3 Gen3 tracking geometry from a TGeo (ROOT) geometry file.
///
/// @param rootFile     path to o2sim_geometry.root
/// @param level        ACTS log level used during construction
/// @param graphvizFile if non-empty, the blueprint tree is dumped there (.dot)
/// @param withMaterial if true (default) every layer volume face is marked as
///                     a binned material receiver (proto material), and thin
///                     passive cylinder layers are inserted at the carbon
///                     supports / services / vessel-wall radii so that their
///                     material can be mapped too. Set false for a pure
///                     geometry without any material surfaces.
/// @param configFile   path to the geometry-config JSON. If empty (default),
///                     gen3_geometry_config.json next to @p rootFile is used.
///                     Missing file / required field -> hard error.
///
/// @return the constructed tracking geometry (surfaces reference detector
///         elements that are kept alive internally for the process lifetime)
/// Navigation policy applied to the layer volumes.
///
/// TryAll is the ACTS blueprint default: it enumerates *every* surface confined
/// in a volume on each navigation initialisation. With ~19k sensors on an outer
/// tracker layer and ~54k on the oTOF that is O(10^6) intersection candidates
/// per event, which dominates CKF time - see docs/GBTS.md section 10.
///
/// SurfaceArray replaces that with a binned lookup for the *sensitive* surfaces
/// and keeps the portals and passives (collected once) alongside it, because
/// SurfaceArrayNavigationPolicy indexes sensitive surfaces only and the material
/// in this geometry sits almost entirely on passive surfaces.
enum class Gen3NavigationPolicy {
  TryAll = 0,       ///< ACTS default; current behaviour
  SurfaceArray = 1  ///< binned sensitive lookup + portals/passives
};

std::shared_ptr<Acts::TrackingGeometry> buildALICE3GeometryGen3(
    const std::string& rootFile,
    Acts::Logging::Level level = Acts::Logging::INFO,
    const std::string& graphvizFile = "", bool withMaterial = true,
    const std::string& configFile = "",
    Gen3NavigationPolicy navigationPolicy = Gen3NavigationPolicy::TryAll);

}  // namespace AliceActsTrk
