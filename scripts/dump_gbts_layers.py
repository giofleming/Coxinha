#!/usr/bin/env python3
"""Dump the ALICE 3 sensitive-layer inventory in the form GBTS needs.

For every (volume, layer) of the ALICE 3 tracking geometry that carries
sensitive surfaces, prints the number of surfaces, the mean reference
coordinate and the true (r, z) envelope, handling cylinders, discs and planes
separately (see docs/GBTS.md section 3).

Usage:
    python3 scripts/dump_gbts_layers.py <o2sim_geometry.root> [--json out.json]
"""

import argparse
import json
import sys

import numpy as np
import acts

from AliceActsPythonBindings import buildALICE3GeometryGen3


def surfaceExtent(surface, gctx):
    """Return (rmin, rmax, zmin, zmax, refCoordR, refCoordZ) of one surface."""
    values = list(surface.bounds.values())
    centre = np.array(surface.center(gctx))
    stype = str(surface.type).split(".")[-1]

    if stype == "Cylinder":
        # CylinderBounds::values() = (R, halfZ, halfPhi, avgPhi, bevelMin, bevelMax)
        radius, halfZ = values[0], values[1]
        return (radius, radius, float(centre[2]) - halfZ, float(centre[2]) + halfZ,
                radius, float(centre[2]))
    if stype == "Disc":
        # RadialBounds::values() = (rMin, rMax, halfPhi, avgPhi)
        return (values[0], values[1], float(centre[2]), float(centre[2]),
                0.5 * (values[0] + values[1]), float(centre[2]))

    # Planar sensor: RectangleBounds::values() = (-hlX, -hlY, hlX, hlY).
    # Push the four corners through localToGlobal; the centre alone badly
    # undersizes a large module.
    momentum = acts.Vector3(1.0, 1.0, 1.0)
    rs, zs = [], []
    for lx, ly in ((values[0], values[1]), (values[2], values[1]),
                   (values[2], values[3]), (values[0], values[3])):
        g = np.array(surface.localToGlobal(gctx, acts.Vector2(lx, ly), momentum))
        rs.append(float(np.hypot(g[0], g[1])))
        zs.append(float(g[2]))
    return (min(rs), max(rs), min(zs), max(zs),
            float(np.hypot(centre[0], centre[1])), float(centre[2]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("geometry", help="path to o2sim_geometry.root")
    parser.add_argument("--json", default=None, help="also write the table as json")
    args = parser.parse_args()

    trackingGeometry = buildALICE3GeometryGen3(
        args.geometry, acts.logging.WARNING, "", True
    )
    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()

    layers = {}

    class Visitor(acts.TrackingGeometryMutableVisitor):
        def visitSurface(self, surface):
            if not surface.isSensitive:
                return
            gid = surface.geometryId
            key = (gid.volume, gid.layer)
            entry = layers.setdefault(key, dict(
                n=0, rmin=1e9, rmax=-1e9, zmin=1e9, zmax=-1e9,
                rsum=0.0, zsum=0.0, type=str(surface.type).split(".")[-1]))
            rmin, rmax, zmin, zmax, rref, zref = surfaceExtent(surface, gctx)
            entry["n"] += 1
            entry["rmin"] = min(entry["rmin"], rmin)
            entry["rmax"] = max(entry["rmax"], rmax)
            entry["zmin"] = min(entry["zmin"], zmin)
            entry["zmax"] = max(entry["zmax"], zmax)
            entry["rsum"] += rref
            entry["zsum"] += zref

    trackingGeometry.apply(Visitor())

    header = (f"{'vol':>4} {'lay':>4} {'n':>6} {'<r>':>8} {'<z>':>9} "
              f"{'rmin':>8} {'rmax':>8} {'zmin':>9} {'zmax':>9}  {'type':<9} kind")
    print(header)
    rows = []
    for (volume, layer) in sorted(layers):
        e = layers[(volume, layer)]
        meanR = e["rsum"] / e["n"]
        meanZ = e["zsum"] / e["n"]
        # GBTS convention: a layer longer in z than deep in r is a barrel
        isBarrel = (e["zmax"] - e["zmin"]) > (e["rmax"] - e["rmin"])
        kind = "Barrel" if isBarrel else "Endcap"
        print(f"{volume:>4} {layer:>4} {e['n']:>6} {meanR:>8.2f} {meanZ:>9.2f} "
              f"{e['rmin']:>8.2f} {e['rmax']:>8.2f} {e['zmin']:>9.2f} "
              f"{e['zmax']:>9.2f}  {e['type']:<9} {kind}")
        rows.append(dict(volume=volume, layer=layer, nSurfaces=e["n"],
                         meanR=meanR, meanZ=meanZ, rmin=e["rmin"], rmax=e["rmax"],
                         zmin=e["zmin"], zmax=e["zmax"], surfaceType=e["type"],
                         gbtsType=kind))

    if args.json is not None:
        with open(args.json, "w") as handle:
            json.dump(rows, handle, indent=1)
        print(f"\nwrote {args.json}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
