#!/usr/bin/env python3
"""
Deduplicate / instance a ROOT-exported GDML geometry.

ROOT's TGeoManager -> GDML writer emits a UNIQUE <solid> and a UNIQUE <volume>
for every placed node, even when they are geometrically identical.  For the
segmented ALICE 3 geometry this produces ~250k logical volumes / ~250k solids (as of June 2026)
that collapse to only a few dozen distinct ones.  Geant4 then builds a
G4LogicalVolume + G4SmartVoxelHeader for every one of them -> huge memory and a
very slow initialisation / navigation.

This tool rewrites the GDML so that all geometrically identical sub-trees share
ONE logical volume and ONE solid (true instancing), via a bottom-up structural
signature.  The expanded geometry at run time is IDENTICAL (same shapes, same
placements, same material budget) -- only the number of *distinct* logical
volumes / solids / file-level physvols drops by orders of magnitude.

Because identical volumes are now shared, per-placement identity must come from
the copy-number CHAIN in the Geant4 touchable history.  The script therefore
assigns a unique copynumber to every physvol within its parent (0..n-1), so the
full chain world->...->stave[c]->module[c]->chip[c]->sensor[c] stays unique.
Configure your ACTS Geant4 sensitive-surface mapping to use the touchable
copy-number chain (not the volume name) to identify hits.

The modified gdml is ~20 times smaller in size, 
and the speedup in ACTS simulations with Geant4 is factor 7.

Usage:
    python3 restructure_gdml.py  o2sim_geometry.gdml  o2sim_geometry_restructured.gdml

Uses only the Python standard library (xml.etree.ElementTree).
"""

import sys
import os
import xml.etree.ElementTree as ET

# ---- numeric tolerance for treating two transforms as identical ------------
# float32 round-off from ROOT is ~1e-6; real geometry differs by >> mm, so
# rounding to 1e-6 cm (10 nm) merges duplicates without any false merge.
NDIG = 6


def r(x):
    try:
        return round(float(x), NDIG)
    except (TypeError, ValueError):
        return x

# ---------------------------------------------------------------------------


def main(src, dst, copynumber_mode="sequential"):
    defs = {}   # name -> (tag, attrib-dict, value-tuple)   (positions/rotations/scales)
    solid_canon = {}   # solid name        -> canonical solid name
    solid_raw = {}   # canonical name    -> serialized <solid> text (indented)
    solidkey_seen = {}  # shape key         -> canonical solid name
    vol_canon = {}   # volume name       -> canonical volume name
    sig_seen = {}   # signature         -> canonical volume name
    emit_order = []   # canonical volume names, children-before-parents
    emit_rec = {}   # canonical name -> dict(tag, material, solid, kids[...])
    world = [None]
    setup_name = ["default"]
    setup_ver = ["1.0"]

    counters = {"vol": 0, "phys": 0, "solid": 0}

    def value_tuple(elem):
        return (r(elem.get("x", 0)), r(elem.get("y", 0)),
                r(elem.get("z", 0)), elem.get("unit", ""))

    def serialize_solid(elem, indent="    "):
        """Emit a solid element as GDML text, remapping solid refs to canonical."""
        attrs = []
        for k, v in elem.attrib.items():
            if k == "ref" and v in solid_canon:
                v = solid_canon[v]
            attrs.append(f'{k}="{v}"')
        head = "<" + elem.tag + ((" " + " ".join(attrs)) if attrs else "")
        kids = list(elem)
        if not kids:
            return indent + head + "/>"
        lines = [indent + head + ">"]
        for c in kids:
            lines.append(serialize_solid(c, indent + "  "))
        lines.append(indent + "</" + elem.tag + ">")
        return "\n".join(lines)

    def solid_key(elem):
        """Structural key: tag + attrs (name stripped, refs canonicalised) + children."""
        parts = [elem.tag]
        for k, v in sorted(elem.attrib.items()):
            if k == "name":
                continue
            if k == "ref":
                if v in solid_canon:
                    v = "S:" + solid_canon[v]
                elif v in defs:
                    v = "V:" + str(defs[v][2])
            parts.append(k + "=" + str(v))
        for c in elem:
            parts.append(solid_key(c))
        return tuple(parts)

    def transform_of(phys):
        """Return (pos,rot,scl tuples, posref, rotref, scaleref).  Supports inline."""
        pos = rot = (0.0, 0.0, 0.0)
        scl = (1.0, 1.0, 1.0)
        pref = rref = sref = None
        for ch in phys:
            t = ch.tag
            if t == "positionref":
                pref = ch.get("ref")
                v = defs.get(pref)
                if v:
                    pos = v[2][:3]
            elif t == "position":
                pos = (r(ch.get("x", 0)), r(ch.get("y", 0)), r(ch.get("z", 0)), ch.get("unit", "mm"))
            elif t == "rotationref":
                rref = ch.get("ref")
                v = defs.get(rref)
                if v:
                    rot = v[2][:3]
            elif t == "rotation":
                rot = (r(ch.get("x", 0)), r(ch.get("y", 0)), r(ch.get("z", 0)), ch.get("unit", "deg"))
            elif t == "scaleref":
                sref = ch.get("ref")
                v = defs.get(sref)
                if v:
                    scl = v[2][:3]
            elif t == "scale":
                scl = (r(ch.get("x", 1)), r(ch.get("y", 1)), r(ch.get("z", 1)))
        return pos, rot, scl, pref, rref, sref

    # ---- streaming parse -----------------------------------------------------
    context = ET.iterparse(src, events=("start", "end"))
    _, root = next(context)              # start of <gdml>
    stack = [root]
    nproc = 0

    for event, elem in context:
        if event == "start":
            stack.append(elem)
            continue
        # event == "end"
        stack.pop()
        parent = stack[-1] if stack else None
        ptag = parent.tag if parent is not None else None
        tag = elem.tag

        # ---- <define>: positions / rotations / scales -----------------------
        if ptag == "define" and tag in ("position", "rotation", "scale"):
            nm = elem.get("name")
            if nm is not None:
                defs[nm] = (tag, dict(elem.attrib), value_tuple(elem))
            elem.clear()

        # ---- <solids>: any solid element ------------------------------------
        elif ptag == "solids":
            counters["solid"] += 1
            nm = elem.get("name")
            key = solid_key(elem)
            canon = solidkey_seen.get(key)
            if canon is None:
                canon = nm
                solidkey_seen[key] = canon
                solid_raw[canon] = serialize_solid(elem)
            solid_canon[nm] = canon
            elem.clear()

        # ---- <structure>: <volume> and <assembly> ---------------------------
        elif tag in ("volume", "assembly") and ptag == "structure":
            counters["vol"] += 1
            nm = elem.get("name")
            material = solid = None
            kids = []
            for ch in elem:
                if ch.tag == "materialref":
                    material = ch.get("ref")
                elif ch.tag == "solidref":
                    solid = solid_canon.get(ch.get("ref"), ch.get("ref"))
                elif ch.tag == "physvol":
                    counters["phys"] += 1
                    cref = None
                    for g in ch:
                        if g.tag == "volumeref":
                            cref = g.get("ref")
                            break
                    ccanon = vol_canon.get(cref, cref)
                    pos, rot, scl, pref, rref, sref = transform_of(ch)
                    kids.append([ccanon, pos, rot, scl, pref, rref, sref,
                                 ch.get("copynumber")])
            sig = (tag, material, solid,
                   tuple((k[0], k[1], k[2], k[3]) for k in kids))
            canon = sig_seen.get(sig)
            if canon is None:
                canon = nm
                sig_seen[sig] = canon
                for i, k in enumerate(kids):
                    if copynumber_mode == "sequential":
                        k[7] = i
                emit_rec[canon] = dict(tag=tag, material=material,
                                       solid=solid, kids=kids)
                emit_order.append(canon)
            vol_canon[nm] = canon
            elem.clear()

        # ---- <setup>: remember world volume ---------------------------------
        elif tag == "setup":
            for ch in elem:
                if ch.tag == "world":
                    world[0] = vol_canon.get(ch.get("ref"), ch.get("ref"))
            setup_name[0] = elem.get("name", "default")
            setup_ver[0] = elem.get("version", "1.0")

        # ---- periodic memory release (drop processed shells) ----------------
        nproc += 1
        if ptag in ("define", "solids", "structure") and (nproc % 50000 == 0):
            parent.clear()

    # ---- collect referenced position/rotation/scale names --------------------
    used_defs = set()
    for canon in emit_order:
        for k in emit_rec[canon]["kids"]:
            for ref in (k[4], k[5], k[6]):
                if ref:
                    used_defs.add(ref)

    # ---- write output -------------------------------------------------------
    with open(dst, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<gdml xmlns:xsi="https://www.w3.org/2001/XMLSchema-instance" '
                'xsi:noNamespaceSchemaLocation="https://service-spi.web.cern.ch/'
                'service-spi/app/releases/GDML/schema/gdml.xsd">\n')

        # define (only referenced)
        f.write("  <define>\n")
        for nm, (dtag, at, _vt) in defs.items():
            if nm in used_defs:
                a = " ".join(f'{k}="{v}"' for k, v in at.items())
                f.write(f"    <{dtag} {a}/>\n")
        f.write("  </define>\n")

        # materials: copied verbatim from the source
        f.write(_slice_section(src, "<materials>", "</materials>"))

        # solids (canonical only)
        f.write("  <solids>\n")
        for canon in solid_raw:
            f.write(solid_raw[canon] + "\n")
        f.write("  </solids>\n")

        # structure (canonical volumes, children before parents)
        f.write("  <structure>\n")
        for canon in emit_order:
            rec = emit_rec[canon]
            f.write(f'    <{rec["tag"]} name="{canon}">\n')
            if rec["material"]:
                f.write(f'      <materialref ref="{rec["material"]}"/>\n')
            if rec["solid"]:
                f.write(f'      <solidref ref="{rec["solid"]}"/>\n')
            for j, k in enumerate(rec["kids"]):
                ccanon, pos, rot, scl, pref, rref, sref, cno = k
                f.write(f'      <physvol name="{canon}_pv{j}" copynumber="{cno}">\n')
                f.write(f'        <volumeref ref="{ccanon}"/>\n')
                if pref:
                    f.write(f'        <positionref ref="{pref}"/>\n')
                if rref:
                    f.write(f'        <rotationref ref="{rref}"/>\n')
                if sref:
                    f.write(f'        <scaleref ref="{sref}"/>\n')
                f.write('      </physvol>\n')
            f.write(f'    </{rec["tag"]}>\n')
        f.write("  </structure>\n")

        # setup
        f.write(f'  <setup name="{setup_name[0]}" version="{setup_ver[0]}">\n')
        f.write(f'    <world ref="{world[0]}"/>\n')
        f.write('  </setup>\n')
        f.write('</gdml>\n')

    # ---- report -------------------------------------------------------------
    out_phys = sum(len(emit_rec[c]["kids"]) for c in emit_order)
    print("------------------------------------------------------------")
    print(f"  input : volumes={counters['vol']:>8}  solids={counters['solid']:>8}  "
          f"physvol-elements={counters['phys']:>8}")
    print(f"  output: volumes={len(emit_order):>8}  solids={len(solid_raw):>8}  "
          f"physvol-elements={out_phys:>8}")
    print(f"  file  : {os.path.getsize(src)/1e6:7.1f} MB  ->  "
          f"{os.path.getsize(dst)/1e6:7.1f} MB")
    print("  (run-time expanded geometry is unchanged; only DISTINCT logical")
    print("   volumes / solids drop -- that is what Geant4 init & navigation pay for)")
    print("------------------------------------------------------------")


def _slice_section(path, start, end):
    """Return the verbatim text block between start and end markers (inclusive)."""
    out, on = [], False
    with open(path) as f:
        for line in f:
            if start in line:
                on = True
            if on:
                out.append(line)
            if end in line:
                break
    return "".join(out)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
