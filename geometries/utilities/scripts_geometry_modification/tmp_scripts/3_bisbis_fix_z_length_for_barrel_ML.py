#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import xml.etree.ElementTree as ET

GDML_PATH = "o2sim_geometry.gdml"
TARGET_Z = "124"
PREFIXES = ["TRKSensor", "TRKMetalStack", "TRKChip", "TRKModule", "TRKStave", "TRKLayer"]
XS = range(0, 5)  # numbers for volume names (barrel Middle Layers)


def indent(elem, level=0):
    i = "\n" + level * "  "
    if len(elem):
        if not (elem.text and elem.text.strip()):
            elem.text = i + "  "
        for e in elem:
            indent(e, level + 1)
        if not (elem.tail and elem.tail.strip()):
            elem.tail = i
    else:
        if level and not (elem.tail and elem.tail.strip()):
            elem.tail = i


tree = ET.parse(GDML_PATH)
root = tree.getroot()

struct = root.find("structure")
solids = root.find("solids")
if struct is None or solids is None:
    raise RuntimeError("Expected <structure> and <solids> sections in GDML.")

# index tubes by name
tube_by_name = {t.get("name"): t for t in solids.findall("tube") if t.get("name")}

target_vol_names = {f"{pref}{x}" for pref in PREFIXES for x in XS}

changed = 0
for vol in struct.findall("volume"):
    vname = vol.get("name")
    if vname not in target_vol_names:
        continue

    sref = vol.find("solidref")
    if sref is None:
        continue
    solid_name = sref.get("ref")
    if not solid_name:
        continue

    tube = tube_by_name.get(solid_name)
    if tube is None:
        continue

    tube.set("z", TARGET_Z)
    changed += 1

indent(root)
tree.write(GDML_PATH, encoding="utf-8", xml_declaration=True)

print(f"Done. Updated z=124 for {changed} tube reference(s). File overwritten: {GDML_PATH}")
