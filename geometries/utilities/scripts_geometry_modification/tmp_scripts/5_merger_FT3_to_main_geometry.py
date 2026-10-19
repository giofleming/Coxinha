#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import copy
import re
import xml.etree.ElementTree as ET
from pathlib import Path

MERGE_SECTIONS = ["define", "materials", "solids", "userinfo"]
STRUCTURE_TAG = "structure"

FT3_WORLD_NAME = "World"
FT3_PHYSVOL_NAME_PREFIX = "FT3Ring_in_World_C"

FT3_ASSEMBLY_NAME = "FT3V"
BARREL_VOLUME_NAME = "barrel"
BARREL_INSERT_PHYSVOL_NAME = "FT3V_2"
BARREL_POSITIONREF = "TRKV_2inbarrelpos"


def indent(elem: ET.Element, level: int = 0):
    """Pretty-print indentation (minimal)."""
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


def index_named_children(parent: ET.Element):
    """Set of (tag, name) for direct children that have @name."""
    s = set()
    for ch in list(parent):
        nm = ch.get("name")
        if nm is not None:
            s.add((ch.tag, nm))
    return s


def ensure_section(root: ET.Element, tag: str) -> ET.Element:
    sec = root.find(tag)
    if sec is None:
        sec = ET.SubElement(root, tag)
    return sec


def merge_section_by_name(dst_root: ET.Element, src_root: ET.Element, tag: str) -> int:
    src_sec = src_root.find(tag)
    if src_sec is None:
        return 0

    dst_sec = ensure_section(dst_root, tag)
    existing = index_named_children(dst_sec)

    added = 0
    for ch in list(src_sec):
        nm = ch.get("name")
        key = (ch.tag, nm) if nm is not None else None
        if key is not None and key in existing:
            continue
        dst_sec.append(copy.deepcopy(ch))
        if key is not None:
            existing.add(key)
        added += 1
    return added


def merge_structure_except_world_and_setup(dst_root: ET.Element, src_root: ET.Element) -> int:
    """
    Merge everything under <structure> from src into dst, but:
      - skip <volume name="World"> from src
      - skip <setup> from src (if any)
    Deduplicate by (tag, name) where possible.
    """
    src_struct = src_root.find(STRUCTURE_TAG)
    if src_struct is None:
        return 0

    dst_struct = ensure_section(dst_root, STRUCTURE_TAG)
    existing = index_named_children(dst_struct)

    added = 0
    for ch in list(src_struct):
        # Skip src setup entirely
        if ch.tag.lower() == "setup":
            continue

        # Skip src World volume
        if ch.tag == "volume" and (ch.get("name") == FT3_WORLD_NAME):
            continue

        nm = ch.get("name")
        key = (ch.tag, nm) if nm is not None else None
        if key is not None and key in existing:
            continue

        dst_struct.append(copy.deepcopy(ch))
        if key is not None:
            existing.add(key)
        added += 1

    return added


def find_volume(structure: ET.Element, name: str) -> ET.Element | None:
    for v in structure.findall("volume"):
        if v.get("name") == name:
            return v
    return None


def find_or_create_assembly(structure: ET.Element, name: str) -> ET.Element:
    for a in structure.findall("assembly"):
        if a.get("name") == name:
            return a
    return ET.SubElement(structure, "assembly", {"name": name})


def collect_ft3_physvols_from_world(src_root: ET.Element) -> list[ET.Element]:
    src_struct = src_root.find(STRUCTURE_TAG)
    if src_struct is None:
        return []

    world = find_volume(src_struct, FT3_WORLD_NAME)
    if world is None:
        return []

    out = []
    rx = re.compile(rf"^{re.escape(FT3_PHYSVOL_NAME_PREFIX)}\d+", re.IGNORECASE)

    for pv in world.findall("physvol"):
        nm = pv.get("name", "")
        if rx.match(nm):
            out.append(copy.deepcopy(pv))

    return out


def ensure_barrel_physvol_insert(dst_root: ET.Element):
    struct = ensure_section(dst_root, STRUCTURE_TAG)
    barrel = find_volume(struct, BARREL_VOLUME_NAME)
    if barrel is None:
        raise RuntimeError(f'No <volume name="{BARREL_VOLUME_NAME}"> in destination.')

    # если уже есть physvol с таким name — не добавляем второй раз
    for pv in barrel.findall("physvol"):
        if pv.get("name") == BARREL_INSERT_PHYSVOL_NAME:
            return False

    pv = ET.SubElement(barrel, "physvol", {"name": BARREL_INSERT_PHYSVOL_NAME})
    ET.SubElement(pv, "volumeref", {"ref": FT3_ASSEMBLY_NAME})
    ET.SubElement(pv, "positionref", {"ref": BARREL_POSITIONREF})
    return True


def check_positionref_exists(dst_root: ET.Element) -> bool:
    define = dst_root.find("define")
    if define is None:
        return False
    for pos in define.findall("position"):
        if pos.get("name") == BARREL_POSITIONREF:
            return True
    return False


def move_barrel_and_cave_to_end(root: ET.Element):
    struct = root.find("structure")
    if struct is None:
        print("No <structure> found")
        return

    # найти нужные элементы
    barrel = None
    cave = None
    ft3v = None

    for ch in list(struct):
        if ch.tag == "volume" and ch.get("name") == "barrel":
            barrel = ch
        elif ch.tag == "volume" and ch.get("name") == "cave":
            cave = ch
        elif ch.tag == "assembly" and ch.get("name") == "FT3V":
            ft3v = ch

    if ft3v is None:
        print("No FT3V assembly found")
        return

    # удалить их из текущего места
    for obj in (barrel, cave):
        if obj is not None and obj in list(struct):
            struct.remove(obj)

    # вставить в конец после FT3V
    children = list(struct)
    idx = children.index(ft3v)

    insert_pos = idx + 1

    if barrel is not None:
        struct.insert(insert_pos, barrel)
        insert_pos += 1

    if cave is not None:
        struct.insert(insert_pos, cave)

    print("Moved barrel and cave after FT3V")


def main():
    # ap = argparse.ArgumentParser(description="Merge ft3 GDML into o2sim geometry and embed FT3 into barrel.")
    # ap.add_argument("dst_gdml", type=Path, help="Base o2sim_geometry.gdml")
    # ap.add_argument("src_gdml", type=Path, help="ft3_multi_rings.gdml")
    # ap.add_argument("out_gdml", type=Path, help="Output merged gdml")
    # args = ap.parse_args()

    SRC = Path("ft3_all_rings_within_maxR_optimized_x_slices.gdml")
    DST = Path("o2sim_geometry.gdml")
    # OUT = Path("o2sim_geometry.gdml")
    OUT = DST

    dst_tree = ET.parse(DST)
    src_tree = ET.parse(SRC)
    dst_root = dst_tree.getroot()
    src_root = src_tree.getroot()

    # 1) Merge standard sections
    for tag in MERGE_SECTIONS:
        n = merge_section_by_name(dst_root, src_root, tag)
        print(f"{tag}: added {n}")

    # 2) Merge structure (skip src World and src setup)
    nstruct = merge_structure_except_world_and_setup(dst_root, src_root)
    print(f"structure: added {nstruct} (excluding src World & setup)")

    # 3) Collect FT3 physvols from src World and put into assembly FT3V
    ft3_physvols = collect_ft3_physvols_from_world(src_root)
    print(f"FT3 physvols collected from src World: {len(ft3_physvols)}")

    dst_struct = ensure_section(dst_root, STRUCTURE_TAG)
    ft3_assembly = find_or_create_assembly(dst_struct, FT3_ASSEMBLY_NAME)

    # avoid re-adding same physvol names if script run twice
    existing_pv_names = {pv.get("name") for pv in ft3_assembly.findall("physvol") if pv.get("name")}
    added_pv = 0
    for pv in ft3_physvols:
        nm = pv.get("name")
        if nm and nm in existing_pv_names:
            continue
        ft3_assembly.append(pv)
        if nm:
            existing_pv_names.add(nm)
        added_pv += 1

    print(f"assembly {FT3_ASSEMBLY_NAME}: appended {added_pv} physvol(s)")

    # 4) Insert FT3V into barrel as FT3V_2 at TRKV_2inbarrelpos
    inserted = ensure_barrel_physvol_insert(dst_root)
    print(f'barrel insert physvol "{BARREL_INSERT_PHYSVOL_NAME}": {"added" if inserted else "already existed"}')

    # sanity check for position
    if not check_positionref_exists(dst_root):
        print(f'WARNING: position "{BARREL_POSITIONREF}" not found in <define>. The reference may be dangling.')

    move_barrel_and_cave_to_end(dst_root)

    # Write out
    indent(dst_root)
    dst_tree.write(OUT, encoding="utf-8", xml_declaration=True)
    print(f"Done. Output: {OUT}")


if __name__ == "__main__":
    main()
