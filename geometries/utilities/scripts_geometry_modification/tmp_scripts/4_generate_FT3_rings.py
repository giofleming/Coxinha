#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GDML generator: tile annuli (cylinders with inner hole) using axis-aligned rectangular tiles in x-y.

Update requested:
- For each vertical column (fixed x), if there's a "hole" near y=0 due to rmin criterion,
  tiles in that column must still be symmetric w.r.t. x-axis, BUT the first tile may be
  at ANY distance from x-axis (as long as rmin constraint is satisfied).
  => Implemented as: for each |x|-column we choose a free starting offset y0 in [0, dy)
     and then place tiles at y = ±(y0 + k*dy). This keeps x-axis symmetry while allowing
     the first tile to "jump over" the inner hole.

Symmetries enforced:
- x-axis symmetry: (x,y) => (x,-y)
- y-axis symmetry: (+x,y) and (-x,y) columns identical (same y0 choice and same y list)

Fill mode (--fill):
  - grid      : y0 = 0 for all columns (still uses symmetric ±(y0+k*dy))
  - stagger   : y0 = 0 for even |ix|, y0 = dy/2 for odd |ix|
  - optimize  : for each |x| choose y0 in [0,dy) (sampled) maximizing kept tiles in that column
               under the current --select criterion, using the symmetric ±(y0+k*dy) construction.

Selection mode (--select):
  - center    : keep if center is within [rmin, rmax]
  - full      : keep if the ENTIRE rectangle is within [rmin, rmax]
  - no_inner  : keep if rectangle does NOT intersect inner hole (r < rmin), and center <= rmax

Multiple cylinders via --cylinders "rmin,rmax,z; rmin,rmax,z; ..."

Naming/size requirements:
- Volumes/solids:
  FT3Sensor_tile, FT3Chip_tile_Ci, FT3Layer_tile_Ci
  boxTileAir, boxTileChip, boxTileSensor
- Defaults:
  tile_x=2.5025 cm, tile_y=3.2025 cm
  tile_x_sensor=2.5 cm, tile_y_sensor=3.2 cm, sensor_z=0.01 cm
- Sensor z-shift: sensor sits on the tile face closer to z=0 (global origin),
  depending on cylinder z_global.

Empty tags are written as <.../> (NO space before '/>') by manual XML writing.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from typing import List, Tuple, Literal, Iterable, Set


SelectMode = Literal["center", "full", "no_inner"]
FillMode = Literal["grid", "stagger", "optimize"]


@dataclass(frozen=True)
class CylinderSpec:
    rmin: float
    rmax: float
    z_global: float  # cm


@dataclass(frozen=True)
class Params:
    # Tile thickness along z (cm)
    z_thickness: float = 0.11

    # Tile box dimensions (cm): FULL lengths for GDML <box x= y= z=>
    tile_x: float = 2.5025
    tile_y: float = 3.2025

    # Sensor sizes (cm): FULL lengths for GDML <box x= y= z=>
    tile_x_sensor: float = 2.5
    tile_y_sensor: float = 3.2
    sensor_z: float = 0.01

    # Materials
    mat_air: str = "FT3_AIR"
    mat_si: str = "FT3_SILICON"
    mat_tech_si: str = "FT3_Technical_SILICON"

    # Output file
    out_gdml: str = "ft3_multi_rings.gdml"

    # World size margin (cm)
    world_margin_r: float = 10.0
    world_margin_z: float = 10.0

    # Tile selection mode
    select: SelectMode = "center"

    # Fill mode
    fill: FillMode = "optimize"

    # For fill=optimize: number of y0 samples in [0, dy)
    y0_samples: int = 128


def fmt_num(x: float) -> str:
    s = f"{x:.6f}".rstrip("0").rstrip(".")
    return "0" if s in ("-0", "") else s


def name_xy(x: float, y: float) -> str:
    def one(prefix: str, v: float) -> str:
        sign = "m" if v < 0 else ""
        s = fmt_num(abs(v)).replace(".", "p")
        return f"{prefix}{sign}{s}"
    return f"{one('x', x)}_{one('y', y)}"


def parse_cylinders(text: str) -> List[CylinderSpec]:
    specs: List[CylinderSpec] = []
    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue
        fields = [f.strip() for f in part.split(",")]
        if len(fields) != 3:
            raise ValueError(f"Bad cylinder spec '{part}'. Expected 'rmin,rmax,z'.")
        rmin, rmax, z = map(float, fields)
        if rmin < 0 or rmax <= rmin:
            raise ValueError(f"Bad radii in '{part}': require 0 <= rmin < rmax.")
        specs.append(CylinderSpec(rmin=rmin, rmax=rmax, z_global=z))
    if not specs:
        raise ValueError("No cylinders provided.")
    return specs


def sensor_z_shift_for_cylinder(z_global: float, tile_z: float, sensor_z: float) -> float:
    dz = (tile_z - sensor_z) / 2.0
    return +dz if z_global < 0 else -dz


def rect_all_corners_radii(xc: float, yc: float, hx: float, hy: float) -> Tuple[float, float]:
    corners = (
        (xc - hx, yc - hy),
        (xc - hx, yc + hy),
        (xc + hx, yc - hy),
        (xc + hx, yc + hy),
    )
    rs = [math.hypot(x, y) for (x, y) in corners]
    return (min(rs), max(rs))


def rect_intersects_inner_disk(xc: float, yc: float, hx: float, hy: float, rmin: float) -> bool:
    xmin, xmax = xc - hx, xc + hx
    ymin, ymax = yc - hy, yc + hy
    x_clamp = min(max(0.0, xmin), xmax)
    y_clamp = min(max(0.0, ymin), ymax)
    return math.hypot(x_clamp, y_clamp) < rmin


def keep_tile(
    x: float,
    y: float,
    rmin: float,
    rmax: float,
    tile_x: float,
    tile_y: float,
    mode: SelectMode,
) -> bool:
    r_center = math.hypot(x, y)
    hx = tile_x / 2.0
    hy = tile_y / 2.0

    if mode == "center":
        return (rmin <= r_center <= rmax)

    if mode == "full":
        rmin_corner, rmax_corner = rect_all_corners_radii(x, y, hx, hy)
        return (rmin_corner >= rmin) and (rmax_corner <= rmax)

    if mode == "no_inner":
        if r_center > rmax:
            return False
        return not rect_intersects_inner_disk(x, y, hx, hy, rmin)

    raise ValueError(f"Unknown select mode: {mode}")


def count_column_symmetric(
    x: float,
    y0: float,
    rmin: float,
    rmax: float,
    dy: float,
    tile_x: float,
    tile_y: float,
    select: SelectMode,
) -> int:
    """
    Count tiles in a column using y = ±(y0 + k*dy), k>=0, and (optionally y=0 if y0==0).
    """
    # How far in y we need to check
    # Using rmax as a hard bound is enough, but add a bit of margin.
    kmax = int(math.ceil((rmax + dy) / dy)) + 2

    cnt = 0
    for k in range(0, kmax + 1):
        y = y0 + k * dy
        if y == 0.0:
            if keep_tile(x, 0.0, rmin, rmax, tile_x, tile_y, select):
                cnt += 1
        else:
            okp = keep_tile(x, +y, rmin, rmax, tile_x, tile_y, select)
            okm = keep_tile(x, -y, rmin, rmax, tile_x, tile_y, select)
            # Selection should be symmetric, but require both to enforce exact symmetry
            if okp and okm:
                cnt += 2
    return cnt


def choose_y0_for_column_absx(
    ix_abs: int,
    dx: float,
    dy: float,
    rmin: float,
    rmax: float,
    tile_x: float,
    tile_y: float,
    select: SelectMode,
    fill: FillMode,
    y0_samples: int,
) -> float:
    """
    Choose y0 in [0,dy) for the column at x=±ix_abs*dx.

    grid    : y0=0
    stagger : y0=0 or dy/2 by |ix| parity
    optimize: scan y0 in [0,dy) and choose argmax of count_column_symmetric()
    """
    if fill == "grid":
        return 0.0
    if fill == "stagger":
        return 0.0 if (ix_abs % 2 == 0) else (dy / 2.0)
    if fill != "optimize":
        raise ValueError(f"Unknown fill mode: {fill}")

    x = ix_abs * dx
    samples = max(2, int(y0_samples))
    best_y0 = 0.0
    best_cnt = -1

    # scan y0 in [0, dy)
    for s in range(samples):
        y0 = (s * dy) / samples
        cnt = count_column_symmetric(x, y0, rmin, rmax, dy, tile_x, tile_y, select)
        if cnt > best_cnt:
            best_cnt = cnt
            best_y0 = y0

    return best_y0


def iter_tiles_for_cylinder_symmetric_columns(
    rmin: float,
    rmax: float,
    dx: float,
    dy: float,
    tile_x: float,
    tile_y: float,
    select: SelectMode,
    fill: FillMode,
    y0_samples: int,
) -> Iterable[Tuple[float, float]]:
    """
    Yield (x,y) centers with:
      - y-axis symmetry: column +x mirrors column -x (same y0 and y list)
      - x-axis symmetry: within each column, y positions are paired ±(...)
    And crucially:
      - y0 is free (for optimize) so the first tile can be at any distance from x-axis
        (as long as rmin criterion is satisfied), handling the inner-hole "gap" cleanly.
    """
    ix_max = int(math.floor(rmax / dx))
    seen: Set[Tuple[str, str]] = set()

    for ix_abs in range(0, ix_max + 1):
        x_pos = ix_abs * dx
        y0 = choose_y0_for_column_absx(
            ix_abs, dx, dy, rmin, rmax, tile_x, tile_y, select, fill, y0_samples
        )

        kmax = int(math.ceil((rmax + dy) / dy)) + 2
        y_list: List[float] = []

        for k in range(0, kmax + 1):
            y = y0 + k * dy
            if y == 0.0:
                if keep_tile(x_pos, 0.0, rmin, rmax, tile_x, tile_y, select):
                    y_list.append(0.0)
            else:
                # enforce x-axis symmetry explicitly: require both +y and -y to pass
                if keep_tile(x_pos, +y, rmin, rmax, tile_x, tile_y, select) and keep_tile(
                    x_pos, -y, rmin, rmax, tile_x, tile_y, select
                ):
                    y_list.append(+y)
                    y_list.append(-y)

        # Place for x=+x_pos and x=-x_pos (if ix_abs>0)
        for x in ([x_pos] if ix_abs == 0 else [x_pos, -x_pos]):
            for y in sorted(set(y_list)):
                key = (fmt_num(x), fmt_num(y))
                if key in seen:
                    continue
                seen.add(key)
                yield (x, y)


def build_gdml(p: Params, cylinders: List[CylinderSpec]) -> Tuple[str, int, List[int]]:
    lines: List[str] = []

    def add(s: str = ""):
        lines.append(s)

    max_r = max(c.rmax for c in cylinders)
    max_abs_z = max(abs(c.z_global) for c in cylinders)
    world_x = (max_r + p.world_margin_r) * 2.0
    world_y = (max_r + p.world_margin_r) * 2.0
    world_z = (max_abs_z + p.z_thickness / 2.0 + p.world_margin_z) * 2.0

    add('<?xml version="1.0" encoding="UTF-8"?>')
    add('<gdml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"')
    add('      xsi:noNamespaceSchemaLocation="http://cern.ch/2001/Schemas/GDML/gdml.xsd">')
    add("")

    # -------- define --------
    add("  <define>")
    add('    <position name="FT3object_zeropos" x="0" y="0" z="0" unit="cm"/>')

    for ci, c in enumerate(cylinders, start=1):
        add(f'    <position name="FT3Ring_C{ci}_pos" x="0" y="0" z="{fmt_num(c.z_global)}" unit="cm"/>')
        zshift = sensor_z_shift_for_cylinder(c.z_global, p.z_thickness, p.sensor_z)
        add(f'    <position name="FT3Sensor_inChip_C{ci}_pos" x="0" y="0" z="{fmt_num(zshift)}" unit="cm"/>')

    total_tiles = 0
    per_cyl_counts: List[int] = []

    for ci, c in enumerate(cylinders, start=1):
        count = 0
        for (x, y) in iter_tiles_for_cylinder_symmetric_columns(
            rmin=c.rmin,
            rmax=c.rmax,
            dx=p.tile_x,
            dy=p.tile_y,
            tile_x=p.tile_x,
            tile_y=p.tile_y,
            select=p.select,
            fill=p.fill,
            y0_samples=p.y0_samples,
        ):
            tag = name_xy(x, y)
            add(f'    <position name="FT3Tile_C{ci}_{tag}_pos" x="{fmt_num(x)}" y="{fmt_num(y)}" z="0" unit="cm"/>')
            count += 1
        per_cyl_counts.append(count)
        total_tiles += count

    add("  </define>")
    add("")

    # -------- materials (minimal; remove if your full geometry already defines these) --------
    add("  <materials>")
    add(f'    <material name="{p.mat_air}" state="gas">')
    add('      <D value="0.001205" unit="g/cm3"/>')
    add('      <fraction n="1.0" ref="G4_AIR"/>')
    add("    </material>")
    add(f'    <material name="{p.mat_si}">')
    add('      <D value="2.33" unit="g/cm3"/>')
    add('      <fraction n="1.0" ref="G4_Si"/>')
    add("    </material>")
    add(f'    <material name="{p.mat_tech_si}">')
    add('      <D value="2.33" unit="g/cm3"/>')
    add('      <fraction n="1.0" ref="G4_Si"/>')
    add("    </material>")
    add("  </materials>")
    add("")

    # -------- solids --------
    add("  <solids>")
    add(f'    <box name="WorldBox" x="{fmt_num(world_x)}" y="{fmt_num(world_y)}" z="{fmt_num(world_z)}" lunit="cm"/>')
    add(f'    <box name="boxTileAir" x="{fmt_num(p.tile_x)}" y="{fmt_num(p.tile_y)}" z="{fmt_num(p.z_thickness)}" lunit="cm"/>')
    add(f'    <box name="boxTileChip" x="{fmt_num(p.tile_x)}" y="{fmt_num(p.tile_y)}" z="{fmt_num(p.z_thickness)}" lunit="cm"/>')
    add(f'    <box name="boxTileSensor" x="{fmt_num(p.tile_x_sensor)}" y="{fmt_num(p.tile_y_sensor)}" z="{fmt_num(p.sensor_z)}" lunit="cm"/>')
    for ci, c in enumerate(cylinders, start=1):
        add(
            f'    <tube name="FT3RingTube_C{ci}" rmin="{fmt_num(c.rmin)}" rmax="{fmt_num(c.rmax)}" '
            f'z="{fmt_num(p.z_thickness)}" startphi="0" deltaphi="360" aunit="deg" lunit="cm"/>'
        )
    add("  </solids>")
    add("")

    # -------- structure --------
    add("  <structure>")
    add("")

    add('    <volume name="FT3Sensor_tile">')
    add(f'      <materialref ref="{p.mat_si}"/>')
    add('      <solidref ref="boxTileSensor"/>')
    add("    </volume>")
    add("")

    for ci, c in enumerate(cylinders, start=1):
        add(f'    <volume name="FT3Chip_tile_C{ci}">')
        add(f'      <materialref ref="{p.mat_tech_si}"/>')
        add('      <solidref ref="boxTileChip"/>')
        add(f'      <physvol name="FT3Sensor_in_Chip_C{ci}" copynumber="1">')
        add('        <volumeref ref="FT3Sensor_tile"/>')
        add(f'        <positionref ref="FT3Sensor_inChip_C{ci}_pos"/>')
        add("      </physvol>")
        add("    </volume>")
        add("")

        add(f'    <volume name="FT3Layer_tile_C{ci}">')
        add(f'      <materialref ref="{p.mat_air}"/>')
        add('      <solidref ref="boxTileAir"/>')
        add(f'      <physvol name="FT3Chip_in_Layer_C{ci}" copynumber="1">')
        add(f'        <volumeref ref="FT3Chip_tile_C{ci}"/>')
        add('        <positionref ref="FT3object_zeropos"/>')
        add("      </physvol>")
        add("    </volume>")
        add("")

        add(f'    <assembly name="FT3V_C{ci}">')
        copyn = 1
        for (x, y) in iter_tiles_for_cylinder_symmetric_columns(
            rmin=c.rmin,
            rmax=c.rmax,
            dx=p.tile_x,
            dy=p.tile_y,
            tile_x=p.tile_x,
            tile_y=p.tile_y,
            select=p.select,
            fill=p.fill,
            y0_samples=p.y0_samples,
        ):
            tag = name_xy(x, y)
            pos_name = f"FT3Tile_C{ci}_{tag}_pos"
            phys_name = f"FT3Layer_tile_C{ci}_{copyn:06d}"
            add(f'      <physvol name="{phys_name}" copynumber="{copyn}">')
            add(f'        <volumeref ref="FT3Layer_tile_C{ci}"/>')
            add(f'        <positionref ref="{pos_name}"/>')
            add("      </physvol>")
            copyn += 1
        add("    </assembly>")
        add("")

        add(f'    <volume name="FT3RingVolume_C{ci}">')
        add(f'      <materialref ref="{p.mat_air}"/>')
        add(f'      <solidref ref="FT3RingTube_C{ci}"/>')
        add(f'      <physvol name="FT3V_in_ring_C{ci}" copynumber="1">')
        add(f'        <volumeref ref="FT3V_C{ci}"/>')
        add('        <positionref ref="FT3object_zeropos"/>')
        add("      </physvol>")
        add("    </volume>")
        add("")

    add('    <volume name="World">')
    add(f'      <materialref ref="{p.mat_air}"/>')
    add('      <solidref ref="WorldBox"/>')
    for ci, _ in enumerate(cylinders, start=1):
        add(f'      <physvol name="FT3Ring_in_World_C{ci}" copynumber="{ci}">')
        add(f'        <volumeref ref="FT3RingVolume_C{ci}"/>')
        add(f'        <positionref ref="FT3Ring_C{ci}_pos"/>')
        add("      </physvol>")
    add("    </volume>")
    add("")
    add("  </structure>")
    add("")

    add('  <setup name="default" version="1.0">')
    add('    <world ref="World"/>')
    add("  </setup>")
    add("")
    add("</gdml>")
    add("")

    return "\n".join(lines), total_tiles, per_cyl_counts


def main():
    ap = argparse.ArgumentParser(
        description="Generate GDML tiling tiles in multiple annular cylinders (symmetric columns).")

    ap.add_argument(
        "--cylinders",
        default="20,70,0.0",
        help='Semicolon-separated list of "rmin,rmax,z_cm". Example: "20,70,0.0; 20,70,1.0; 25,60,-2.5"',
    )
    ap.add_argument("--out", default=Params().out_gdml, help="Output GDML filename.")

    ap.add_argument("--zthick", type=float, default=Params().z_thickness, help="Tile/ring thickness along z (cm).")
    ap.add_argument("--tile-x", type=float, default=Params().tile_x, help="Tile size along x (cm).")
    ap.add_argument("--tile-y", type=float, default=Params().tile_y, help="Tile size along y (cm).")

    ap.add_argument("--sensor-x", type=float, default=Params().tile_x_sensor, help="Sensor size along x (cm).")
    ap.add_argument("--sensor-y", type=float, default=Params().tile_y_sensor, help="Sensor size along y (cm).")
    ap.add_argument("--sensor-z", type=float, default=Params().sensor_z, help="Sensor thickness along z (cm).")

    ap.add_argument(
        "--select",
        choices=["center", "full", "no_inner"],
        default=Params().select,
        help=(
            "Tile selection mode in annulus: "
            "'center' (center in [rmin,rmax]); "
            "'full' (entire rectangle inside [rmin,rmax]); "
            "'no_inner' (rectangle does not intersect inner hole, center within outer)."
        ),
    )

    ap.add_argument(
        "--fill",
        choices=["grid", "stagger", "optimize"],
        default=Params().fill,
        help=(
            "Fill mode under enforced symmetries: "
            "'grid' (y0=0), "
            "'stagger' (y0=0/dy/2 by |ix| parity), "
            "'optimize' (choose y0 in [0,dy) to maximize tiles per |x| column)."
        ),
    )

    ap.add_argument(
        "--y0-samples",
        type=int,
        default=Params().y0_samples,
        help="For --fill optimize: number of y0 samples in [0,dy) per |x| column (default 128).",
    )

    args = ap.parse_args()

    cylinders = parse_cylinders(args.cylinders)
    p = Params(
        out_gdml=args.out,
        z_thickness=args.zthick,
        tile_x=args.tile_x,
        tile_y=args.tile_y,
        tile_x_sensor=args.sensor_x,
        tile_y_sensor=args.sensor_y,
        sensor_z=args.sensor_z,
        select=args.select,
        fill=args.fill,
        y0_samples=max(2, args.y0_samples),
    )

    gdml, total_tiles, per_cyl = build_gdml(p, cylinders)
    with open(p.out_gdml, "w", encoding="utf-8", newline="\n") as f:
        f.write(gdml)

    print(f"Wrote: {p.out_gdml}")
    print(f"Selection mode: {p.select}")
    print(f"Fill mode: {p.fill} (columns use y=±(y0+k*dy); y0 is free for optimize; symmetry about x and y enforced)")
    for i, (c, n) in enumerate(zip(cylinders, per_cyl), start=1):
        zshift = sensor_z_shift_for_cylinder(c.z_global, p.z_thickness, p.sensor_z)
        print(f"  C{i}: r=[{c.rmin},{c.rmax}] cm, z_global={c.z_global} cm -> tiles={n}, sensor_z_shift={fmt_num(zshift)} cm")
    print(f"Total tiles placed: {total_tiles}")
    print(f"Creating .root based on this .gdml...")

    # import ROOT
    # ROOT.TGeoManager.Import(p.out_gdml)
    # from pathlib import Path
    # path = Path(p.out_gdml)
    # ROOT.gGeoManager.Export(path.stem+".root")


if __name__ == "__main__":
    main()
