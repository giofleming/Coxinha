#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import random
import re
from pathlib import Path


# 1 micron = 1e-4 cm (since 1 cm = 10,000 microns)
MICRON_TO_CM = 1e-4
input_gdml = "o2sim_geometry.gdml"


def build_name_matcher():
    patterns = [
        re.compile(r"^TRKLayer\d+_1inTRKVpos$"),
        re.compile(
            r"^PETALCASE0_LAYER\d+_TRKSensor\d+_1inPETALCASE0_LAYER\d+_TRKChip\d+pos$"
        ),
        re.compile(r"^OTOFLayer_1inIOTOFVpos$"),
        # ITOFLayer is not shifted: it ends flush against the A-side
        # TRK_MIDBARCONN service disk at z = 64.5 cm, so any +z shift overlaps.
    ]

    def matches(name: str) -> bool:
        return any(p.match(name) for p in patterns)

    return matches


def main():
    ap = argparse.ArgumentParser(
        description="Randomize z for selected <position .../> tags in GDML (Gaussian, sigma in microns)."
    )
    # ap.add_argument("input", type=Path, default=Path(input_gdml), help="Input GDML file")
    input = Path(input_gdml)
    output = Path(input_gdml)
    # ap.add_argument("output", type=Path, default=Path(input_gdml), help="Output GDML file")
    ap.add_argument("--seed", type=int, default=None, help="Random seed (optional)")
    ap.add_argument(
        "--sigma-micron", type=float, default=20.0, help="Gaussian sigma in microns"
    )
    ap.add_argument(
        "--only-unit-cm",
        action="store_true",
        help='If set, modify only tags with unit="cm" (recommended).',
    )
    args = ap.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    sigma_cm = args.sigma_micron * MICRON_TO_CM

    name_matches = build_name_matcher()

    # Regex: capture prefix up to z=" , then current z value, then closing quote
    # Also capture name="...".
    # Works even if attributes are in different order, as long as name and z are on the same line.
    pos_re = re.compile(
        r'(<position\b[^>]*\bname="(?P<name>[^"]+)"[^>]*\bz=")(?P<z>[^"]*)(")',
        re.IGNORECASE,
    )

    # Optional check for unit="cm"
    unit_cm_re = re.compile(r'\bunit="cm"\b', re.IGNORECASE)

    text = input.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)

    changed = 0
    for i, line in enumerate(text):
        m = pos_re.search(line)
        if not m:
            continue

        name = m.group("name")
        if not name_matches(name):
            continue

        if args.only_unit_cm and not unit_cm_re.search(line):
            continue

        new_z = random.gauss(0.0, sigma_cm)
        new_z_str = f"{new_z:.10f}"

        # Replace only the z value (preserve rest of line)
        new_line = line[: m.start("z")] + new_z_str + line[m.end("z"):]
        text[i] = new_line
        changed += 1

    # Regex to match the tube line with TRKSensor0/1/2 and z="50"
    pattern = re.compile(
        r'(<tube[^>]*name="[^"]*TRKSensor[012][^"]*"[^>]*\bz=")50(")',
    )
    # Replace z="50" to z="49.9" for VD sensor tubes
    new_content = pattern.sub(r"\g<1>49.9\2", "".join(text))

    # Replace all http:// with https://
    new_content = new_content.replace("http://", "https://")

    output.write_text(new_content, encoding="utf-8")
    print(f"Done. Updated z in {changed} <position> lines. Output: {output}")


if __name__ == "__main__":
    main()
