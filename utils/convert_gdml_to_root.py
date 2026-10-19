#!/usr/bin/env python3


import argparse
import pathlib
import ROOT

parser = argparse.ArgumentParser(
    description="Simple converter gdml->root")

parser.add_argument(
    "--gdml",
    "-g",
    help="Input gdml file",
    type=pathlib.Path,
    default=pathlib.Path.cwd() / "o2sim_geometry.gdml"
)

parser.add_argument(
    "--root",
    "-r",
    help="Output root file",
    type=pathlib.Path,
    default=pathlib.Path.cwd() / "o2sim_geometry.root"
)


args = parser.parse_args()

# from math import *
ROOT.TGeoManager.Import(str(args.gdml))
ROOT.gGeoManager.Export(str(args.root))

# exit()
