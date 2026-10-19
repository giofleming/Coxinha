#!/usr/bin/env python3


# from math import *
import ROOT
print("converting .root to .gdml...")
ROOT.TGeoManager.Import("o2sim_geometry.root")
ROOT.gGeoManager.Export("o2sim_geometry.gdml")
exit()
