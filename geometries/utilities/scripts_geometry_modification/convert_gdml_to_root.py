#!/usr/bin/env python3


# from math import *
import ROOT
print("converting .gdml to .root...")
ROOT.TGeoManager.Import("o2sim_geometry.gdml")
ROOT.gGeoManager.Export("o2sim_geometry.root")
exit()
