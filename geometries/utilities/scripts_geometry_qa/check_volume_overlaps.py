import ROOT
ROOT.TGeoManager.Import("o2sim_geometry.root")
ROOT.gGeoManager.CheckOverlaps(0.001)
ROOT.gGeoManager.PrintOverlaps()
exit()
