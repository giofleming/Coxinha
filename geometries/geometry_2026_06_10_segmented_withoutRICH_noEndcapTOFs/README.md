# About this geometry

- Geometry was created with the following O2 command:
```
o2-sim-serial-run5 -n 1 -g pythia8hi -m A3IP TRK FT3 TF3   --configKeyValues "TRKBase.layoutVD=kIRISFullCyl;TRKBase.layoutMLOT=kSegmented;FT3Base.layoutFT3=kSegmentedStave;IOTOFBase.segmentedInnerTOF=true;IOTOFBase.segmentedOuterTOF=true;IOTOFBase.enableBackwardTOF=false;IOTOFBase.enableForwardTOF=false"
```

- TOF endcaps are removed.

- RICH is absent

- Digitization (`digi-segm_config_pitch20um_with_barrel_TOFs_RICH.json`) is now done with segmentation into pixels (tracker VD, ML, OT layers, iTOF, oTOF), which is more realistic that the simple hit smearing. 

- Material is mapped on sensitive surfaces.

- GDML version of the geometry was "restructred" (with the [restructure_gdml.py]([actso2/geometries/utilities/scripts_geometry_modification/restructure_gdml.py](https://gitlab.cern.ch/alice3-trackers/wp1-simulationsandperformances/actso2/-/blob/master/geometries/utilities/scripts_geometry_modification/restructure_gdml.py?ref_type=heads)) script), which allows:
  - to reduce the .gdml size by ~20 times,
  - to speedup Geant4 simulations by a factor ~7.
  - The modified .gdml was converted back to .root file (it was checked that sim+rec results are the same as those obtained with original gdml and root files).

