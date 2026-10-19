# About this geometry

- Geometry (.root, .gdml) is the same as in `geometry_2026_06_10_segmented_withoutRICH_noEndcapTOFs`, was created with the following O2 command:
```
o2-sim-serial-run5 -n 1 -g pythia8hi -m A3IP TRK FT3 TF3   --configKeyValues "TRKBase.layoutVD=kIRISFullCyl;TRKBase.layoutMLOT=kSegmented;FT3Base.layoutFT3=kSegmentedStave;IOTOFBase.segmentedInnerTOF=true;IOTOFBase.segmentedOuterTOF=true;IOTOFBase.enableBackwardTOF=false;IOTOFBase.enableForwardTOF=false"
```
- The cave volume is increased in `z` to avoid cutting of the last OT disks, services, oTOF.
  
- *The main difference*: the ACTS tracking geometry is now implemented using the "Generation 3" (*Gen3*) approach (earlier geometries were produced with *Gen1*)

- Digitization (`digi-segm_config_pitch20um_with_barrel_TOFs-gen3.json`) is done with segmentation into pixels (tracker VD, ML, OT layers, iTOF, oTOF), which is more realistic that the simple hit smearing. 

- Material is mapped on sensitive **AND passive** surfaces, see figures in `geom/plots/material_mapping_validation`.
