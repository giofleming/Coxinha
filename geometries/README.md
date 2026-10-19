# About this folder

This folder contains several sub-directories with *various ALICE 3 geometries prepared to for tests with the ACTS sim+rec*.

More precisely, each geometry is actually a set of files: starting from .root and corresponding .gdml, ACTS tracking geometry was created (via alice3.py), and material was scanned (via Geantino scan) and mapped (material-map.json).


Geometries were created with commands:
1) *barrel*: simple cylinders (TGeoTubes); *disks*: O2 disks are tiled in modules
```
o2-sim-serial-run5 -n 1 -g pythia8hi -m A3IP TRK FT3 TF3 --configKeyValues "TRKBase.layoutML=kTurboStaves;TRKBase.layoutOL=kStaggered;"
```
(segmentation for disks: O2 implementation is replaced by simple tiles: .root --> .gdml --> modifying .gdml --> creating new .root)

**--> FOLDER: [ to be uploaded ]**

2) *barrel*: segmentation of ML and OT as in O2; *disks*: O2 disks are replaced by simple tiles (via changing .gdml-->.root):
```
o2-sim-serial-run5 -n 1 -g pythia8hi -m A3IP TRK FT3 TF3 --configKeyValues "TRKBase.layoutML=kCylinder;TRKBase.layoutOL=kCylinder;"
```
(segmentation for disks: O2 implementation is replaced by simple tiles: .root --> .gdml --> modifying .gdml --> creating new .root)

**--> FOLDER: geometry_2026_01_30_mockupTiledDisks_for_Geant_tests**



## Examples of running with various arguments:

The `full_chain_Jan_2026.py` script in each folder runs the full sim+rec chain
(help on script arguments: `python full_chain_Jan_2026.py --help`).

Runnable with ACTS 43.0.1 (Jan-Feb 2026).

 ### Muon gun (`--gunPID 13 --gunRandCharge`) in `--gunEtaRange 0.0 3.0` with Fatras detector response sim:
```
 python full_chain_Jan_2026.py   --out_dir_prefix TEST_Fatras_gun_muons\
    -n1000 --nThreads 1 --field 2.0  \
    --gunMult 1 --gunPtRange 0.9999 1.0  --gunEtaRange 0.0 3.0 --gunPhiRange 0 6.2831853  --gunPID 13 --gunRandCharge \
    --detSim Fatras  --seedingLayers VD --minSeedPt 0.07  \
    --nMeasurementsMin 7 --ckfChi2Measurement 45 --ckfMeasPerSurf 1 --ckfChi2Outlier 100
```

### $\Lambda^0$ gun (`--gunPID 3122`) in `--gunEtaRange 0.0 2.0` with Geant4 response:
```
python full_chain_Jan_2026.py   --out_dir_prefix TEST_Geant_gun_Lambda0\
    -n10 --nThreads 1 --field 2.0  \
    --gunMult 1 --gunPtRange 0.9999 1.0  --gunEtaRange 0.0 2.0 --gunPhiRange 0 6.2831853  --gunPID 3122 \
    --detSim Geant4  --seedingLayers ML3 --minSeedPt 0.07  \
    --impParForSeeds 400 --seedingAlgo Default \
    --nMeasurementsMin 5 --ckfChi2Measurement 45 --ckfMeasPerSurf 1 --ckfChi2Outlier 100      
```
(Note relaxed impact parameter for seeds, MLs for seeds, descreased nMeasurementsMin per track.)


### Pythia pp events with Geant4 detector response:
```
 python full_chain_Jan_2026.py   --out_dir_prefix TEST_Geant_PYTHIA_pp\
    -n200 --nThreads 1 --field 2.0 --usePythia --system pp --pileup 1 \
    --detSim Geant4  --seedingLayers VD --minSeedPt 0.07  \
    --nMeasurementsMin 7 --ckfChi2Measurement 45 --ckfMeasPerSurf 1 --ckfChi2Outlier 100 
```



*Note*: if you have Magnetic Field (MF) map, adding the argument `--fieldMap` and specifying the path, (e.g. `--fieldMap ../magneticFieldMaps/solenoid_R1625_L5500_B2T_scaled.txt`) allows to use it, otherwise, a constant MF (along z) will be used.