#!/bin/bash

# visualization of ACTS tracking geometry
python geometry.py
python  ~/alice/ACTS/Examples/Scripts/MaterialMapping/GeometryVisualisationAndMaterialHandling.py --geometry geometry-map.json

# record material hits via geantino scan
python material_recording.py

# tune binning, map material
python ~/alice/ACTS/Examples/Scripts/MaterialMapping/writeMapConfig.py geometry-map.json config-map.json
python modifyConfigMap.py config-map.json 80 80
python  ~/alice/ACTS/Examples/Scripts/MaterialMapping/configureMap.py geometry-map.json config-map-binned.json
python material_mapping.py

# validation
python material_validation.py
mkdir Validation
root -l -b -q ~/alice/ACTS/Examples/Scripts/MaterialMapping/Mat_map.C'("propagation-material.root","material-map_tracks.root","Validation")'