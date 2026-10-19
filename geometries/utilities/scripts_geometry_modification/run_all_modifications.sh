#!/bin/bash

#### python scripts for geom modifications (Feb 2026)

# keep the initial version of O2 geometry
python convert_root_to_gdml.py
DIR_ORIGINAL_GEOM=original_O2_geometry
mkdir ${DIR_ORIGINAL_GEOM}
mv o2sim_geometry.root ${DIR_ORIGINAL_GEOM}
cp o2sim_geometry.gdml ${DIR_ORIGINAL_GEOM}

echo ">>> running geometry modifications...."
python 0_shrink_cave_volume.py
#python 1_rename_ITOFSensor_make_TOF_sensors_thinner.py
# python 2_change_material_by_pattern.py
python 3_shift_z_of_layers_slightly_randomly_for_Geant.py

#echo ">>> generating FT3 disks...."

#python 4_generate_FT3_rings.py --cylinders "8, 35,77.0; 8, 35,100.0; 8, 35,122.0;    8, 35,-77.0; 8, 35,-100.0; 8, 35,-122.0;     15,71,150.0; 15,71,180.0; 15,71,220;  15,71,-151.0; 15,71,-180.0; 15,71,-220" \
#	--out ft3_all_rings_within_maxR_optimized_x_slices.gdml  \
#	--select full --fill optimize --y0-samples 256 \
#	--tile-x 5.005 --tile-y 12.81  --sensor-x 5.0 --sensor-y 12.8 

#echo ">>> merging FT3 disks...."

#cp o2sim_geometry.gdml ${DIR_ORIGINAL_GEOM}/o2sim_geometry_TMP_BEFORE_MERGING_OF_DISKS_for_QA.gdml
#python 5_merger_FT3_to_main_geometry.py
#cp ft3_all_rings_within_maxR_optimized_x_slices.gdml ${DIR_ORIGINAL_GEOM}

python3 restructure_gdml.py  o2sim_geometry.gdml  o2sim_geometry_restructured.gdml

mv o2sim_geometry_restructured.gdml o2sim_geometry.gdml

python convert_gdml_to_root.py

### TMP: not needed after fixes in O2! (Feb 13, 2026)
# python 3_bis_fix_y_for_TRKSensor_cylinder_positions.py
# python 3_bisbis_fix_z_length_for_barrel_ML.py
