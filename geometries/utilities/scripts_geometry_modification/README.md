This set of scripts modifies the .root geometry, tuning it for the ACTS (see [details in WP1 report](https://indico.cern.ch/event/1649872/#2-update-on-geometry-steps-tow) ).

How to run: `bash run_all_modifications.sh`

The current chain (Feb 13, 2026, can be changed in the future):
* Convert `.root` to `.gdml`, and modify `.gdml` via a set of python scripts:
*  `rename_ITOFSensor_make_TOF_sensors_thinner.py`  - to "absorb" iTOF into the barrel tracking volume; make TOF sensors thinner (for better hits in Geant4)
*  `change_material_names_by_pattern.py` -- for Geant4, to form hits from Silicon only in Sensors (not chips, modules, …)
*  `shift_z_of_cylindrical_layers_slightly_randomly_for_Geant.py` -- to let Geant distinguish the cylinders (and match with ACTS surfaces)
*  `generate_FT3_rings.py` -- tiling of disks with “modules” (a temporary solution for FT3-for-ACTS; to be changed, namely, FT3 should be updated in O2)
*  `merger_FT3_to_main_geometry.py` -- merge gdml with FT3 rings to the main gdml, then run `convert_root_to_gdml.py`

