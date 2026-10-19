#!/bin/bash

### run full sim+rec chain (ok for simple geometries, but slow for complex geoms, especially with Geant4)
python alice3_full_chain.py --out_dir_prefix sim_rec_Pythia_PU30 -n100  --nThreads 1 --config job_configs/pp.py --runMode full_chain

### run simulation only (slow, but can be done one time); `--nThreads 1` allows to avoid issues with Geant4 multi-threading
# python alice3_full_chain.py --out_dir_prefix sim_Pythia_PU30 -n100  --nThreads 1 --config job_configs/pp.py --runMode sim_only

### run reconstruction using avaiable simulation output (faster, but requires sim output to be available)
# python alice3_full_chain.py --out_dir_prefix rec_Pythia_PU30 -n100  --config job_configs/pp.py --runMode reco_only  --simdir output/sim_Pythia_PU30

### run multiple reconstruction passes with varying configurations, see example list of varied parameters in alice3_multi_reco_passes.py (sim output must be available)
# python alice3_multi_reco_passes.py --out_dir_prefix rec_MULTI_PASSES_Pythia_PU30 -n100  --config job_configs/pp.py --runMode reco_only  --simdir output/sim_pp_PU30_nEv100
