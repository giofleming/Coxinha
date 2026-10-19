#!/bin/bash

# This script generates the events and saves them to a file.

reset

# Check O2DPG_MC_CONFIG_ROOT
if [ -z "$O2DPG_MC_CONFIG_ROOT" ]; then
    echo "O2DPG_MC_CONFIG_ROOT is not set. Please set it to the O2DPG MC config root directory."
    exit 1
fi

OUT_DIR="SimOutput"
if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
    ln -s $PWD/configuration.json "$OUT_DIR"/configuration.json
fi

pushd "$OUT_DIR"

# Pb-Pb
# o2-sim-dpl-eventgen --generator external --nEvents 1000 --aggregate-timeframe 50 --configFile $O2DPG_MC_CONFIG_ROOT/MC/config/ALICE3/ini/pythia8_PbPb.ini -b --configuration json://configuration.json |
#     o2-sim-mctracks-to-aod --tf-offset 391403000000 -b --configuration json://configuration.json --aod-writer-keep dangling
# mv AnalysisResults_trees.root AnalysisResults_trees_PbPb.root

# pp
o2-sim-dpl-eventgen --generator external --nEvents 100000 --aggregate-timeframe 50 --configFile $O2DPG_MC_CONFIG_ROOT/MC/config/ALICE3/ini/pythia8_pp.ini -b --configuration json://configuration.json |
    o2-sim-mctracks-to-aod --tf-offset 391403000000 -b --configuration json://configuration.json --aod-writer-keep dangling
mv AnalysisResults_trees.root AnalysisResults_trees_pp.root

popd
