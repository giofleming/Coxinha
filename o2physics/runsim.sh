#!/bin/bash

# This script runs the simulation and analysis for the O2 framework.
# It should be run after rungen.sh to generate the events first but you can pipe the two workflows if needed.

# $HOME/alice/ninja.sh

sleep 1

reset

pushd SimOutput

OutputTag=${1:-""}

OPT="--configuration json://configuration.json"
OPT="-b $OPT"

# Tracking
# o2-analysis-on-the-fly-detector-geometry-provider --configuration json://configuration.json --severity fatal |
#     o2-analysis-alice3-tracking-performance --configuration json://configuration.json --severity fatal |
#     o2-analysis-onthefly-richpid $OPT |
#     o2-analysis-onthefly-tofpid $OPT |
#     o2-analysis-alice3-efficiency --configuration json://configuration.json |
#     o2-analysis-onthefly-tracker --configuration json://configuration.json --aod-file AnalysisResults_trees.root --severity info 

o2-analysis-on-the-fly-detector-geometry-provider $OPT |
    o2-analysis-onthefly-tracker $OPT |
    o2-analysis-alice3-tracking-performance $OPT |
    o2-analysis-alice3-efficiency $OPT --aod-file AnalysisResults_trees.root

mv AnalysisResults.root AnalysisResults${OutputTag}.root

popd
