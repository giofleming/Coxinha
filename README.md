[![run_test](https://github.com/AliceO2Group/actsO2/actions/workflows/run_test.yml/badge.svg)](https://github.com/AliceO2Group/actsO2/actions/workflows/run_test.yml)
[![geometry_test](https://github.com/AliceO2Group/actsO2/actions/workflows/geo_build_materialmap.yml/badge.svg)](https://github.com/AliceO2Group/actsO2/actions/workflows/geo_build_materialmap.yml)


# Standalone ACTS workflow for ALICE 3 studies

## Instructions how to run the full chain simulation and reconstruction.

Full ACTS documentation is located [here](https://acts.readthedocs.io/en/latest/).

# ACTS lxplus installation: 

ACTS installation has been tested with LCG version 108 and 107 (see CI)

## LCG Setup

```
source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc14-opt/setup.sh
```

## Compile and install ACTS

The following instructions provide the minimal installation of acts to run with Alice3 full reconstruction

```
git clone https://github.com/acts-project/acts.git
cd acts
mkdir build
cmake -S acts -B acts/build -DACTS_BUILD_PLUGIN_GEANT4=on \
-DACTS_BUILD_PLUGIN_JSON=on -DACTS_BUILD_PLUGIN_ROOT=on \
-DACTS_BUILD_FATRAS=on -DACTS_BUILD_FATRAS_GEANT4=on \
-DACTS_BUILD_EXAMPLES_GEANT4=on \
-DACTS_BUILD_EXAMPLES=on \
-DACTS_BUILD_PYTHON_BINDINGS=on \
-DACTS_BUILD_ANALYSIS_APPS=on \
-DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
-DACTS_BUILD_EXAMPLES_PYTHIA8=on \
-DCMAKE_INSTALL_PREFIX=acts/install
cd acts/build
make -j install
```

## Setup ACTS 

Setup acts to have it available for python bindings and compilation of this package.
```
source acts/install/bin/this_acts.sh
```


## ACTS aliBuild installation - OLD To be tested:
```
aliBuild init
git clone git@github.com:AliceO2Group/acts.git
aliBuild -d build ACTS
```
### Fixing potential issues with the aliBuild installation
- Issue: `ImportError: libHepMC3.so.4: cannot open shared object file: No such file or directory` (M. Faggin)
  - Fix: `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/path-to-HepMC3/3.3.0-44/lib"`
- Issue: Missing or not found Boost package
  - Fix: Install Boost or edit the $PATH: `export PATH=$PATH:"/path-to-boost/v1.83.0-alice2-local1/include"`
- Issue: `Imported target "HepMC3::rootIO" includes non-existent path`
  - Fix: Declare ROOT as a development package `aliBuild init ROOT`

When the installation is completed you can pull the tutorial repository and place it in ```actsdir```.

# Material maps:
The material-map is provided in .zip format. The full-chain job will unzip it automatically at the first run. 


# Compilation of this package

```
source setup.sh
mkdir build && cd build
cmake -DCMAKE_PREFIX_PATH=<path/to/>acts/install ../
```

# Running the full chain simulation and reconstruction on lxplus

If running with Geant4, only single-thread is supported

```
source setup.sh
python3 alice3_full_chain.py --nEv 100 --nThreads 1
```

The output files will be located in a folder located in ```./output```

# Get the $\it{p}_{\rm{T}}$ resolution plot

We provide a simple macro to run the momentum resolution plots in ```macros```

```
root -l -q 'getPtResolution.C+g("<outputFolder>")'
```

This will create the folders ```Plots``` and ```treeoutput``` in ```reco_output_pythia``` which contain the figure and the ROOT with the $\it{p}_{\rm{T}}$ resolution.

