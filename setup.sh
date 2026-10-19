#!/usr/bin/env bash

# Source the ACTS setup script to set up the environment
if [ -z "$ACTS_ROOT" ]; then
    echo "ACTS_ROOT is not set. Please set ACTS_ROOT to the root directory of your ACTS installation."
else
    . $ACTS_ROOT/bin/this_acts.sh
fi

# Absolute path to this setup.sh
export MAIN_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"

# Add the directory to the python path
export PYTHONPATH="${MAIN_DIR}:${PYTHONPATH}"

# Add the local python bindings to the python path
if [ -z "$ACTSO2_ROOT" ]; then
    echo "ACTSO2_ROOT is not set. Using local build directory for Python bindings."
    export PYTHONPATH="${MAIN_DIR}/build/ActsAlgorithms/python:${PYTHONPATH}"
else
    echo "Using ACTSO2_ROOT=${ACTSO2_ROOT} for Python bindings."
    export PYTHONPATH="${ACTSO2_ROOT}/python:${PYTHONPATH}"
fi

# Set DYLD_LIBRARY_PATH on macOS
if [[ "$(uname)" == "Darwin" ]]; then
    export DYLD_LIBRARY_PATH="$ACTS_ROOT/lib:$ACTSO2_ROOT/lib:$DYLD_LIBRARY_PATH"
    export SDKROOT=$(xcrun --show-sdk-path)
    export CPATH="$SDKROOT/usr/include/c++/v1"
fi

# Disable FPE errors
export ACTS_SEQUENCER_DISABLE_FPEMON=1

cd $MAIN_DIR
