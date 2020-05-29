#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

export PYTHONPATH=$DIR/..:$PYTHOPATH

# Read PARAVIEW variable from config.sh
. $DIR/../config_paths.sh

COM="$PARAVIEW --script=earth.py"

echo "Executing $COM"

$COM