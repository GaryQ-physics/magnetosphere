#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

export PYTHONPATH=$DIR/..:$PYTHONPATH

# Read PARAVIEW variable from config.sh
. $DIR/../config_paths.sh

COM="$PARAVIEW --script=$DIR/earth.py"

echo "Executing $COM"

$COM