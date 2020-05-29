#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

export PYTHONPATH=$DIR/..:$PYTHOPATH

# On OS-X, need to first execute
# cd /usr/local/bin/
# sudo ln -s /Applications/ParaView-5.4.1.app/Contents/MacOS/paraview

paraview --script=earth.py