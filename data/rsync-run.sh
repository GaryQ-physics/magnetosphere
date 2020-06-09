#!/bin/bash
# Download MHD run files from mag.gmu.edu
# Usage:
#   ./rsync-run.sh USERNAME [test]
#   where USERNAME is account name on mag.gmu.edu
# Example:
#   ./rsync weigel test
#   ./rsync weigel
#
# Alternative:
# wget -m -nP --cut-dirs=2 http://mag.gmu.edu/git-data/sblake/

RUN=SCARR5_GM_IO2

if [ "$2" == "test" ]; then
    dry="--dry-run"
    echo "Dry run"
fi

rsync -avz $dry --exclude "*.out" $1@mag.gmu.edu:/media/disk/git-data/sblake/$RUN $RUN

