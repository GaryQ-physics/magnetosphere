#!/bin/bash
# Push or pull files derived from MHD run files to mag.gmu.edu
# Usage:
#   ./rsync.sh [push|pull] USERNAME [test|]
#   where USERNAME is account name on mag.gmu.edu
# Example:
#   ./rsync.sh push weigel test # Test what push would do
#   ./rsync.sh push weigel      # Do actual push

# Master directory is 
#   /media/disk/git-data/sblake
# on mag.gmu.edu. Directory is visible at
#   http://mag.gmu.edu/git-data/sblake

if [ "$3" == "test" ]; then
    dry="--dry-run"
    echo "Dry run"
fi

# Push files from current directory to master directory
if [ "$1" == "push" ]; then
    rsync -avz -L $dry --filter "exclude *~" --exclude SCARR5_GM_IO2 --filter="exclude *.sh" \
	. $2@mag.gmu.edu:/media/disk/git-data/sblake/SCARR5_GM_IO2-derived
fi

# Pull files from master directory to current directory
# Use the command --delete-before to delete any files in current directory
# that do not exist on remote directory
if [ "$1" == "pull" ]; then
    rsync -avz $dry --filter "exclude *.out SCARR5_GM_IO2" \
	$2@mag.gmu.edu:/media/disk/git-data/sblake/SCARR5_GM_IO2-derived .
fi

