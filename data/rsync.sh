#!/bin/bash
# Usage:
#   ./rsync.sh [push|pull] USERNAME [test|]
#   where USERNAME is account name on mag.gmu.edu
# Example:
#   ./rsync.sh push weigel test # Test what push would do
#   ./rsync.sh push weigel      # Do actual push

# Master directory is 
#   /media/disk/git-data/GaryQ-Physics/magnetosphere
# on mag.gmu.edu. Directory is visible at
#   http://mag.gmu.edu/git-data/GaryQ-Physics/magnetosphere/

if [ "$3" == "test" ]; then
    dry="--dry-run"
    echo "Dry run"
fi

# Push files from current directory to master directory
if [ "$1" == "push" ]; then
    rsync -avz $dry --filter "exclude *~" --exclude SCARR5_GM_IO2 --exclude="rsync.sh" \
	. $2@mag.gmu.edu:/media/disk/git-data/GaryQ-Physics/magnetosphere
fi

# Pull files from master directory to current directory
# Use the command --delete-before to delete any files in current directory
# that do not exist on remote directory
if [ "$1" == "pull" ]; then
    rsync -avz $dry --filter "exclude *.out SCARR5_GM_IO2" \
	$2@mag.gmu.edu:/media/disk/git-data/GaryQ-Physics/magnetosphere/ .
fi

