#!/bin/bash
# NOTE : Quote it else use array to avoid problems #
FILES="./sbatch_*.sh"
for f in $FILES
do
  echo "Launching $f ..."
  # take action on each file. $f store current file name
  sbatch "$f"
done
