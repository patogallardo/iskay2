#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --nodes=1
source $HOME/code/niagara_scripts/load_python.sh
export NCORES=79
echo 'Hello, world!'

python $HOME/code/iskay2/test/srun_pairwiseksz.py
