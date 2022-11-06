#!/bin/bash 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=3:00:00
#SBATCH --job-name=ap_photo
#SBATCH --output=mpi_output_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pgallardo@uchicago.edu
 
module load NiaEnv/2019b python/3.8.5
source ~/.virtualenvs/env/bin/activate

export PYTHONPATH=$PYTHONPATH:$HOME/code/iskay2
export PATH=$PATH:$HOME/code/iskay2/misc
export PATH=$PATH:$HOME/.local/bin

python $HOME/code/iskay2/misc/iskay2_extract_ap_photo.py
