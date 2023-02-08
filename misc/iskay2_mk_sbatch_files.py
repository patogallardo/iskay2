#! /usr/bin/env python
import glob
import iskay2
import os
from iskay2 import paramfile

fnames = glob.glob('params_*.json')

def mk_sbatchfile(param_fname):
    param = paramfile.load_paramfile(param_fname)

    string = ['#!/bin/bash']
    string.append('#SBATCH --nodes=1')
    string.append('#SBATCH --ntasks-per-node=40')
    string.append('#SBATCH --time=3:00:00')
    string.append('#SBATCH --job-name=%s' % param["NAME"])
    string.append('#SBATCH --output=mpi_output_%j.txt')

    string.append('module load NiaEnv/2019b python/3.8.5')
    string.append('source ~/.virtualenvs/env/bin/activate')

    string.append('export PYTHONPATH=$PYTHONPATH:$HOME/code/iskay2')
    string.append('export PATH=$PATH:$HOME/code/iskay2/misc')
    string.append('export PATH=$PATH:$HOME/.local/bin')

    iskay2dir = os.path.split(os.path.split(iskay2.__file__)[0])[0]
    miscdir = os.path.join(iskay2dir, 'misc')

    string.append('python %s/iskay2_compute_pwksz_and_bootstrap.py %s' % (miscdir, param_fname))

    string = "\n".join(string)
    return string

for param_fname in fnames:
    string = mk_sbatchfile(param_fname)
    param = paramfile.load_paramfile(param_fname)
    print("Making batch file for %s" % param['NAME'])
    with open('sbatch_%s.sh' % param['NAME'], 'w') as f:
        f.write(string)
