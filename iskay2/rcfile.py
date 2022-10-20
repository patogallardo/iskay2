#!/gpfs/fs1/home/r/rbond/gallardo/.virtualenvs/env/bin/python
import json
from pathlib import Path
import os
import numpy as np
from multiprocessing import cpu_count

def write_default_rcfile():
    data = {
            "MAPS_PATH": os.path.join(os.environ.get("SCRATCH"), 'data/maps'),
            "CATS_PATH": os.path.join(os.environ.get("SCRATCH"), 'data/catalogs'),
            "OVERSAMPLE":2,
            "R_ARCMIN_SUBMAP":10,
            "R_RAD_SUBMAP":np.deg2rad(10./60),
            "NPROC_AP_PHOTO": cpu_count(),
            "NPROC_PAIRWISE":cpu_count(),
            "RANDOM_SEED":0,
            }
    fname_out = os.path.join(str(Path.home()), ".iskay2rc")

    with open(fname_out, 'w') as f:
        json.dump(data, f, indent=4)


def load_rcfile():
    fname = os.path.join(str(Path.home()), ".iskay2rc")
    with open(fname, 'r') as f:
        data = json.load(f)
    print("\nRC File Parameters:")
    for element in data:
        print("%s: %s" % (element, data[element]))
    print("\n")
    return data
