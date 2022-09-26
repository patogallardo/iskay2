import json
import os
import numpy as np

def write_paramfile():
    '''Write a default param file in current directory.'''
    data = {"NAME": "DR5_150GHz",
            "R_BINS_MPC":[0.01, 10, 20, 30, 40, 50, 60, 70, 80,
                          90, 100, 110, 120, 130, 140, 150, 200,
                          250, 315, 395],
            "REDSHIFT_CORRECTION": True,
            "SIGMA_Z": 0.01,
            "R_DISK_ARCMIN_PAIRWISEKSZ": 2.1,
            "R_DISKS_ARCMIN": [2.1, 2.5, 3.0],
            "R_RING_OVER_R_DISK": np.sqrt(2.0),
            "MAP_FNAME": "act_planck_s08_s18_cmb_f150_night_map.fits",
            "DIVMAP_FNAME": "act_planck_s08_s18_cmb_f150_night_ivar.fits",
            "CAT_FNAME": "V20_DR15_Catalog.csv",
            "QUERY": 'lum > 7.962650518080777e10 and S18coadd == True',
            "NGAL": 1000000,
            "SORTBY": 'lum',
            "OMEGA_M": 0.315,
            "LITTLE_H": 0.6731,
            "ESTIMATOR_NAME": 'ferreira99',
            
            "N_BOOTSTRAP_ITERATIONS": 100
            }
    fname_out = os.path.join(os.getcwd(), "params.json")
    with open(fname_out, 'w') as f:
        json.dump(data, f, indent=4)


def load_paramfile(fname='params.json'):
    '''Load a param file located in the current directory.
    fname: specifies the filename in current dir.'''
    print("Opening param file: %s" % fname)
    fname = os.path.join(os.getcwd(), fname)
    with open(fname, 'r') as f:
        params = json.load(f)
    print("\nParams:")
    for element in params:
        print("%s: %s" % (element, params[element]))
    return params
