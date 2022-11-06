import json
import os
import numpy as np

def write_paramfile():
    '''Write a default param file in current directory.'''
    data = {"NAME": "DR6_150GHz",
            "R_BINS_MPC":[0.01, 10, 20, 30, 40, 50, 60, 70, 80,
                          90, 100, 110, 120, 130, 140, 150, 200,
                          250, 315, 395],
            "REDSHIFT_CORRECTION": True,
            "SIGMA_Z": 0.01,
            "R_DISK_ARCMIN_PAIRWISEKSZ": 2.1,
            "R_DISKS_ARCMIN": [2.1, 2.5, 3.0],
            "R_RING_OVER_R_DISK": np.sqrt(2.0),
            "MAP_FNAME": "20211219_act_planck_night_f150_map.fits",
            "DIVMAP_FNAME": "20211219_act_planck_night_f150_ivar.fits",
            "CAT_FNAME": "V20_DR15_Catalog_v3.csv",
            "MASK0_FNAME": "act_mask_20220316_GAL020_rms_70.00_downgrade_None.fits",
            "MASK1_FNAME": "act_mask_20220316_GAL040_rms_70.00_downgrade_None.fits",
            "MASK2_FNAME": "act_mask_20220316_GAL060_rms_70.00_downgrade_None.fits",
            "MASK3_FNAME": "srcfind_mask_f150.fits",
            "MASK4_FNAME": "None",
            "MASK5_FNAME": "None",
            "MASKED_RADIUS_ARCMIN": 3.0,
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


def get_mask_fnames(params):
    fnames = []
    names = []
    for j in range(5):
        tag = "MASK%i_FNAME" % j
        if ".fits" in params[tag]:
            fnames.append(params[tag])
            names.append(tag.split("_")[0])
    return fnames, names


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
