from Corrfunc.mocks.DDsmu_mocks import DDsmu_mocks
import numpy as np
from . import tzav
import os
import pandas as pd
import time

RECIPES_NUM_DEN = {'ferreira99': {'num': ['cdT'], 'den': ['c2']}}

class Recipe:
    '''Contains all information to define a
    recipe for computing the pairwise ksz
    estimator.'''
    def __init__(self, recipe_name):
        self.recipe_name = recipe_name
        self.recipe_num_den = RECIPES_NUM_DEN[recipe_name]
        weight_names = []
        for key in RECIPES_NUM_DEN[recipe_name].keys():
            weight_names += RECIPES_NUM_DEN[recipe_name][key]
        self.weight_names = weight_names
# every recipe contains a numerator and a denominator.
# Final estimator is the sum of the numerators divided in
# the sum of the denominators

def pw_ksz(df, params, rc):
    ra_deg = df.ra.values
    dec_deg = df.dec.values
    d_mpc_over_h = df["d_mpc_over_h"].values
    dtype = ra_deg.dtype
    
    ap_photo_sel_string = "dT_%1.2f_arcmin" % params["R_DISK_ARCMIN_PAIRWISEKSZ"]
    dT = df[ap_photo_sel_string].values
    r_bins_mpc_over_h = (np.array(params["R_BINS_MPC"]) * params["LITTLE_H"]).astype(dtype)
    
    autocorr = 1
    cosmology = 2
    nthreads = 2 # this should come from rc
    mumax = 1.0
    mubins = 1
    isa='avx512f'
    verbose=False
    is_comoving_dist = True
    
    if params["REDSHIFT_CORRECTION"]:
        sigma_z = params["SIGMA_Z"]
        dT_ksz = tzav.correct_dT_tzav(dT, df.z.values, sigma_z)
    else:
        dT_ksz = dT

    res1 = DDsmu_mocks(autocorr, cosmology, nthreads,
                      mumax, mubins, r_bins_mpc_over_h,
                      ra_deg, dec_deg, d_mpc_over_h,
                      weights1=dT_ksz,
                      weight_type = 'cdT',
                      is_comoving_dist=is_comoving_dist,
                      c_api_timer=False, isa=isa, verbose=True)
    res2 = DDsmu_mocks(autocorr, cosmology, nthreads,
                     mumax, mubins, r_bins_mpc_over_h,
                     ra_deg, dec_deg, d_mpc_over_h,
                     weights1=dT_ksz,
                     weight_type='c2',
                     is_comoving_dist=is_comoving_dist,
                     c_api_timer=False, isa=isa)
    npairs = res1['npairs']
    s = 0.5 * (res1['smin'] + res2['smax'])
    pksz = res1['weightavg']/res2['weightavg']
    return s, pksz, npairs


def pw_compute_ksz(df, params, rc):
    '''Computes the sums specified in the
    recipe.
    '''
    if rc['VERBOSE']:
        t1 = time.time()
    ra_deg = df.ra.values
    dec_deg = df.dec.values
    d_mpc_over_h = df['d_mpc_over_h'].values
    dtype = ra_deg.dtype

    ap_photo_sel_string = "dT_%1.2f_arcmin" % params["R_DISK_ARCMIN_PAIRWISEKSZ"]
    dT = df[ap_photo_sel_string].values
    r_bins_mpc_over_h = (np.array(params["R_BINS_MPC"]) * params["LITTLE_H"]).astype(dtype)

    # params for corrfunc
    autocorr = 1
    cosmology = 2
    if 'login' in os.uname().nodename:
        nthreads = 2
    else:
        nthreads = rc['NPROC_PAIRWISE']
    mumax = 1.0
    mubins = 1
    is_comoving_dist = True
    isa = "avx512f" # fastest
    verbose=False
    # end params for corrfunc
    
    if params["REDSHIFT_CORRECTION"]:
        sigma_z = params["SIGMA_Z"]
        dT_ksz = tzav.correct_dT_tzav(dT, df.z.values, sigma_z)
    else:
        dT_ksz = dT

    recipe = Recipe(params["ESTIMATOR_NAME"])
    sums = []
    for weight_type in recipe.weight_names:
        res = DDsmu_mocks(autocorr, cosmology, nthreads,
                          mumax, mubins, r_bins_mpc_over_h,
                          ra_deg, dec_deg, d_mpc_over_h,
                          weights1=dT_ksz,
                          weight_type = weight_type,
                          is_comoving_dist=is_comoving_dist,
                          verbose=verbose,
                          c_api_timer=False,
                          isa=isa)
        sums.append(res['weightavg'])
    npairs = res['npairs']
    sums.append(npairs)
    s = 0.5 * (res['smin'] + res['smax'])
    sums.append(s)
    names = recipe.weight_names + ['npairs', 'r_mp_over_h']
    
    l = {}
    for j, name in enumerate(names):
        l[name] = sums[j]
    df_pw = pd.DataFrame(l)
    
    kszcurve = compute_curve(df_pw, params)
    df_pw['ksz_curve'] = kszcurve
    df_pw['r_mp'] = s/params['LITTLE_H']
    if rc['VERBOSE']:
        t2 = time.time()
        print("BS time: %1.2f" %(t2-t1))
    return df_pw


def compute_curve(df_pw, params):
    recipe = Recipe(params["ESTIMATOR_NAME"])
    num_names = recipe.recipe_num_den['num']
    den_names = recipe.recipe_num_den['den']

    num = np.zeros(len(df_pw))
    den = np.zeros(len(df_pw))

    for name in num_names:
        num += df_pw[name].values
    for name in den_names:
        den += df_pw[name].values
    return num/den
