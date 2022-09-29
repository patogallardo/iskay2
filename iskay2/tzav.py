'''
This module contains the redshift average correction.
'''

import numba
import math as mt
import numpy as np
from scipy import interpolate


# get tzav
@numba.guvectorize(['float64[:],float64[:],float64[:],float64[:],'
                    'float64[:],float64[:]'],
                   '(n),(n),(),()->(),()', target='parallel')
def get_tzav_and_w_nb(dT, z, zj, sigma_z, res1, res2):
    '''Launched by get_tzav to compute formula in parallel '''
    for i in range(dT.shape[0]):
        res1 += dT[i] * mt.exp(-(zj[0]-z[i])**2.0/(2.0*sigma_z[0]**2))
        res2 += mt.exp(-(zj[0]-z[i])**2/(2.0*sigma_z[0]**2))


def get_tzav_fast(dTs, zs, sigma_z):
    '''Subsample and interpolate Tzav to make it fast.
    dTs: entire list of dT decrements
    zs: entire list of redshifts
    sigma_z: width of the gaussian kernel we want to apply.
    '''
    N_samples_in_sigmaz = 15  # in one width of sigmaz use Nsamples
    zmin, zmax = zs.min(), zs.max()
    delta_z = zmax - zmin

    # evaluates Tzav N times
    N_samples = int(round(delta_z/sigma_z)) * N_samples_in_sigmaz
    z_subsampled = np.linspace(zmin, zmax, N_samples)

    #now compute tzav as we usually do.
    res1 = np.zeros(z_subsampled.shape[0])
    res2 = np.zeros(z_subsampled.shape[0])
    get_tzav_and_w_nb(dTs, zs, z_subsampled, sigma_z, res1, res2)
    tzav_subsampled = res1/res2
    #interpolate
    f = interpolate.interp1d(z_subsampled, tzav_subsampled, kind='cubic')
    tzav_fast = f(zs)
    return tzav_fast


def correct_dT_tzav(dT, z, sigma_z):
    '''Return a corrected temperature decrement.
    dT: vector with temperature decrements
    z: vector with redshifts
    sigma_z: window in which to average decrements.
    returns: dT - tzav(z)
    '''
    tzav = get_tzav_fast(dT, z, sigma_z)
    dT_ksz = dT - tzav
    if dT_ksz.dtype != dT.dtype:
        dT_ksz = dT_ksz.astype(dT.dtype)
    return dT_ksz
