import numba
import math as mt
import itertools
import numpy as np
from scipy.integrate import quad
from scipy import constants

clight = constants.c/1000. # km/s

class cosmo:
    def __init__(self,
                 Omega_m=0.315,
                 Omega_b=0.048,
                 Little_h=0.6731,
		 Sigma_8 = 0.828,
                 Ns=0.96):
        self.Omega_m = Omega_m
        self.Omega_b = Omega_b 
        self.Omega_l = 1.0 - Omega_m
        self.Little_h = Little_h
        self.Big_h = Little_h * 100.
        self.Sigma_8 = Sigma_8
        self.Ns = Ns


@numba.jit(nopython=True)
def InvEz(z, H0, om, ol):
    '''Cosmological 1/E(z) function to integrate.'''
    return 1.0/(mt.sqrt(om*(1.0+z)**3 + ol)*H0)


def Dc(z, c, distance='Mpc'):
    '''Loop over z and integrate InvEz
    z: redshifts
    c: cosmology,
    distance in ["Mpc", "Mpc_over_little_h"]'''
    assert distance in ['Mpc', 'Mpc_over_little_h']
    n = len(z)
    zeros = itertools.repeat(0.0, n)
    f_toIntegrate = itertools.repeat(InvEz, n)
    Hub, mat, lam = c.Big_h, c.Omega_m, c.Omega_l
    args = itertools.repeat((Hub, mat, lam), n)
    res = map(quad, f_toIntegrate, zeros, z, args)
    res = np.array(list(res))
    res = res[:, 0]*clight
    if distance == 'Mpc':
        return res
    elif distance == 'Mpc_over_little_h':
        return res * c.Little_h
