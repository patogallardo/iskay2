import numpy as np
import healpy as hp


def classify_grid(df, Nside):
    '''Receives a dataframe with ra, dec and bins it
    in a healpix grid.'''
    decs = df.dec.values
    ras = df.ra.values
    indx = hp.pixelfunc.ang2pix(Nside,
                                np.radians(-decs+90.),
                                np.radians(360.-ras))
    df1 = df.copy()
    df1['JK_index'] = indx
    return df1


def healpix_histogram_catalog(df1, Nside):
    '''Returns a vector in healpix format with
    per bin counts for the given catalog df.'''
    Npix = hp.nside2npix(Nside)
    df = df1.copy()
    m = np.zeros(Npix)
    df = classify_grid(df, Nside)
    cnts = df['JK_index'].value_counts()
    indices = cnts.index
    counts = cnts.values
    m[indices] = counts
    return m
