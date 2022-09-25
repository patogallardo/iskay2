import numba
from .rcfile import load_rcfile
from .pairwiser import angle_jit_rad, vecdiff_jit, inWhatBinIsIt
import math as mt
import numpy as np
from concurrent import futures
import itertools


def gen_groups_start_end_pw_fromrowtorow(Ngal, Ngroups=500):
    '''Divide the rows of a vector to use with
    pairwise_from_rowtorow.
    returns the row_start and row_end numbers.
    Both numbers are inclusive.'''
    if Ngal/4 < Ngroups: # how many groups of 2 can be formed?
        Ngroups = 1
    rowend = row_end(Ngal)
    groups = [ [a[0], a[-1]] for a in 
               np.array_split(np.arange(0, rowend+1), Ngroups)]
    return groups


def row_end(Ngal):
    '''Use this with pairwise_from_rowtorow.
    Computes the end row index that marks the end
    of the calculation if we are doing mirrored
    rows. Note this index is inclusive, and range()
    needs to go to row_end + 1 to be inclusive due
    to zero-indexing.'''
    return mt.ceil((Ngal-1)/2) - 1


def get_nbins(r_max, r_step):
    return mt.ceil(r_max/r_step)


@numba.jit(nopython=True, nogil=True)
def pairwise_from_rowtorow(row_start, row_end,
                           Dc, ra_rad, dec_rad,
                           Tmapsc, r_max, r_step,
                           dTw, w2, w, dT, Npairs):
    '''Compute the pairwise sum balancing the number of computations
       by computing the sum at row i and row N-2-i.
       For the complete matrix row_end = ceil((Ngal-1)/2)
    '''
    Ngal = len(ra_rad)
    for row in range(row_start, row_end + 1):
        pairwise_one_row(row, Ngal, Dc, ra_rad, dec_rad,
                         Tmapsc, r_max, r_step,
                         dTw, w2, w, dT, Npairs)
        mirror_row = Ngal - 2 - row
        if mirror_row > row:
            pairwise_one_row(mirror_row, Ngal, Dc, ra_rad, dec_rad,
                             Tmapsc, r_max, r_step,
                             dTw, w2, w, dT, Npairs)


@numba.jit(nopython=True, nogil=True)  # row_end = roundup(N/2)
def pairwise_from_rowtorow_onlyonerowatatime(row_start, row_end,
                                             Dc, ra_rad, dec_rad,
                                             Tmapsc, r_max, r_step,
                                             dTw, w2, w, dT, Npairs):
    '''Use for diagnostics. Computes the pairwise sum
       inside a for loop from row_start to row_end.'''
    Ngal = len(ra_rad)
    for row in range(row_start, row_end+1):
        pairwise_one_row(row, Ngal, Dc, ra_rad, dec_rad,
                         Tmapsc, r_max, r_step,
                         dTw, w2, w, dT, Npairs)


@numba.jit(nopython=True, nogil=True)
def pairwise_one_row(row, Ngal, Dc, ra_rad, dec_rad,
                     Tmapsc, r_max, r_step,
                     dTw, w2, w, dT, Npairs):
    '''Computes one pairwise row
       inputs: row, Ngal, Dc, ... Tmapsc,
               r_max, r_sep: maximum separation and uniform step size
               dTw, w2, w, dT, Npairs are the partial sum vectors, 
               which should be zero for one row, otherwise the
               result will just add up.
    '''
    for j in range(row+1, Ngal):
        ang_ij = angle_jit_rad(dec_rad[row], ra_rad[row],
                               dec_rad[j], ra_rad[j])
        vecdiff_ij = vecdiff_jit(Dc[row], Dc[j], ang_ij)
        if (vecdiff_ij > 0) and (vecdiff_ij < r_max):
            binval_ij = mt.floor(vecdiff_ij/r_step)
            dT_ij = (Tmapsc[row]) - (Tmapsc[j])
            cij = (Dc[row]-Dc[j])*(1.0 + mt.cos(ang_ij)) / (2.0*vecdiff_ij)
            dTw[binval_ij] += dT_ij * cij
            w2[binval_ij] += cij**2.
            dT[binval_ij] += dT_ij
            w[binval_ij] += cij
            Npairs[binval_ij] += 1


def gen_output_pairwise_one_row(r_max, r_step):
    N_bins = get_nbins(r_max, r_step)
    dTw = np.zeros(N_bins)
    w2 = np.zeros(N_bins)
    w = np.zeros(N_bins)
    dT = np.zeros(N_bins)
    Npairs = np.zeros(N_bins, dtype=np.int64)
    return dTw, w2, w, dT, Npairs


def pairwise_full_matrix_singlecore(Dc, ra_rad, dec_rad,
                                    Tmapsc,
                                    r_max, r_step):
    group_limits = gen_groups_start_end_pw_fromrowtorow(len(Dc), 400)
    dTws, w2s, ws, dTs, Npairss = [], [], [], [], []
    for group in group_limits:
        start, end = group[0], group[1]
        blank_output = gen_output_pairwise_one_row(r_max,
                                                   r_step)
        dTw = blank_output[0]
        w2 = blank_output[1]
        w = blank_output[2]
        dT = blank_output[3]
        Npairs = blank_output[4]
        
        pairwise_from_rowtorow(start, end,
                                  Dc, ra_rad, dec_rad,
                                  Tmapsc, r_max, r_step,
                                  dTw, w2, w, dT, Npairs)
        dTws.append(dTw)
        w2s.append(w2)
        ws.append(w)
        dTs.append(dT)
        Npairss.append(Npairs)

    dTw = np.sum(np.array(dTws), axis=0)
    w2 = np.sum(np.array(w2s), axis=0)
    w = np.sum(np.array(ws), axis=0)
    dT = np.sum(np.array(dTs), axis=0)
    Npairs = np.sum(np.array(Npairss), axis=0)
    return dTw, w2, w, dT, Npairs


def gen_output_pairwise_ksz(r_max, r_step, Ngroups):
    N_bins = get_nbins(r_max, r_step)

    dTw = [np.zeros(N_bins, dtype=np.float64) for j in range(Ngroups)]
    w2 = [np.zeros(N_bins, dtype=np.float64) for j in range(Ngroups)]
    w = [np.zeros(N_bins, dtype=np.float64) for j in range(Ngroups)]
    dT = [np.zeros(N_bins, dtype=np.float64) for j in range(Ngroups)]
    Npairs = [np.zeros(N_bins, dtype=np.int64) for j in range(Ngroups)]
    return dTw, w2, w, dT, Npairs


def pairwise_ksz(Dc, ra_deg, dec_deg, Tmapsc, r_max, r_step, 
                 Nthreads=1, Ngroups=400):
    assert len(Dc) == len(ra_deg)
    assert len(Dc) == len(dec_deg)
    assert len(Dc) == len(Tmapsc)

    Ngal = len(Dc)
    ra_rad = np.deg2rad(ra_deg)
    dec_rad = np.deg2rad(dec_deg)

    group_limits = gen_groups_start_end_pw_fromrowtorow(Ngal, Ngroups)

    group_starts = [group[0] for group in group_limits]
    group_ends = [group[1] for group in group_limits]

    dTw, w2, w, dT, Npairs = gen_output_pairwise_ksz(r_max,
                                           r_step, Ngroups)

    Dcs = itertools.repeat(Dc, Ngroups)
    ras_rad = itertools.repeat(ra_rad, Ngroups)
    decs_rad = itertools.repeat(dec_rad, Ngroups)
    Tmapscs = itertools.repeat(Tmapsc, Ngroups)
    r_maxs = itertools.repeat(r_max, Ngroups)
    r_steps = itertools.repeat(r_step, Ngroups)

    with futures.ThreadPoolExecutor(Nthreads) as ex:
        ex.map(pairwise_from_rowtorow,
               group_starts, group_ends,
               Dcs, ras_rad, decs_rad, Tmapscs, r_maxs, r_steps,
               dTw, w2, w, dT, Npairs)
    dTw = np.sum(dTw, axis=0)
    w2 = np.sum(w2, axis=0)
    w = np.sum(w, axis=0)
    dT = np.sum(dT, axis=0)
    Npairs = np.sum(Npairs, axis=0)

    return dTw, w2, w, dT, Npairs
