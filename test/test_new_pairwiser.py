from iskay2.new_pairwiser import pairwise_one_row
from iskay2.new_pairwiser import pairwise_from_rowtorow_onlyonerowatatime
from iskay2.new_pairwiser import pairwise_from_rowtorow
from iskay2.new_pairwiser import row_end
from iskay2.new_pairwiser import get_nbins
from iskay2.new_pairwiser import pairwise_full_matrix_singlecore
from iskay2.new_pairwiser import pairwise_ksz
import numpy as np
import math as mt

def setup_inputs(Ngal):
    Dc = np.random.uniform(10, 500, size=Ngal)
    ra_rad = np.random.uniform(0, np.pi/2, Ngal)
    dec_rad = np.random.uniform(0, np.pi/6, Ngal)
    Tmapsc = np.random.normal(0, 5, Ngal)
    return Dc, ra_rad, dec_rad, Tmapsc


def setup_outputs():
    r_max, r_step = 300, 10
    Nbins = int(300/10)
    dTw = np.zeros(Nbins)
    w2 = np.zeros(Nbins)
    w = np.zeros(Nbins)
    dT = np.zeros(Nbins)
    Npairs = np.zeros(Nbins, dtype=np.int64)
    return r_max, r_step, dTw, w2, w, dT, Npairs


def test_pairwise_one_row():
    row = 1
    Ngal = 1000
    Dc, ra_rad, dec_rad, Tmapsc = setup_inputs(Ngal)
    r_max, r_step, dTw, w2, w, dT, Npairs = setup_outputs()
    pairwise_one_row(row, Ngal, Dc, ra_rad, dec_rad,
                     Tmapsc, r_max, r_step, dTw, w2, w, dT, Npairs)


def test_pairwise_from_rowtorow():
    Ngal = 1000
    Dc, ra_rad, dec_rad, Tmapsc = setup_inputs(Ngal)
    
    r_max, r_step, dTw_1, w2_1, w_1, dT_1, Npairs_1 = setup_outputs()
    r_max, r_step, dTw_2, w2_2, w_2, dT_2, Npairs_2 = setup_outputs()
    pairwise_from_rowtorow_onlyonerowatatime(0, Ngal-1,
                                             Dc, ra_rad, dec_rad, Tmapsc,
                                             r_max, r_step, dTw_1, w2_1, w_1, dT_1, Npairs_1)
    rowend = row_end(Ngal)
    pairwise_from_rowtorow(0, rowend,
                           Dc, ra_rad, dec_rad, Tmapsc,
                           r_max, r_step,
                           dTw_2, w2_2, w_2, dT_2, Npairs_2)
    assert np.allclose(dTw_1, dTw_2)

    Ngal = 1001
    Dc, ra_rad, dec_rad, Tmapsc = setup_inputs(Ngal)
    
    r_max, r_step, dTw_1, w2_1, w_1, dT_1, Npairs_1 = setup_outputs()
    r_max, r_step, dTw_2, w2_2, w_2, dT_2, Npairs_2 = setup_outputs()
    pairwise_from_rowtorow_onlyonerowatatime(0, Ngal-1,
                                             Dc, ra_rad, dec_rad, Tmapsc,
                                             r_max, r_step, dTw_1, w2_1, w_1, dT_1, Npairs_1)
    rowend = row_end(Ngal) # will end at ceil(Ngal-1)/2-1
    pairwise_from_rowtorow(0, rowend,
                           Dc, ra_rad, dec_rad, Tmapsc,
                           r_max, r_step,
                           dTw_2, w2_2, w_2, dT_2, Npairs_2)
    assert np.allclose(dTw_1, dTw_2)


def test_row_end():
    Ngal = 100
    assert mt.ceil((Ngal-1)/2) - 1 == row_end(Ngal)


def test_getnbins():
    r_max = 300
    r_step = 10
    assert get_nbins(r_max, r_step) == 30

    r_max = 45
    assert get_nbins(r_max, r_step) == 5


def test_pw_fullmat_singlecore():
    Ngal = 5000
    Dc, ra_rad, dec_rad, Tmapsc = setup_inputs(Ngal)
    r_max, r_step, dTw_1, w2_1, w_1, dT_1, Npairs_1 = setup_outputs()
    
    pairwise_from_rowtorow_onlyonerowatatime(0, Ngal-1,
                                             Dc, ra_rad, dec_rad, Tmapsc,
                                             r_max, r_step, dTw_1, w2_1, w_1, dT_1, Npairs_1)
    dTw_2, w2_2, w_2, dT_2, Npairs_2 = pairwise_full_matrix_singlecore(Dc,
                                  ra_rad, dec_rad, 
                                  Tmapsc,
                                  r_max, r_step)
    assert np.allclose(Npairs_1, Npairs_2)


def test_pwksz():
    Ngal = 5000
    Dc, ra_rad, dec_rad, Tmapsc = setup_inputs(Ngal)
    ra_deg, dec_deg = np.rad2deg(ra_rad), np.rad2deg(dec_rad)
    r_max, r_step, dTw_1, w2_1, w_1, dT_1, Npairs_1 = setup_outputs()
    pairwise_from_rowtorow_onlyonerowatatime(0, Ngal-1,
                                             Dc, ra_rad, dec_rad, Tmapsc,
                                             r_max, r_step, dTw_1, w2_1, w_1, dT_1, Npairs_1)
    
    dTw, w2, w, dT, Npairs = pairwise_ksz(Dc, ra_deg, dec_deg,
                                          Tmapsc, r_max, r_step,
                                          Nthreads=2,
                                          Ngroups=400)
    assert np.allclose(dTw_1, dTw)
    assert np.allclose(w2_1, w2)
    assert np.allclose(w_1, w)
    assert np.allclose(dT_1, dT)
    assert np.allclose(Npairs_1, Npairs)

    Ngal = 5001
    Dc, ra_rad, dec_rad, Tmapsc = setup_inputs(Ngal)
    ra_deg, dec_deg = np.rad2deg(ra_rad), np.rad2deg(dec_rad)
    r_max, r_step, dTw_1, w2_1, w_1, dT_1, Npairs_1 = setup_outputs()
    pairwise_from_rowtorow_onlyonerowatatime(0, Ngal-1,
                                             Dc, ra_rad, dec_rad, Tmapsc,
                                             r_max, r_step, dTw_1, w2_1, w_1, dT_1, Npairs_1)
    
    dTw, w2, w, dT, Npairs = pairwise_ksz(Dc, ra_deg, dec_deg,
                                          Tmapsc, r_max, r_step,
                                          Nthreads=2,
                                          Ngroups=400)
    assert np.allclose(dTw_1, dTw)
    assert np.allclose(w2_1, w2)
    assert np.allclose(w_1, w)
    assert np.allclose(dT_1, dT)
    assert np.allclose(Npairs_1, Npairs)
