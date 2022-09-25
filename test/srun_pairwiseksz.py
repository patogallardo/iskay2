from iskay2.new_pairwiser import pairwise_ksz
import numpy as np
import os

Ngal = 300000

np.random.seed(1)
Dc = np.random.uniform(0, 500, Ngal)
ra_deg = np.random.uniform(0, 180, Ngal)
dec_deg = np.random.uniform(0, 40, Ngal)
Tmapsc = np.random.normal(0, 1, Ngal)

r_max = 300
r_step = 10

Nthreads = int(os.environ.get("NCORES"))
Ngroups = 400

dTw, w2, w, dT, Npairs = pairwise_ksz(Dc, ra_deg, dec_deg,
                                      Tmapsc, r_max, r_step,
                                      Nthreads=Nthreads,
                                      Ngroups=400)
a = np.array([dTw, w2, w, dT, Npairs])
np.savetxt('/scratch/r/rbond/gallardo/test/a%02i.txt' % Nthreads,
 a)
