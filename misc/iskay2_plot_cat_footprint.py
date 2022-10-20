#! /usr/bin/env python
from iskay2 import tile_tools
import healpy as hp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
from iskay2.paramfile import load_paramfile

show = False
param_fname = sys.argv[1]
NSIDE = 64
params = load_paramfile(fname=param_fname)

print("Approximate resolution at NSIDE {} is {:.2} deg".format(
      NSIDE, hp.nside2resol(NSIDE, arcmin=True) / 60))
res_el = hp.nside2resol(NSIDE, arcmin=True) / 60


df = pd.read_hdf('./ApPhotoResults/ap_photo.hdf')

df = tile_tools.classify_grid(df, Nside=NSIDE)
hist = tile_tools.healpix_histogram_catalog(df, NSIDE)
number_of_nonzero_pixs = np.sum(hist > 0)

area_el = res_el**2
cat_area = area_el * number_of_nonzero_pixs

print("number of non-zero pixs: %i" % number_of_nonzero_pixs)
print("covered_area is: %1.2f sq deg" % (cat_area))

DIRNAME = './cat_plots'
if not os.path.exists(DIRNAME):
    os.mkdir(DIRNAME)
plot_fname = os.path.join(DIRNAME, params['NAME'])

hp.mollview(hist)
plt.title("Area: %1.2f sq. deg" % cat_area)

if show:
    plt.show()
else:
    plt.savefig(plot_fname + ".pdf")
    plt.savefig(plot_fname + ".png", dpi=120)
    plt.close()
