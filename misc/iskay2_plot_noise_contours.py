#! /usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d
import numpy as np
import os

SHOW = False
PLOTDIR = "cat_plots"
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

df = pd.read_hdf('./ApPhotoResults/ap_photo.hdf')

res = binned_statistic_2d(df.dec,
                          df.ra,
                          df.noise_uk, bins=100)

hist_2d = res[0]
y_edge = res[1]
y_edge = 0.5 * (y_edge[1:] + y_edge[:-1])
x_edge = res[2]
x_edge = 0.5 * (x_edge[1:] + x_edge[:-1])
hist_2d[np.isnan(hist_2d)] = 500


plt.figure(figsize=[8, 4.5])
plt.contourf(x_edge, y_edge,
             hist_2d,
             levels=[0, 25, 45, 75, 100])
plt.colorbar(label='noise [uK]')
plt.xlabel('ra [deg]')
plt.ylabel('dec [deg]')
if SHOW:
    plt.show()
else:
    plt.savefig(os.path.join(PLOTDIR, "noise_map.png"), dpi=120)
    plt.savefig(os.path.join(PLOTDIR, "noise_map.pdf"))
    plt.close()

nbins = 1000
fig, ax = plt.subplots(figsize=[8, 4.5])
ax.hist(df.noise_uk, histtype='step',
        color='black',
        range=[0, 100],
        bins=nbins,
        label='histogram')
ax.set_xlabel('noise [uK]')
ax.set_xticks(np.arange(0, 110, 10))
ax.grid(axis='x')

ax2 = ax.twinx()
ax2.hist(df.noise_uk, range=[0, 100], histtype='step',
         cumulative=True, density=True,
         bins=nbins,
         label='cumulative')
ax2.set_yticks(np.arange(0, 1.1, .10))
ax2.grid()

if SHOW:
    plt.show()
else:
    plt.savefig(os.path.join(PLOTDIR, "noise_hist.png"), dpi=120)
    plt.savefig(os.path.join(PLOTDIR, 'noise_hist.pdf'))
    plt.close()
