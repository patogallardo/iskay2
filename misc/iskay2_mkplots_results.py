import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import glob
import json

font = {'size'   : 18}

matplotlib.rc('font', **font)

fnames = glob.glob('./params_*.json')
for j in range(len(fnames)):
    print(fnames[j])
    with open(fnames[j], 'r') as f:
        paramfile = json.load(f)
    resultfilename = './results/' + paramfile['NAME'] + '.hdf'

    df_curve = pd.read_hdf(resultfilename, 'df_ksz_err')
    df_cov = pd.read_hdf(resultfilename, 'df_cov')


    plt.figure(figsize=[6, 3])
    plt.axhline(0, color='black')
    plt.errorbar(df_curve.r_mp,
             df_curve.ksz_curve,
             yerr=df_curve.errorbar,
             capsize=5,
             marker='.',
             ls='')

    plt.xlabel('r [Mpc]')
    plt.ylabel('$p_{est}$')
    plt.title(paramfile['NAME'], size=9)
    plt.ylim([-0.15, 0.02])
    plt.tight_layout()
    plt.show()
