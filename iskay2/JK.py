from iskay2 import pwksz
import progressbar
import pandas as pd
import numpy as np
import os
import time

def get_bs(df, params, rc, seed=1):
    '''Compute bootstrap iterating on 
    random catalogs.'''
    Nit = params['N_BOOTSTRAP_ITERATIONS']
    ap_photo_sel_string = ("dT_%1.2f_arcmin" % params["R_DISK_ARCMIN_PAIRWISEKSZ"]).replace('.', 'p')
    df_pws = []
    
    df_bs = df.copy()
    dTs = df_bs[ap_photo_sel_string].values
    Nel = len(df)
    np.random.seed(seed)
    for j in progressbar.progressbar(range(Nit)):
        choose = np.random.choice(Nel, Nel)
        df_bs[ap_photo_sel_string] = dTs[choose]

        df_pw = pwksz.pw_compute_ksz(df_bs,
                                     params, rc)
        df_pw['realization'] = j
        df_pws.append(df_pw)
#    df_pws = pd.concat(df_pws)
    return df_pws


class BS:
    def __init__(self, df, params, rc, save=False):
        df_pw_full_dataset = pwksz.pw_compute_ksz(df, params, rc)
        t1 = time.time()
        df_pws = get_bs(df, params, rc)
        t2 = time.time()
        curves = np.array([df['ksz_curve'].values
                          for df in df_pws])
        errorbar_std = curves.std(axis=0)
        errorbar_percentiles = (np.percentile(curves, 50+68/2, axis=0) - 
                                np.percentile(curves, 50-68/2, axis=0))/2
        cov = np.cov(curves.T)
        corr = np.corrcoef(curves.T)
        
        self.r_mp = df_pw_full_dataset.r_mp.values
        self.r_mp_over_h = df_pw_full_dataset.r_mp_over_h.values
        self.ksz_curve_full_dataset = df_pw_full_dataset.ksz_curve.values
        self.ksz_curve_full_dataset_object = df_pw_full_dataset
        self.curves = curves
        self.df_pws = df_pws
        self.ksz_curves = curves
        self.errorbar_std = errorbar_std
        self.errorbar_percentiles = errorbar_percentiles
        self.cov = cov
        self.corr = corr

        if save:
            SAVEDIR = "./results/"
            if not os.path.exists(SAVEDIR):
                os.mkdir(SAVEDIR)
            hdf_out_fname = os.path.join(SAVEDIR, params["NAME"] + '.hdf')
            df_curve_err = pd.DataFrame({'r_mp': df_pw_full_dataset.r_mp.values ,
                                         'ksz_curve': df_pw_full_dataset.ksz_curve,
                                         'errorbar': errorbar_std})
            df_curve_err.to_hdf(hdf_out_fname, key='df_ksz_err')

            df_pw_full_dataset.to_hdf(hdf_out_fname,
                                      key='df_pw')

            labels = ["%i" % val for val in  np.round(df_pw_full_dataset.r_mp.values, 0)]

            df_corr = pd.DataFrame(corr, index=labels, columns=labels)
            df_corr.to_hdf(hdf_out_fname,
                           key='df_corr')

            df_cov = pd.DataFrame(cov, index=labels, columns=labels)
            df_cov.to_hdf(hdf_out_fname,
                          key='df_cov')
            s = pd.Series({'Ngal': len(df), 
                           'runtime': t2-t1,
           'N_BOOTSTRAP_ITERATIONS': params["N_BOOTSTRAP_ITERATIONS"]})
            s.to_hdf(hdf_out_fname, 
                     key='s_additional')
            df_curves = pd.DataFrame(curves)
            df_curves.to_hdf(hdf_out_fname, key='bs_curves')
