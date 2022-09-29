from iskay2 import pwksz
import progressbar
import pandas as pd
import numpy as np

def get_bs(df, params, rc):
    '''Compute bootstrap iterating on 
    random catalogs.'''
    Nit = params['N_BOOTSTRAP_ITERATIONS']
    ap_photo_sel_string = "dT_%1.2f_arcmin" % params["R_DISK_ARCMIN_PAIRWISEKSZ"]
    df_pws = []
    
    df_bs = df.copy()
    dTs = df_bs[ap_photo_sel_string].values
    Nel = len(df)
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
    def __init__(self, df, params, rc):
        df_pw_full_dataset = pwksz.pw_compute_ksz(df, params, rc)
        df_pws = get_bs(df, params, rc)
        curves = np.array([df['ksz_curve'].values
                          for df in df_pws])
        errorbar_std = curves.std(axis=0)
        errorbar_percentiles = (np.percentile(curves, 50+68/2, axis=0) - 
                                np.percentile(curves, 50-68/2, axis=0))/2
        cov = np.cov(curves.T)
        corr = np.corrcoef(curves.T)
        
        self.ksz_curve_full_dataset = df_pw_full_dataset
        self.curves = curves
        self.df_pws = df_pws
        self.ksz_curves = curves
        self.errorbar_std = errorbar_std
        self.errorbar_percentiles = errorbar_percentiles
        self.cov = cov
        self.corr = corr
