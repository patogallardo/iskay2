from iskay2 import pwksz
import progressbar
import pandas as pd

def get_bs(df, params, rc):
    '''Compute bootstrap iterating on 
    random catalogs.'''
    Nit = params['N_BOOTSTRAP_ITERATIONS']
    
    df_pws = []
    
    for j in progressbar.progressbar(range(Nit)):
        df_sampled = df.sample(n=len(df),
                               replace=True,
                               weights=None,
                               random_state = j)
        df_pw = pwksz.pw_compute_ksz(df_sampled,
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
        cov = np.cov(curves)
        corr = np.corrcoef(curves)
        
