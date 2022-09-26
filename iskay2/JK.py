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
