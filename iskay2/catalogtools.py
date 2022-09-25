import pandas as pd
from os.path import join
from . import rcfile

def load_catalog(params, rc=None):
    '''Opens catalog.
    Uses default directory in rc file.
    fname: filename.csv
    '''
    if rc is None:
        rc = rcfile.load_rcfile()
    cat_path = rc['CATS_PATH']

    cat_fname = params["CAT_FNAME"]
    print("Loading catalog: %s" %cat_fname)
    cat_fname = join(cat_path, cat_fname)
    df = pd.read_csv(cat_fname, comment="#")
    return df
