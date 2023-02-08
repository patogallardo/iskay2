#! /usr/bin/env python
'''Receive a parameter file and compute pairwise ksz.
Saves result in a hdf file and quits.
'''
from iskay2 import pwksz
from iskay2 import paramfile
from iskay2 import rcfile
from iskay2 import JK
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Runs pairwise ksz estimator on param file')
parser.add_argument('paramfname', type=str, help='The json filename to process: param.json')
args = parser.parse_args()


ap_photo_fname = "./ApPhotoResults/ap_photo.hdf"

params = paramfile.load_paramfile(args.paramfname)
rc = rcfile.load_rcfile()

df = pd.read_hdf(ap_photo_fname)

query = params["QUERY"]
df_queried = df.query(query, inplace=False)

print("N: %i" %len(df_queried))
# now compute it with the recipe, params, and rc files.
df_pw1 = pwksz.pw_compute_ksz(df_queried, params, rc)
bs1 = JK.BS(df_queried, params, rc, save=True)
