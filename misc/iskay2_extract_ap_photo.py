#! /usr/bin/env python
from iskay2 import paramfile
from iskay2 import maptools
from iskay2 import rcfile
from iskay2 import catalogtools
from iskay2 import ap_photo
import pandas as pd

params = paramfile.load_paramfile()
rc = rcfile.load_rcfile()

themap = maptools.load_map(params, rc=rc)
the_noise_map = maptools.load_map(params,
                                  rc=rc,
                                  noisemap=True)
df = catalogtools.load_catalog(params, rc=rc)
#df = df.head(1000)

df_ap_photo = ap_photo.get_ap_photo_in_catalog(df,
                                               themap,
                                               the_noise_map,
                                               params, rc)

df_masks = ap_photo.get_mask_tags(df, params, rc)
df_tosave = pd.concat([df_ap_photo, df_masks], 
                      axis='columns')

ap_photo.save_ap_photo(df_tosave)
