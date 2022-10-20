from iskay2 import paramfile
from iskay2 import maptools
from iskay2 import rcfile
from iskay2 import catalogtools
from iskay2 import ap_photo

params = paramfile.load_paramfile()
rc = rcfile.load_rcfile()

themap = maptools.load_map(params, rc=rc)
the_noise_map = maptools.load_map(params,
                                  rc=rc,
                                  noisemap=True)
df = catalogtools.load_catalog(params, rc=rc)
# df = df.head(1000)

ap_photo.get_ap_photo_in_catalog_and_save(df,
                                          themap,
                                          the_noise_map,
                                          params, rc)
