from pixell import enmap
from os.path import join
from . import rcfile

def load_map(params, rc=None, noisemap=False):
    '''Opens map in fname.
    Uses default directory in rc file.'''
    if noisemap:
        fname = params["DIVMAP_FNAME"]
    else:
        fname = params["MAP_FNAME"]

    if fname == "None":
        return None

    print("Loading map: %s" % fname)
    if rcfile is None:
        rc = rcfile.load_rcfile()
    map_path = rc['MAPS_PATH']
    map_path = join(map_path, fname)
    themap = enmap.read_fits(map_path)
    if themap.shape[0] == 3:
        themap = themap[0]
    return themap
