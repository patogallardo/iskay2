'''
Aperture photometry legacy code from v1.
'''
import numpy as np
import pandas as pd

from pixell import enmap
from pixell import reproject
from pixell import utils
from multiprocessing import Pool
from multiprocessing import cpu_count
from itertools import repeat
from . import cosmology
from . import ap_photo
import os


def get_legacy_ap_photo_in_catalog_and_save(df_cat, themap,
    params, rc):
    '''Computes ap_photo for all rows in df_cat. Needs to
    suppply themap and, rc and params files.
    df_cat: dataframe containing ra, dec in degrees as a col.
    themap: map to apply aperture photo on.
    params: param file with parameters for this run.
    rc: rc file with values for this run.
    '''
    ras_deg = df_cat.ra.values
    decs_deg = df_cat.dec.values
    r_disks_arcmin = params["R_DISKS_ARCMIN"]
    r_rad_submap = rc["R_RAD_SUBMAP"]
    oversample = rc["OVERSAMPLE"]
    r_ring_over_r_disk = params["R_RING_OVER_R_DISK"]
    Nproc = rc["NPROC_AP_PHOTO"]

    omega_m = params["OMEGA_M"]
    little_h = params["LITTLE_H"]

    df_ap_photo = get_ap_photo_legacy_full_cat(ras_deg, decs_deg,
                                               themap,
                                               r_disks_arcmin,
                    r_rad_submap=r_rad_submap,
                    oversample=oversample,
                    r_ring_over_r_disk=r_ring_over_r_disk,
                    Nproc=Nproc)
    # define cosmology and compute distances
    z = df_cat.z.values
    c = cosmology.cosmo(Omega_m = omega_m, Little_h = little_h)
    d_mpc = cosmology.Dc(z, c, distance='Mpc')
    d_mpc_over_h = cosmology.Dc(z, c, distance='Mpc_over_little_h')
    # end cosmology distance calculation
    df_cat['d_mpc'] = d_mpc
    df_cat['d_mpc_over_h'] = d_mpc_over_h
    save_ap_photo_legacy(df_cat, df_ap_photo)


def save_ap_photo_legacy(df_cat, df_ap_photo):
    if not os.path.exists("./ApPhotoResults"):
        os.mkdir("ApPhotoResults")
    df_ap_photo.index = df_cat.index
    df = pd.concat([df_cat, df_ap_photo],
                   axis='columns')
    df.to_hdf("ApPhotoResults/ap_photo_legacy.hdf",
              key='df_ap_photo',
              mode='w')


def get_ap_photo_legacy_full_cat(ras_deg, decs_deg,
                                 the_map,
                                 r_disks_arcmin,
                                 r_rad_submap=np.deg2rad(10./60.),
                                 oversample=2,
                                 r_ring_over_r_disk=np.sqrt(2),
                                 Nproc=2):
    '''Receives a list of ra, dec in degrees, the map and
    a list of radii for the disks and ratios of disk and ring.
    ras, decs_deg: ra, dec in degrees
    themap: the map!
    r_disks_arcmin: list of r_disks in arcmin
    R_RING_OVER_R_DISK: ratio of r ring over r disk.
    '''
    r_disks_arcmin = np.array(r_disks_arcmin)
    global themap # use this to use themap within get_ap_photo_in_one_coordinate
    global R_RAD_SUBMAP
    global OVERSAMPLE
    global R_RING_OVER_R_DISK
    
    themap = the_map
    R_RAD_SUBMAP = r_rad_submap
    R_RING_OVER_R_DISK = r_ring_over_r_disk
    OVERSAMPLE = oversample
    
    ras_rad = np.deg2rad(ras_deg)
    decs_rad = np.deg2rad(decs_deg)
    coords_decs_ras_rad = np.vstack([decs_rad, ras_rad]).T

    args = zip(repeat(r_disks_arcmin), coords_decs_ras_rad)

    with Pool(processes=Nproc) as pool:
        res = pool.map(get_ap_photo_legacy_in_one_coordinate, args)

#    res = list(map(get_ap_photo_legacy_in_one_coordinate, args))

#    res = []
#    for j in range(len(coords_decs_ras_rad)):
#        print(j)
#        res.append(get_ap_photo_legacy_in_one_coordinate([r_disks_arcmin, coords_decs_ras_rad[j]]))

    T_disks, T_rings = np.array(res)[:, 0], np.array(res)[:, 1]
    dTs = T_disks - T_rings

    r_rings_arcmin = r_disks_arcmin * R_RING_OVER_R_DISK
    titles_disks = ["T_disk_%1.2f_arcmin" % r for r in r_disks_arcmin]
    titles_rings = ["T_ring_%1.2f_arcmin" % r for r in r_disks_arcmin]
    titles_dTs = ["dT_%1.2f_arcmin" % r for r in r_disks_arcmin]
    
    vals = np.hstack([T_disks, T_rings, dTs])
    df = pd.DataFrame(vals,
               columns=titles_disks + titles_rings + titles_dTs)

    return df


def gen_box(ra_rad, dec_rad, semiwidth_deg):
    semiwidth_rad = np.deg2rad(semiwidth_deg)
    box = np.empty([2, 2])
    box[0, 0] = dec_rad - semiwidth_rad
    box[0, 1] = ra_rad - semiwidth_rad
    box[1, 0] = dec_rad + semiwidth_rad
    box[1, 1] = ra_rad + semiwidth_rad
    return box


def get_ap_photo_legacy_in_one_coordinate(r_disks_arcmin_decra_rad):
    '''This needs a variable called "themap" of type global.
    Bad practice, but it is necesary to use subprocess on shared
    memory.'''
    r_disks_arcmin, decra_rad = r_disks_arcmin_decra_rad
    semi_width_deg_for_submap = np.rad2deg(R_RAD_SUBMAP) # degree
    box = gen_box(decra_rad[1], decra_rad[0], semi_width_deg_for_submap)
    submap = enmap.submap(themap, box) # original pixelization
    
    submap_repix = enmap.resample(submap,
                                  10*np.array(submap.shape))
    OVERSAMPLE = 0.5 # this will not oversample and inherit sampling from previous line
    submap_reprojected = ap_photo.get_reprojection(submap_repix,
                                                   decra_rad,
                                                   R_RAD_SUBMAP,
                                                   OVERSAMPLE)
    ap_photo_res = ap_photo.get_ap_photo(r_disks_arcmin, submap_reprojected,
                                  R_RING_OVER_R_DISK=R_RING_OVER_R_DISK)
    #R_RING_OVER_R_DISK has to exist in a higher scope.
    return ap_photo_res
