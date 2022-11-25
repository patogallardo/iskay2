'''
Aperture photometry using pixell.

###
by: P. Gallardo
'''
from pixell import enmap
import numpy as np
from pixell import reproject
from multiprocessing import Pool
from multiprocessing import cpu_count
from itertools import repeat
import pandas as pd
from . import cosmology
from . import paramfile
from . import maptools
import os


def get_reprojection(themap, decra_rad, R_RAD_SUBMAP, OVERSAMPLE):
    '''Gets a submap from themap.
    themap: complete map
    decra_rad: declination, right ascencion in radians
    R_RAD_SUBMAP: radius of the postage stamp in radians
    OVERSAMPLE: typically 2 or None for speed.
    '''
    if OVERSAMPLE is None:
        res = None
    else:
        res = min(np.abs(themap.wcs.wcs.cdelt))*enmap.utils.degree/(2*OVERSAMPLE)
    submap = reproject.thumbnails(themap, decra_rad,r=R_RAD_SUBMAP, res=res)
    return submap


def get_ap_photo(r_disks, submap, R_RING_OVER_R_DISK=np.sqrt(2.)):
    '''For one submap, get aperture photometry for a set of
    r_disks.
    r_disks: radius in arcmins for the aperture in the disk
    submap: the postage stamp to where to apply the aperture photo.
    R_RING_OVER_R_DISK: the ratio r_ring/r_disk = 1.4142... for
    equal area.
    '''
    r_rings = r_disks * R_RING_OVER_R_DISK
    r_arcmin_map = np.rad2deg(submap.modrmap()) * 60

    t_rings, t_disks = np.zeros_like(r_disks), np.zeros_like(r_disks)

    for j, (r_disk, r_ring) in enumerate(zip(r_disks, r_rings)):
        sel_disk = r_arcmin_map < r_disk
        sel_ring = np.logical_and(r_arcmin_map >= r_disk, r_arcmin_map <= r_ring) 
        t_disks[j] = np.mean(submap[sel_disk])
        t_rings[j] = np.mean(submap[sel_ring])
    return t_disks, t_rings


def get_ap_photo_in_one_coordinate(r_disks_arcmin_decra_rad):
    '''This needs a variable called "themap" of type global.
    Bad practice, but it is necesary to use subprocess on shared
    memory.'''
    r_disks_arcmin, decra_rad = r_disks_arcmin_decra_rad
    submap = get_reprojection(themap, decra_rad, R_RAD_SUBMAP, OVERSAMPLE)
    
    ap_photo_res = get_ap_photo(r_disks_arcmin, submap, R_RING_OVER_R_DISK=R_RING_OVER_R_DISK)
    #R_RING_OVER_R_DISK has to exist in a higher scope.
    return ap_photo_res


def get_submap_native_pix(themap, decra_rad, R_RAD_SUBMAP):
    '''Receives themap, dec_ra in radians and a semiwidth.
    returns a square pixel
    '''
    dec_rad, ra_rad = decra_rad
    box = np.array([[dec_rad - R_RAD_SUBMAP, ra_rad - R_RAD_SUBMAP],
                    [dec_rad + R_RAD_SUBMAP, ra_rad + R_RAD_SUBMAP]])
    submap = themap.submap(box)
    return submap


def is_submap_masked(r_disk_arcmin, submap):
    '''Check if submap has True values within r_disk_arcmin'''
    r_arcmin_map = np.rad2deg(submap.modrmap()) * 60
    sel_disk = r_arcmin_map < r_disk_arcmin
    is_true = np.any(submap[sel_disk]==False) # if any is zero, where False is do not use
    return is_true


def get_mask_status_in_one_coordinate(r_disk_arcmin_decra_rad):
    '''needs a global variable called themap
    '''
    r_disk_arcmin, decra_rad = r_disk_arcmin_decra_rad
    submap = get_submap_native_pix(themap, decra_rad, R_RAD_SUBMAP)
    is_masked = is_submap_masked(r_disk_arcmin, submap)
    return is_masked


def get_ap_mask_full_cat(ras_deg, decs_deg,
                         the_map, r_disk_arcmin,
                         r_rad_submap=np.deg2rad(10./60.),
                         Nproc=2,
                         mask_name="masked"): 
    '''Determines if galaxies are masked for the full catalog.
    Identical params than get_ap_photo_full_cat
    '''
    global themap
    global R_RAD_SUBMAP
    themap = the_map
    R_RAD_SUBMAP = r_rad_submap
    ras_rad = np.deg2rad(ras_deg)
    decs_rad = np.deg2rad(decs_deg)
    coords_decs_ras_rad = np.vstack([decs_rad, ras_rad]).T

    args = zip(repeat(r_disk_arcmin), coords_decs_ras_rad)
    
    with Pool(processes=Nproc) as pool:
        res = pool.map(get_mask_status_in_one_coordinate, args)
    masked = np.array(res)

    vals = np.hstack([masked])
    title = [mask_name]
    df = pd.DataFrame(vals,
                      columns=title)
    return df


   

def get_ap_photo_full_cat(ras_deg, decs_deg,
                          the_map,
                          r_disks_arcmin,
                          
                          r_rad_submap=np.deg2rad(10./60.),
                          oversample=2,
                          r_ring_over_r_disk=np.sqrt(2),
                          Nproc=2,
                          noisemap=False):
    '''Receives a list of ra, dec in degrees, the map and
    a list of radii for the disks and ratios of disk and ring.
    ras, decs_deg: ra, dec in degrees
    themap: the map!
    r_disks_arcmin: list of r_disks in arcmin
    R_RING_OVER_R_DISK: ratio of r ring over r disk.
    '''
    if noisemap:
        r_disks_arcmin = np.array([2.1])
    else:
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
        res = pool.map(get_ap_photo_in_one_coordinate, args)

    T_disks, T_rings = np.array(res)[:, 0], np.array(res)[:, 1]
    dTs = T_disks - T_rings

    r_rings_arcmin = r_disks_arcmin * R_RING_OVER_R_DISK

    if noisemap:
        titles_disks = ['noise_uk']
    else:
        titles_disks = [("T_disk_%1.2f_arcmin" % r).replace('.', 'p') for r in r_disks_arcmin]
        titles_rings = [("T_ring_%1.2f_arcmin" % r).replace('.', 'p') for r in r_disks_arcmin]
        titles_dTs = [("dT_%1.2f_arcmin" % r).replace('.', 'p') for r in r_disks_arcmin]

    if noisemap:
        T_disks = 1/np.sqrt(T_disks)
        vals = np.hstack([T_disks])
        df = pd.DataFrame(vals,
                          columns=titles_disks)
    else:
        vals = np.hstack([T_disks, T_rings, dTs])
        df = pd.DataFrame(vals,
               columns=titles_disks + titles_rings + titles_dTs)

    return df


def save_ap_photo(df_res):
    if not os.path.exists("./ApPhotoResults"):
        os.mkdir("ApPhotoResults")
    df_res.to_hdf("ApPhotoResults/ap_photo.hdf",
                  key='df_ap_photo',
                  mode='w')


def get_mask_tags(df_cat, params, rc):
    mask_fnames, mask_names = paramfile.get_mask_fnames(params)
    ras_deg, decs_deg = df_cat.ra.values, df_cat.dec.values
    r_disk_arcmin = params["MASKED_RADIUS_ARCMIN"]
    Nproc = rc["NPROC_AP_PHOTO"]
    dfs = []
    for j, fname in enumerate(mask_fnames):
        the_mask = maptools.load_mask(fname, rc)
        df_out = get_ap_mask_full_cat(ras_deg, decs_deg,
                         the_mask, r_disk_arcmin,
                         r_rad_submap=rc["R_RAD_SUBMAP"],
                         Nproc=Nproc,
                         mask_name=mask_names[j])
        dfs.append(df_out)
    dfs = pd.concat(dfs, axis='columns')
    return dfs


def get_ap_photo_in_catalog(df_cat, themap,
    themap_noise,
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

    df_ap_photo = get_ap_photo_full_cat(ras_deg, decs_deg,
                                        themap,
                                        r_disks_arcmin,
                    r_rad_submap=r_rad_submap,
                    oversample=oversample,
                    r_ring_over_r_disk=r_ring_over_r_disk,
                    Nproc=Nproc)
    df_noise = get_ap_photo_full_cat(ras_deg, decs_deg,
                                     themap_noise,
                                     r_disks_arcmin,
                    r_rad_submap=r_rad_submap,
                    oversample=oversample,
                    r_ring_over_r_disk=r_ring_over_r_disk,
                    Nproc=Nproc,
                    noisemap=True)
    # define cosmology and compute distances
    z = df_cat.z.values
    c = cosmology.cosmo(Omega_m = omega_m, Little_h = little_h)
    d_mpc = cosmology.Dc(z, c, distance='Mpc')
    d_mpc_over_h = cosmology.Dc(z, c, distance='Mpc_over_little_h')
    # end cosmology distance calculation
    df_cat['d_mpc'] = d_mpc
    df_cat['d_mpc_over_h'] = d_mpc_over_h
    df_res = pd.concat([df_cat, df_ap_photo, df_noise],
                       axis='columns')
    return df_res
