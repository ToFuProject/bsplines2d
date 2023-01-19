# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 20:14:40 2023

@author: dvezinet
"""


import numpy as np
import scipy.stats as scpst
import astropy.units as asunits
import datastock as ds


# ############################################################
# ############################################################
#               interpolate spectral
# ############################################################


def binning(
    coll=None,
    key=None,
    ref_key=None,
    bins=None,
):
    """ Return the spectrally interpolated coefs
    
    Either E xor Ebins can be provided
    - E: return interpolated coefs
    - Ebins: return binned (integrated) coefs
    """
    
    # ----------
    # checks
    
    # keys
    isbs, key, ref_key, axis = _binning_check_keys(
        coll=coll,
        key=key,
        ref_key=ref_key,
    )
        
    # ----------
    # trivial
    
    if not isbs:
        return ds._class1_binning.binning(
            coll=coll,
            key=key,
            ref_key=ref_key,
            bins=bins,
        )
    
    # ---------
    # keys
    
    # units
    keym = coll.dobj[coll._which_bsplines][ref_key][coll._which_mesh]
    kknots = coll.dobj[coll._which_mesh][keym]['knots'][0]
    
    units = coll.ddata[key]['units']
    units_ref = coll.ddata[kknots]['units']
    
    # vect
    vect = coll.ddata[kknots]['data']
    
    # ----------
    # bins
    
    # bins
    bins, units_bins, db, npts = _binning_check_bins(
        coll=coll,
        bins=bins,
        vect=vect,
        deg=coll.dobj[coll._which_bsplines][ref_key]['deg'],
    )
    
    # units
    units = ds._class1_binning._units(
        key=key,
        units=units,
        units_ref=units_ref,
        units_bins=units_bins,
    )
        
    # ------------
    # bsplines
    
    coefs = coll.ddata[key]['data']
    clas = coll.dobj[coll._which_bsplines][ref_key]['class']

    # sample mesh
    res0 = np.abs(np.min(np.diff(vect)))
    xx = coll.get_sample_mesh(keym, res=res0 / npts, mode='abs')[0]

    # interpolate
    val = clas(
        coefs=coefs,
        xx=xx,
        axis=axis,
        val_out=0.,
    )
    
    # bin
    val = ds._class1_binning._bin(
        bins=bins,
        db=db,
        vect=xx,
        data=val,
        axis=axis[0],
    )
    
    return val, units
    
    
# ###################################
#       check
# ####################################


def _binning_check_keys(
    coll=None,
    key=None,
    ref_key=None,
):

    # ---------
    # keys
    
    # dict of bsplines data
    dbs = coll.get_dict_bsplines()[0]
    
    # binning only for non-bsplines or 1d bsplines
    lok_nobs = [k0 for k0, v0 in coll.ddata.items() if k0 not in dbs.keys()]
    lok_bs = [
        k0 for k0, v0 in coll.ddata.items()
        if (
            k0 in dbs.keys()
            and any([len(v1) == 1 for v1 in dbs[k0].values()])
        )
    ]
    
    # ----------
    # key
    
    key = ds._generic_check._check_var(
        key, 'key',
        types=str,
        allowed=lok_nobs + lok_bs,
    )

    isbs = key in lok_bs

    # ------------
    # ref_key

    axis = None
    if isbs:
        
        lax = np.concatenate(tuple([v0 for v0 in dbs[key].values()]))
        lok_ref = [
            k0 for k0 in coll.ddata[key]['ref']
            if coll.ddata[key]['ref'].index(k0) not in lax
        ]
        lok_bs = [k0 for k0, v0 in dbs[key].items() if len(v0) == 1]
        
        ref_key = ds._generic_check._check_var(
            ref_key, 'ref_key',
            types=str,
            allowed=lok_ref + lok_bs,
        )

        isbs = ref_key in lok_bs
        
        if isbs:
            axis = dbs[key][ref_key]

    return isbs, key, ref_key, axis


def _binning_check_bins(
    coll=None,
    bins=None,
    vect=None,
    deg=None,
):
    
    # ---------
    # bins
    
    if isinstance(bins, str):
        lok = [k0 for k0, v0 in coll.ddata.items() if v0['monot'] == (True,)]
        bins = ds._generic_check._check_var(
            bins, 'bins',
            types=str,
            allowed=lok,
        )

        bins_units = coll.ddata[bins]['units']        
        bins = coll.ddata[bins]['data']
        
    else:
        bins_units = None
    
    bins = ds._generic_check._check_flat1darray(
        bins, 'bins',
        dtype=float,
        unique=True,
        can_be_None=False,
    )
    
    # -----------------
    # check uniformity
    
    db = np.abs(np.diff(bins))
    if not np.allclose(db[0], db):
        msg = (
            "Arg bins must be a uniform bin edges vector!"
            f"Provided diff(bins) = {db}"
        )
        raise Exception(msg)
    db = db[0]
        
    # ----------
    # bin method
    
    dv = np.abs(np.mean(np.diff(vect)))
    npts = (deg + 3) * max(1, dv / db)
        
    return bins, bins_units, db, npts