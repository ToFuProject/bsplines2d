# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 20:14:40 2023

@author: dvezinet
"""


import itertools as itt


import numpy as np
import datastock as ds


# ############################################################
# ############################################################
#               interpolate spectral
# ############################################################


def binning(
    coll=None,
    keys=None,
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
    isbs, keys, ref_key, daxis, dunits, units_ref = _binning_check_keys(
        coll=coll,
        keys=keys,
        ref_key=ref_key,
        only1d=True,
    )
    
    # because 1d only
    if not isbs:
        ref_key = ref_key[0]
        for k0, v0 in daxis.items():
            daxis[k0] = v0[0]
        units_ref = units_ref[0]
        
    # ----------
    # trivial
    
    if not isbs:
        return ds._class1_binning.binning(
            coll=coll,
            keys=keys,
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
    keys=None,
    ref_key=None,
    only1d=None,
):
    
    # only1d
    only1d = ds._generic_check._check_var(
        only1d, 'only1d',
        types=bool,
        default=True,
    )
    
    maxd = 1 if only1d else 2

    # ---------------
    # keys vs ref_key
    
    # ref_key
    dbs = coll.get_dict_bsplines()[0]
    if ref_key is not None:
        
        # basic checks
        if isinstance(ref_key, str):
            ref_key = (ref_key,)
            
        lref = list(coll.dref.keys())
        ldata = list(coll.ddata.keys())
        lbs = list(coll.dobj[coll._which_bsplines].keys())
        
        ref_key = list(ds._generic_check._check_var_iter(
            ref_key, 'ref_key',
            types=(list, tuple),
            types_iter=str,
            allowed=lref + ldata + lbs,
        ))
        
        lc = [
            all([rr in lref + ldata for rr in ref_key]),
            all([rr in lbs for rr in ref_key]),
        ]
        if np.sum(lc) != 1:
            msg = (
                "Arg ref_key must refer to (ref or vector) xor bsplines!\n"
                f"Provided: {ref_key}"
            )
            raise Exception(msg)
        
        # check vs maxd
        if (lc[0] and len(ref_key) > maxd) or (lc[1] and len(ref_key) != 1):
            msg = (
                "Arg ref_key shall not more than {maxd} ref!\n"
                "And can only contain one single bsplines!\n"
                f"Provided: {ref_key}"
            )
            raise Exception(msg)
        
        # bs vs refs
        if lc[1]:
            ref_key = ref_key[0]
            lok_bs = [
                k0 for k0, v0 in coll.ddata.items()
                if ref_key in v0[coll._which_bsplines]
            ]
            lok_nobs = []
            
        else:
        
            # check vs valid vectors
            for ii, rr in enumerate(ref_key):
                if rr in lref:
                    kwd = {'ref': rr}
                else:
                    kwd = {'key': rr}
                hasref, hasvect, ref, ref_key[ii] = coll.get_ref_vector(**kwd)[:4]
                
                if not (hasref and hasvect):
                    msg = (
                        f"Provided ref_key[{ii}] not a valid ref or ref vector!\n"
                        "Provided: {rr}"
                    )
                    raise Exception(msg)
                
        
            if not (hasref and hasvect):
                msg = (
                    "Provided ref_key not a valid ref or ref vector!\n"
                    "Provided: {ref_key}"
                )
                raise Exception(msg)
                
            lok_nobs = [
                k0 for k0, v0 in coll.ddata.items()
                if all([coll.ddata[rr]['ref'][0] in v0['ref'] for rr in ref_key])
            ]
            lok_bs = []
            
        if keys is None:
            keys = lok_nobs + lok_bs
                
    else:
        # binning only for non-bsplines or 1d bsplines
        lok_nobs = [k0 for k0, v0 in coll.ddata.items() if k0 not in dbs.keys()]
        lok_bs = [
            k0 for k0, v0 in coll.ddata.items()
            if (
                k0 in dbs.keys()
                and any([len(v1) <= maxd for v1 in dbs[k0].values()])
            )
        ]
    
    # ---------
    # keys
    
    if isinstance(keys, str):
        keys = [keys]

    keys = ds._generic_check._check_var_iter(
        keys, 'keys',
        types_iter=str,
        types=list,
        allowed=lok_nobs + lok_bs,
    )

    libs = [
        all([k0 in lok_bs for k0 in keys]),
        all([k0 not in lok_bs for k0 in keys]),
    ]
    if np.sum(libs) != 1:
        msg = (
            "Either all keys must refer to bsplines or to non-bsplines!\n"
            f"Provided: {keys}"
        )
        raise Exception(msg)
        
    isbs = libs[0]

    # ------------
    # ref_key

    if isbs:
        
        # check ref_key
        wbs = coll._which_bsplines
        lbs = [
            coll.ddata[k0][wbs] for k0 in keys
            if len(coll.dobj[wbs][coll.ddata[k0][wbs]]['ref']) <= maxd
        ]
        lbsu = sorted(set(itt.chain.from_iterable(lbs)))
        lbsok = [
            bs for bs in lbsu
            if all([bs in bb for bb in lbs])
        ]
        
        # ref_key
        ref_key = ds._generic_check._check_var(
            ref_key, 'ref_key',
            types=str,
            allowed=lbsok,
        )
        
        # daxis
        daxis = {
            k0: coll.ddata[k0]['ref'].index(coll.dobj[wbs][ref_key]['ref'][0])
            for k0 in keys
        }
        
        # dunits
        dunits = {
            k0: coll.ddata[coll.dobj[wbs][ref_key]['knots'][0]]['units']
            for k0 in keys
        }
        
        # units_ref
        wbs = coll._which_bsplines
        units_ref = coll.ddata[coll.dobj[wbs][ref_key]['knots'][0]]['units']
                
    # ref_key
    elif ref_key is None:
        hasref, ref, ref_key, val, dkeys = coll.get_ref_vector_common(
            keys=keys,
        )
        if ref_key is None:
            msg = f"No matching ref vector found for {keys}"
            raise Exception(msg)
        ref_key = (ref_key,)

        # daxis
        daxis = {
            k0: [
                coll.ddata[k0]['ref'].index(coll.ddata[rr]['ref'][0])
                for rr in ref_key
            ]
            for k0 in keys
        }
        
        # dunits
        dunits = {k0: coll.ddata[k0]['units'] for k0 in keys}
        
        # units_ref
        units_ref = [coll.ddata[rr]['units'] for rr in ref_key]

    return isbs, keys, ref_key, daxis, dunits, units_ref


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