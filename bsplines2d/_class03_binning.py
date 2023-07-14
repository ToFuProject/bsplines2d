# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 20:14:40 2023

@author: dvezinet
"""


import warnings


import numpy as np
import datastock as ds


# ############################################################
# ############################################################
#               interpolate spectral
# ############################################################


def binning(
    coll=None,
    data=None,
    data_units=None,
    axis=None,
    # binning
    bins0=None,
    bins1=None,
    bin_data0=None,
    bin_data1=None,
    bin_units0=None,
    # kind of binning
    integrate=None,
    statistic=None,
    # options
    safety_ratio=None,
    dref_vector=None,
    ref_vector_strategy=None,
    verb=None,
    returnas=None,
    # storing
    store=None,
    store_keys=None,
):
    """ Return the spectrally interpolated coefs

    Either E xor Ebins can be provided
    - E: return interpolated coefs
    - Ebins: return binned (integrated) coefs
    """

    # ----------
    # checks

    # keys
    isbs, bin_data0 = _check_bs(
        coll=coll, 
        bin_data0=bin_data0,
        bin_data1=bin_data1,
    )

    # ----------
    # trivial

    if isbs:
        
        # add ref and data
        kr, kd, ddatan = _interpolate(
            coll=coll,
            data=data,
            data_units=data_units,
            # binning
            bins0=bins0,
            bin_data0=bin_data0,
            # options
            dref_vector=dref_vector,
            verb=verb,
            store=store,
        )
        
        lk = list(ddatan.keys())
        data = [ddatan[k0]['data'] for k0 in lk]
        bin_data0 = [ddatan[k0]['bin_data'] for k0 in lk]
        
    dout = ds._class1_binning.binning(
        coll=coll,
        data=data,
        data_units=data_units,
        axis=axis,
        # binning
        bins0=bins0,
        bins1=bins1,
        bin_data0=bin_data0,
        bin_data1=bin_data1,
        bin_units0=bin_units0,
        # kind of binning
        integrate=integrate,
        statistic=statistic,
        # options
        safety_ratio=safety_ratio,
        dref_vector=dref_vector,
        ref_vector_strategy=ref_vector_strategy,
        verb=verb,
        returnas=returnas,
        # storing
        store=store,
        store_keys=store_keys,
    )

    # --------------------------------
    # remove intermediate ref and data

    if isbs is True:
        for dd in data + bin_data0 + [kd]:
            if dd in coll.ddata.keys():
                coll.remove_data(dd)
        if kr in coll.dref.keys():
            coll.remove_ref(kr)
            
        for k0 in data:
            k1 = [k1 for k1, v1 in ddatan.items() if v1['data'] == k0][0]
            dout[k1] = dict(dout[k0])
            del dout[k0]

    # ----------
    # return

    return dout


# ######################################################
# ######################################################
#                   check
# ######################################################


def _check_bs(
    coll=None,
    bin_data0=None,
    bin_data1=None,
):
    
    wbs = coll._which_bsplines
    lok_bs = [
        k0 for k0, v0 in coll.dobj.get(wbs, {}).items()
        if len(v0['ref']) == 1
    ]
    lok_dbs = [
        k0 for k0, v0 in coll.ddata.items()
        if v0.get(wbs) is not None
        and len(v0[wbs]) == 1
        and v0[wbs][0] in coll.dobj.get(wbs, {}).keys()
        and len(coll.dobj[wbs][v0[wbs][0]]['ref']) == 1
    ]
        
    c0 = (
        isinstance(bin_data0, str)
        and bin_data1 is None
        and bin_data0 in lok_dbs + lok_bs
    )
    
    if bin_data0 in lok_bs:
        bin_data0 = coll.dobj[wbs][bin_data0]['apex'][0]
    
    return c0, bin_data0
        

# ######################################################
# ######################################################
#                   interpolate
# ######################################################


def _interpolate(
    coll=None,
    data=None,
    data_units=None,
    # binning
    bins0=None,
    bin_data0=None,
    # options
    dref_vector=None,
    verb=None,
    store=None,
):

    # ---------
    # sampling

    # mesh knots
    wm = coll._which_mesh
    wbs = coll._which_bsplines
    key_bs = coll.ddata[bin_data0][wbs][0]
    keym = coll.dobj[wbs][key_bs][wm]
    kknots = coll.dobj[wm][keym]['knots'][0]

    # resolution
    vect = coll.ddata[kknots]['data']
    res0 = np.abs(np.min(np.diff(vect)))

    # ---------
    # sampling
    
    ddata = ds._class1_binning._check_data(
        coll=coll,
        data=data,
        data_units=data_units,
        store=True,
    )
    lkdata = list(ddata.keys())
    
    # --------------------
    # bins

    dbins0 = ds._class1_binning._check_bins(
        coll=coll,
        lkdata=lkdata,
        bins=bins0,
        dref_vector=dref_vector,
        store=store,
    )

    # ----------------------
    # npts for interpolation
    
    dv = np.abs(np.diff(vect))
    dvmean = np.mean(dv) + np.std(dv)
    db = np.mean(np.diff(dbins0[lkdata[0]]['edges']))
    npts = (coll.dobj[wbs][key_bs]['deg'] + 3) * max(1, dvmean / db)

    # sample mesh, update dv
    xx = coll.get_sample_mesh(
        keym,
        res=res0 / npts,
        mode='abs',
        Dx0=[bins0[0], bins0[-1]],
    )['x0']['data']

    # -------------------
    #  add ref
    
    kr = "ntemp"
    kd = "xxtemp"
    
    coll.add_ref(kr, size=xx.size)
    coll.add_data(kd, data=xx, ref=kr, units=coll.ddata[kknots]['units'])

    ddata_new = {}
    for ii, (k0, v0) in enumerate(ddata.items()):

        # interpolate bin_data
        kbdn = f"kbdn{ii}_temp"
        coll.interpolate(
            keys=bin_data0,
            ref_key=key_bs,
            x0=kd,
            val_out=0.,
            returnas=False,
            store=True,
            inplace=True,
            store_keys=kbdn,
        )
        
        # interpolate_data
        kdn = f"kbd{ii}_temp"
        coll.interpolate(
            keys=k0,
            ref_key=key_bs,
            x0=kd,
            val_out=0.,
            returnas=False,
            store=True,
            inplace=True,
            store_keys=kdn,
        )
        ddata_new[k0] = {'bin_data': kbdn, 'data': kdn}

    return kr, kd, ddata_new