# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 20:14:40 2023

@author: dvezinet
"""


import numpy as np
import datastock as ds


# specific
from . import _class02_interpolate as _interpolate


# ############################################################
# ############################################################
#               interpolate spectral
# ############################################################


def binning(
    coll=None,
    keys=None,
    ref_key=None,
    bins=None,
    # optional storing
    verb=None,
    returnas=None,
    store=None,
    key_store=None,
):
    """ Return the spectrally interpolated coefs

    Either E xor Ebins can be provided
    - E: return interpolated coefs
    - Ebins: return binned (integrated) coefs
    """

    # ----------
    # checks

    refb, bins, verb, store, returnas = _check(
        coll=coll,
        ref_key=ref_key,
        bins=bins,
        verb=verb,
        store=store,
        returnas=returnas,
    )

    # keys
    (
        isbs, keys, ref_key,
        daxis, dunits, units_ref,
        _, _,
    ) = _interpolate._check_keys(
        coll=coll,
        keys=keys,
        ref_key=ref_key,
        only1d=False,
        details=False,
    )

    # ----------
    # trivial

    if not isbs:
        dout = ds._class1_binning.binning(
            coll=coll,
            keys=keys,
            ref_key=ref_key[0],
            bins=bins,
        )

    else:
        dout = _binning(
            coll=coll,
            keys=keys,
            ref_key=ref_key,
            bins=bins,
            dunits=dunits,
            units_ref=units_ref,
            daxis=daxis,
        )

    # ------------
    # adjust ref

    if refb is not None:
        for k0, v0 in dout.items():
            ref = list(coll.ddata[k0]['ref'])
            ref[daxis[k0][0]] = refb
            dout[k0]['ref'] = tuple(ref)

    # ----------------
    # optional verb

    if verb is True:
        # TBD
        pass

    # ----------------
    # optional storing

    if store is True:
        _store(
            coll=coll,
            dout=dout,
            key_store=key_store,
        )

    # ----------
    # return

    if returnas is True:
        return dout


# ######################################################
# ######################################################
#                   check
# ######################################################


def _check(
    coll=None,
    keys=None,
    bins=None,
    ref_key=None,
    verb=None,
    returnas=None,
    store=None,
):

    # catch bin as str
    wb = coll._which_bins
    refb = None
    if isinstance(bins, str):
        lok = list(coll.dobj.get(wb, {}).keys())
        if bins in lok:
            refb = coll.dobj[wb][bins]['ref'][0]
            bins = coll.dobj[wb][bins]['edges']
        else:
            msg = (
                "Arg bins refers to an unknown bins vector!\n"
                f"\t- Provided: '{bins}'\n"
                f"\t- Available: {lok}"
            )
            raise Exception(msg)

    # verb
    verb = ds._generic_check._check_var(
        verb, 'verb',
        types=bool,
        default=True,
    )

    # store
    lok = [False]
    if refb is not None:
        lok.append(True)
    store = ds._generic_check._check_var(
        store, 'store',
        types=bool,
        allowed=lok,
        default=False,
    )

    # returnas
    lok = [False, dict]
    returnas = ds._generic_check._check_var(
        returnas, 'returnas',
        types=bool,
        default=store is False,
    )

    return refb, bins, verb, store, returnas


# ######################################################
# ######################################################
#                   binning
# ######################################################


def _binning(
    coll=None,
    keys=None,
    ref_key=None,
    bins=None,
    dunits=None,
    units_ref=None,
    daxis=None,
):

    # ---------
    # sampling

    # mesh knots
    wm = coll._which_mesh
    wbs = coll._which_bsplines
    keym = coll.dobj[wbs][ref_key][wm]
    kknots = coll.dobj[wm][keym]['knots'][0]

    # resolution
    vect = coll.ddata[kknots]['data']
    res0 = np.abs(np.min(np.diff(vect)))

    # ----------
    # bins

    # bins
    bins, units_bins, _, npts = ds._class1_binning._check_bins(
        coll=coll,
        keys=keys,
        ref_key=ref_key,
        bins=bins,
        vect=vect,
        strict=False,
        deg=coll.dobj[wbs][ref_key]['deg'],
    )

    # sample mesh, update dv
    xx = coll.get_sample_mesh(keym, res=res0 / npts, mode='abs')['x0']['data']
    dv = np.abs(np.diff(xx))

    dv = np.append(dv, dv[-1])

    # units
    dout = ds._class1_binning._units(
        dunits=dunits,
        units_ref=units_ref,
        units_bins=units_bins,
    )

    # --------------
    # actual binning

    for k0, v0 in dout.items():

        # interpolate
        val = coll.interpolate(
            keys=k0,
            ref_key=ref_key,
            x0=xx,
            val_out=0.,
        )[k0]['data']

        # bin
        dout[k0]['data'] = ds._class1_binning._bin(
            bins=bins,
            dv=dv,
            vect=xx,
            data=val,
            axis=daxis[k0][0],
        )

    return dout


# ######################################################
# ######################################################
#                   storing
# ######################################################


def _store(
    coll=None,
    dout=None,
    key_store=None,
):

    # ---------
    # key_store

    if key_store is None:
        key_store = [f'{key}_bin' for key in dout.keys()]
    if isinstance(key_store, str):
        key_store = [key_store]

    lout = list(coll.ddata.keys())
    key_store = ds._generic_check._check_var_iter(
        key_store, 'key_store',
        types=(list, tuple),
        types_iter=str,
        excluded=lout,
    )

    if len(key_store) != len(dout):
        msg = (
            "Arg key_store must be of same length as keys!\n"
            f"\t- keys: {list(dout.keys())}\n"
            f"\t- key_store: {key_store}\n"
        )
        raise Exception(msg)

    # --------
    # store

    for ii, (k0, v0) in enumerate(dout.items()):
        coll.add_data(
            key=key_store[ii],
            **v0,
        )
