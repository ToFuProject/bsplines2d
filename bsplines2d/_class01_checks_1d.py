# -*- coding: utf-8 -*-


# Common
import numpy as np
import datastock as ds


# specific
from . import _generic_mesh


# #############################################################################
# #############################################################################
#                           mesh generic check
# #############################################################################


def check(
    coll=None,
    key=None,
    # knots
    knots=None,
    uniform=None,
    # defined from pre-existing bsplines
    subkey=None,
    # additional attributes
    **kwdargs,
):

    # --------
    # keys

    # key
    key = ds._generic_check._obj_key(
        d0=coll._dobj.get(coll._which_mesh, {}),
        short='m',
        key=key,
    )

    # ------------
    # knots vector

    knots, res = _check_knots(
        key=key,
        knots=knots,
        uniform=uniform,
    )

    # ----------------
    # angles handdling

    isangle = str(kwdargs.get('units')) == 'rad'
    if isangle:
        knots, cents = _knots_angle(knots)

    else:
        cents = 0.5*(knots[1:] + knots[:-1])

    # ------------------------
    # depend on other bsplines

    submesh, subbs, kwdargs = _defined_from(
        coll=coll,
        subkey=subkey,
        # parameters
        kwdargs=kwdargs,
    )

    # --------------
    # to dict

    dref, ddata, dobj = _to_dict(
        coll=coll,
        key=key,
        knots=knots,
        cents=cents,
        res=res,
        # sub quantity
        subkey=subkey,
        subbs=subbs,
        submesh=submesh,
        # attributes
        **kwdargs,
    )

    return key, dref, ddata, dobj


# ##################################################################
# ##################################################################
#                       knots
# ##################################################################


def _check_knots(
    key=None,
    knots=None,
    uniform=None,
):

    # ------------
    # check input

    uniform = ds._generic_check._check_var(
        uniform, 'uniform',
        types=bool,
        default=False,
    )

    # ---------
    # check x

    knots = ds._generic_check._check_flat1darray(
        knots, 'knots',
        dtype=float,
        unique=True,
        can_be_None=False,
    )

    # resolution
    res = np.diff(knots)

    # -----------------
    # check uniformity

    if np.allclose(res, np.mean(res), atol=1e-14, rtol=1e-6):
        res = res[0]

    elif uniform is True:
        msg = (
            f"Non-uniform resolution for user-provided mesh '{key}'\n"
            f"\t- unique res: {np.unique(res)}\n"
            f"\t- diff res: {np.diff(np.unique(res))}\n"
            f"\t- res: {res}\n"
            )
        raise NotImplementedError(msg)

    return knots, res


def _knots_angle(
    knots=None,
    res=None,
):

    # knots in ]-pi; pi]
    knots_temp = np.unique(np.arctan2(np.sin(knots), np.cos(knots)))
    if not np.allclose(knots, knots_temp):
        msg = (
            "Angle knots must be in ]-pi; pi]!\n"
            f"Provided: {knots}"
        )
        raise Exception(msg)

    # cents - handle discontinuity at -pi
    cents = 0.5*(knots[1:] + knots[:-1])
    mid = 0.5*(knots[-1] + (2.*np.pi + knots[0]))
    mid = np.arctan2(np.sin(mid), np.cos(mid))
    if mid < cents[0]:
        cents = np.r_[mid, cents]
    else:
        cents = np.r_[cents, mid]

    return knots, cents


# #############################################################################
# #############################################################################
#                        defined_from
# #############################################################################


def _defined_from(
    coll=None,
    subkey=None,
    nd=None,
    # parameters
    kwdargs=None,
):

    # ------------
    # trivial

    if subkey is None:
        return None, None, kwdargs

    # ------------
    # check key_on

    wm = coll._which_mesh
    wbs = coll._which_bsplines
    if coll.dobj.get(wbs) is not None:
        dbs, dref = coll.get_dict_bsplines()
        lok =[
            k0 for k0, v0 in dbs.items()
            if len(v0) == 1
        ]

        if nd is not None:
            lok2 = []
            for k0 in lok:
                kbs = list(dbs[k0].keys())[0]
                if coll.dobj[wm][coll.dobj[wbs][kbs]['mesh']]['nd'] == nd:
                    lok2.append(k0)

            lok = lok2

    else:
        lok = []

    subkey = ds._generic_check._check_var(
        subkey, 'subkey',
        types=str,
        allowed=lok,
    )

    # ----------------
    # complete kwdargs

    lq = ['dim', 'quant', 'name', 'units']
    for k0 in lq:
        if kwdargs.get(k0) is None:
            kwdargs[k0] = str(coll.ddata[subkey][k0])

    # --------------
    # key_submesh

    subbs = coll.ddata[subkey][wbs][0]
    submesh = coll.dobj[wbs][subbs]['mesh']

    return submesh, subbs, kwdargs


# #############################################################################
# #############################################################################
#                           to_dict
# #############################################################################


def _to_dict(
    coll=None,
    key=None,
    knots=None,
    cents=None,
    res=None,
    # submesh
    subkey=None,
    subbs=None,
    submesh=None,
    # attributes
    **kwdargs,
):

    # ---------
    # prepare

    # keys
    # keys of knots and cents
    kkr, kcr, kk, kc = _generic_mesh.names_knots_cents(key=key)

    # variable
    variable = not np.isscalar(res)

    # attributes
    latt = ['dim', 'quant', 'name', 'units']
    dim, quant, name, units = [kwdargs.get(ss) for ss in latt]

    # subkey
    if subkey is not None:
        subkey = (subkey,)
    else:
        subkey = None

    # -------------
    # prepare dict

    # dref
    dref = {
        kkr: {
            'size': knots.size,
        },
        kcr: {
            'size': cents.size,
        },
    }

    # ddata
    ddata = {
        kk: {
            'data': knots,
            'units': units,
            # 'source': None,
            'dim': dim,
            'quant': quant,
            'name': name,
            'ref': kkr,
        },
        kc: {
            'data': cents,
            'units': units,
            # 'source': None,
            'dim': dim,
            'quant': quant,
            'name': name,
            'ref': kcr,
        },
    }

    # dobj
    dobj = {
        coll._which_mesh: {
            key: {
                'nd': '1d',
                'type': None,
                'knots': (kk,),
                'cents': (kc,),
                'shape-c': (cents.size,),
                'shape-k': (knots.size,),
                'variable': variable,
                'subkey': subkey,
                'subbs': subbs,
                'submesh': submesh,
                'crop': False,
            },
        },
    }

    # additional attributes
    for k0, v0 in kwdargs.items():
        if k0 not in latt:
            dobj[coll._which_mesh][key][k0] = v0

    return dref, ddata, dobj
