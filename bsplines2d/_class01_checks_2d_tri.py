# -*- coding: utf-8 -*-


# Built-in
import warnings


# Common
import numpy as np
from matplotlib.tri import Triangulation as mplTri
import scipy.spatial as scpsp
import datastock as ds


# specific
from . import _generic_mesh


# ################################################################
# ################################################################
#                           mesh generic check
# ################################################################


def check(
    coll=None,
    key=None,
    # knots
    knots=None,
    indices=None,
    # from pts
    pts_x0=None,
    pts_x1=None,
    # names of coords
    knots0_name=None,
    knots1_name=None,
    # attributes
    **kwdargs,
):

    # key
    key = ds._generic_check._obj_key(
        d0=coll._dobj.get(coll._which_mesh, {}),
        short='m',
        key=key,
    )

    # -------------
    # knots vs pts

    lc = [
        knots is not None and indices is not None,
        pts_x0 is not None and pts_x1 is not None,
    ]
    if np.sum(lc) != 1:
        msg = (
            "Please provide args (knots, indices) xor (pts_x0, pts_x1)!\n"
            f"\t- knots: {knots}\n"
            f"\t- indices: {indices}\n"
            f"\t- pts_x0: {pts_x0}\n"
            f"\t- pts_x1: {pts_x1}\n"
        )
        raise Exception(msg)

    # --------------------
    # from pts and polygon

    if lc[1]:
        knots, indices = _from_pts_poly(
            pts_x0=pts_x0,
            pts_x1=pts_x1,
        )

    # -------------
    # knots vectors

    knots, indices, ntri = _check_knotscents(
        key=key,
        knots=knots,
        indices=indices,
    )

    cents0 = np.mean(knots[indices, 0], axis=1)
    cents1 = np.mean(knots[indices, 1], axis=1)

    # --------------
    # to dict

    dref, ddata, dobj = _to_dict(
        coll=coll,
        key=key,
        knots=knots,
        indices=indices,
        cents0=cents0,
        cents1=cents1,
        ntri=ntri,
        knots0_name=knots0_name,
        knots1_name=knots1_name,
        **kwdargs,
    )

    return key, dref, ddata, dobj


# ################################################################
# ################################################################
#                      from pts and polygon
# ################################################################


def _from_pts_poly(
    pts_x0=None,
    pts_x1=None,
):

    pts_x0 = ds._generic_check._check_flat1darray(
        pts_x0, 'pts_x0',
        dtype=float,
        can_be_None=False,
    )

    pts_x1 = ds._generic_check._check_flat1darray(
        pts_x1, 'pts_x1',
        dtype=float,
        size=pts_x0.size,
        can_be_None=False,
    )

    # --------
    # Delauney

    knots = np.array([pts_x0, pts_x1]).T

    delau = scpsp.Delaunay(
        knots,
        furthest_site=False,
        incremental=False,
        qhull_options='QJ',
    )

    return knots, delau.simplices


# ################################################################
# ################################################################
#                          check knots and indices
# ################################################################


def _check_knotscents(
    key=None,
    knots=None,
    indices=None,
):

    # ---------------------
    # check mesh conformity

    indices, knots = _mesh2DTri_conformity(knots=knots, indices=indices, key=key)

    # ---------------------------------------------
    # define triangular mesh and trifinder function

    # triangular mesh
    if indices.shape[1] == 3:

        # check clock-wise triangles
        indices = _mesh2DTri_ccw(knots=knots, indices=indices, key=key)
        ntri = 1

    # Quadrangular mesh => convert to triangular
    elif indices.shape[1] == 4:

        ind2 = np.empty((indices.shape[0]*2, 3), dtype=int)
        ind2[::2, :] = indices[:, :3]
        ind2[1::2, :-1] = indices[:, 2:]
        ind2[1::2, -1] = indices[:, 0]
        indices = ind2

        # Re-check mesh conformity
        indices, knots = _mesh2DTri_conformity(knots=knots, indices=indices, key=key)
        indices = _mesh2DTri_ccw(knots=knots, indices=indices, key=key)
        ntri = 2

    # ----------------------------
    # Check on trifinder function

    out = 'not done'
    try:
        trifind = mplTri(knots[:, 0], knots[:, 1], indices).get_trifinder()
        out = trifind(np.r_[0.], np.r_[0])
        assert isinstance(out, np.ndarray)
    except Exception as err:
        msg = str(err) + (
            "\nArg trifind must return an array of indices when fed with arrays "
            "of (R, Z) coordinates!\n"
            f"\ttrifind(np.r_[0], np.r_[0.]) = {out}\n"
            f"\t- ntri = {ntri}\n"
            f"\t- indices = {indices}\n"
            f"\t- knots = {knots}"
        )
        raise Exception(msg)

    return knots, indices, ntri


# ################################################################
# ################################################################
#                           Mesh2DTri
# ################################################################


def _mesh2DTri_conformity(knots=None, indices=None, key=None):

    # ---------------------------------
    # make sure np.ndarrays of dim = 2

    knots = np.atleast_2d(knots).astype(float)
    indices = np.atleast_2d(indices).astype(int)

    # --------------
    # check shapes

    c0 = (
        knots.shape[1] == 2
        and knots.shape[0] >= 3
        and indices.shape[1] in [3, 4]
        and indices.shape[0] >= 1
        and indices.dtype == int
    )
    if not c0:
        msg = (
            "Arg knots must be of shape (nknots>=3, 2) and "
            "arg indices must be of shape (nind>=1, 3 or 4) and dtype = int\n"
            "Provided:\n"
            f"\t- knots.shape: {knots.shape}\n"
            f"\t- indices.shape: {indices.shape}\n"
            f"\t- indices.dtype: {indices.dtype}\n"
        )
        raise Exception(msg)

    nknots = knots.shape[0]
    nind = indices.shape[0]

    # -------------------
    # Test for duplicates

    # knots (floats => distance)
    dist = np.full((nknots, nknots), np.nan)
    ind = np.zeros(dist.shape, dtype=bool)
    for ii in range(nknots):
        dist[ii, ii+1:] = np.sqrt(
            (knots[ii+1:, 0] - knots[ii, 0])**2
            + (knots[ii+1:, 1] - knots[ii, 1])**2
        )
        ind[ii, ii+1:] = True

    ind[ind] = dist[ind] < 1.e-6
    if np.any(ind):
        iind = np.any(ind, axis=1).nonzero()[0]
        lstr = [f'\t\t- {ii}: {ind[ii, :].nonzero()[0]}' for ii in iind]
        msg = (
            f"Non-valid mesh {key}: \n"
            f"  Duplicate knots: {ind.sum()}\n"
            f"\t- knots.shape: {indices.shape}\n"
            f"\t- duplicate indices:\n"
            + "\n".join(lstr)
        )
        raise Exception(msg)

    # cents (indices)
    indu = np.unique(indices, axis=0)
    if indu.shape[0] != nind:
        msg = (
            f"Non-valid mesh {key}: \n"
            f"  Duplicate indices: {nind - indu.shape[0]}\n"
            f"\t- indices.shape: {indices.shape}\n"
            f"\t- unique shape: {indu.shape}"
        )
        raise Exception(msg)

    # -------------------------------
    # Test for unused / unknown knots

    indu = np.unique(indu)
    c0 = np.all(indu >= 0) and indu.size == nknots

    # unused knots
    ino = (~np.in1d(
        range(0, nknots),
        indu,
        assume_unique=False,
        invert=False,
    )).nonzero()[0]

    # unknown knots
    unknown = np.setdiff1d(indu, range(nknots), assume_unique=True)

    if ino.size > 0 or unknown.size > 0:
        msg = f"Knots non-conformity identified for triangular mesh '{key}':\n"
        if ino.size > 0:
            msg += f"\t- Unused knots indices: {ino}\n"
        if unknown.size > 0:
            msg += f"\t- Unknown knots indices: {unknown}\n"
        raise Exception(msg)

    if indu.size < nknots:
        msg = (
            f"Unused knots in {key}:\n"
            f"\t- unused knots indices: {ino}"
        )
        warnings.warn(msg)

    elif indu.size > nknots or indu.max() != nknots - 1:
        unknown = np.setdiff1d(indu, range(nknots), assume_unique=True)
        msg = (
            "Unknown knots refered to in indices!\n"
            f"\t- unknown knots: {unknown}"
        )
        raise Exception(msg)

    return indices, knots


def _mesh2DTri_ccw(knots=None, indices=None, key=None):

    x, y = knots[indices, 0], knots[indices, 1]
    orient = (
        (y[:, 1] - y[:, 0])*(x[:, 2] - x[:, 1])
        - (y[:, 2] - y[:, 1])*(x[:, 1] - x[:, 0])
    )

    clock = orient > 0.
    if np.any(clock):
        msg = (
            "Some triangles are not counter-clockwise\n"
            "  (necessary for matplotlib.tri.Triangulation)\n"
            f"\t=> {clock.sum()} / {indices.shape[0]} triangles reshaped\n"
            f"\t indices: {clock.nonzero()[0]}"
        )
        warnings.warn(msg)
        indices[clock, 1], indices[clock, 2] = indices[clock, 2], indices[clock, 1]
    return indices


def _to_dict(
    coll=None,
    key=None,
    knots=None,
    indices=None,
    cents0=None,
    cents1=None,
    ntri=None,
    # names of coords
    knots0_name=None,
    knots1_name=None,
    # attributes
    **kwdargs,
):

    # ---------
    # check

    knots0_name = ds._generic_check._check_var(
        knots0_name, 'knots0_name',
        types=str,
        default='x',
    )

    knots1_name = ds._generic_check._check_var(
        knots1_name, 'knots1_name',
        types=str,
        default='y',
        excluded=[knots0_name],
    )

    # -----------------
    # keys

    kk = f"{key}-nk"
    kc = f"{key}-nc"
    ki = f"{key}-nind"

    _, _, kk0, kc0 = _generic_mesh.names_knots_cents(
        key=key,
        knots_name='0',
    )
    _, _, kk1, kc1 = _generic_mesh.names_knots_cents(
        key=key,
        knots_name='1',
    )

    kii = f"{key}-ind"

    # attributes
    latt = ['dim', 'quant', 'name', 'units']
    dim, quant, name, units = _generic_mesh._get_kwdargs_2d(kwdargs, latt)

    # -----------------
    # dict

    # dref
    dref = {
        kk: {
            'size': knots.shape[0],
        },
        kc: {
            'size': indices.shape[0],
        },
        ki: {
            'size': 3,
        },
    }

    # ddata
    ddata = {
        kk0: {
            'data': knots[:, 0],
            'units': units[0],
            'quant': quant[0],
            'dim': dim[0],
            'name': name[0],
            'ref': kk,
        },
        kk1: {
            'data': knots[:, 1],
            'units': units[1],
            'quant': quant[1],
            'dim': dim[1],
            'name': name[1],
            'ref': kk,
        },
        kc0: {
            'data': cents0,
            'units': units[0],
            'quant': quant[0],
            'dim': dim[0],
            'name': name[0],
            'ref': kc,
        },
        kc1: {
            'data': cents1,
            'units': units[1],
            'quant': quant[1],
            'dim': dim[1],
            'name': name[1],
            'ref': kc,
        },
        kii: {
            'data': indices,
            # 'units': '',
            'quant': 'indices',
            'dim': 'indices',
            'ref': (kc, ki),
        },
    }

    # dobj
    dobj = {
        coll._which_mesh: {
            key: {
                'nd': '2d',
                'type': 'tri',
                'ntri': ntri,
                'cents': (kc0, kc1),
                'knots': (kk0, kk1),
                'ind': kii,
                # 'ref-k': (kk,),
                # 'ref-c': (kc,),
                'shape-c': (indices.shape[0],),
                'shape-k': (knots.shape[0],),
                'crop': False,
            },
        }
    }

    # additional attributes
    for k0, v0 in kwdargs.items():
        if k0 not in latt:
            dobj[coll._which_mesh][key][k0] = v0

    return dref, ddata, dobj
