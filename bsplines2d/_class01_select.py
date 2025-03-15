# -*- coding: utf-8 -*-


# Built-in
import warnings

# Common
import numpy as np
from scipy.spatial import ConvexHull
from matplotlib.path import Path
from contourpy import contour_generator
import datastock as ds


# tofu
from . import _generic_mesh
from . import _utils_bsplines
# from . import _class02_checks as _checks
from . import _class02_bsplines_rect
from . import _class02_bsplines_tri
from . import _class02_bsplines_polar
from . import _class02_bsplines_1d


# #############################################################################
# #############################################################################
#                           index
# #############################################################################


def _select_ind(
    coll=None,
    key=None,
    ind=None,
    elements=None,
    returnas=None,
    crop=None,
    # spatial vs spectral
):
    """ return indices
    ind can be:
        - None
        - tuple: (R, Z), possibly 2d
        - 'tuple-flat': (R, Z) flattened
        - np.ndarray: array of unique indices
        - 'array-flat': flattened ordered array of unique
    """

    # ------------
    # check inputs

    km, keybs, cat = _generic_mesh._get_key_mesh_vs_bplines(
        coll=coll,
        key=key,
    )

    meshtype = coll.dobj[coll._which_mesh][km]['type']

    shape2d = None
    if cat == coll._which_mesh:
        key = km
        shape2d = len(coll.dobj[cat][key]['shape_c']) == 2
    else:
        key = keybs
        shape2d = len(coll.dobj[cat][key]['shape']) == 2

    # ind, elements, ...
    # elements = cents or knots
    ind, elements, returnas, crop = _generic_mesh._select_ind_check(
        ind=ind,
        elements=elements,
        returnas=returnas,
        crop=crop,
        shape2d=shape2d,
    )

    elem = f'{elements}' if cat == coll._which_mesh else 'ref'

    if shape2d:
        if cat == coll._which_mesh:
            ke = f'shape_{elem[0]}'
            nR, nZ = coll.dobj[cat][key][ke]
        else:
            nR, nZ = coll.dobj[cat][key]['shape']
    else:
        if cat == coll._which_mesh:
            ke = f'shape_{elem[0]}'
            nelem = coll.dobj[cat][key][ke][0]
        else:
            nelem = coll.dobj[cat][key]['shape'][0]

    # ------------
    # ind to tuple

    if shape2d:
        ind_bool = np.zeros((nR, nZ), dtype=bool)
        if ind is None:
            # make sure R is varying in dimension 0
            ind_tup = (
                np.repeat(np.arange(0, nR)[:, None], nZ, axis=1),
                np.tile(np.arange(0, nZ), (nR, 1)),
            )
            ind_bool[...] = True

        elif isinstance(ind, tuple):
            c0 = (
                np.all((ind[0] >= 0) & (ind[0] < nR))
                and np.all((ind[1] >= 0) & (ind[1] < nZ))
            )
            if not c0:
                msg = (
                    f"Non-valid values in ind (< 0 or >= size ({nR}, {nZ}))"
                )
                raise Exception(msg)
            ind_tup = ind
            ind_bool[ind_tup[0], ind_tup[1]] = True

        else:
            if np.issubdtype(ind.dtype, np.integer):
                c0 = np.all((ind >= 0) & (ind < nR*nZ))
                if not c0:
                    msg = (
                        f"Non-valid values in ind (< 0 or >= size ({nR*nZ}))"
                    )
                    raise Exception(msg)
                ind_tup = (ind % nR, ind // nR)
                ind_bool[ind_tup[0], ind_tup[1]] = True

            elif np.issubdtype(ind.dtype, np.bool_):
                if ind.shape != (nR, nZ):
                    msg = (
                        f"Arg ind, if bool, must have shape {(nR, nZ)}\n"
                        f"Provided: {ind.shape}"
                    )
                    raise Exception(msg)
                # make sure R varies first
                ind_tup = ind.T.nonzero()[::-1]
                ind_bool = ind

            else:
                msg = f"Unknown ind dtype!\n\t- ind.dtype: {ind.dtype}"
                raise Exception(msg)

        if ind_tup[0].shape != ind_tup[1].shape:
            msg = (
                "ind_tup components do not have the same shape!\n"
                f"\t- ind_tup[0].shape = {ind_tup[0].shape}\n"
                f"\t- ind_tup[1].shape = {ind_tup[1].shape}"
            )
            raise Exception(msg)

    # triangular + polar1d case
    else:
        ind_bool = np.zeros((nelem,), dtype=bool)
        if ind is None:
            ind_bool[...] = True
        elif np.issubdtype(ind.dtype, np.integer):
            c0 = np.all((ind >= 0) & (ind < nelem))
            if not c0:
                msg = (
                    f"Arg ind has non-valid values (< 0 or >= size ({nelem}))"
                )
                raise Exception(msg)
            ind_bool[ind] = True
        elif np.issubdtype(ind.dtype, np.bool_):
            if ind.shape != (nelem,):
                msg = (
                    f"Arg ind, when array of bool, must have shape {(nelem,)}"
                    f"\nProvided: {ind.shape}"
                )
                raise Exception(msg)
            ind_bool = ind
        else:
            msg = (
                "Non-valid ind format!"
            )
            raise Exception(msg)

    # ------------
    # optional crop

    crop = (
        crop is True
        and coll.dobj[cat][key].get('crop') not in [None, False]
        and bool(np.any(~coll.ddata[coll.dobj[cat][key]['crop']]['data']))
    )
    if crop is True:
        cropi = coll.ddata[coll.dobj[cat][key]['crop']]['data']
        if meshtype == 'rect':
            if cat == coll._which_mesh and elements == 'knots':
                cropiknots = np.zeros(ind_bool.shape, dtype=bool)
                cropiknots[:-1, :-1] = cropi
                cropiknots[1:, :-1] = cropiknots[1:, :-1] | cropi
                cropiknots[1:, 1:] = cropiknots[1:, 1:] | cropi
                cropiknots[:-1, 1:] = cropiknots[:-1, 1:] | cropi

                ind_bool = ind_bool & cropiknots

                # ind_tup is not 2d anymore
                ind_tup = ind_bool.T.nonzero()[::-1]  # R varies first
                # warnings.warn("ind is not 2d anymore!")

            elif ind_tup[0].shape == cropi.shape:
                ind_bool = ind_bool & cropi
                # ind_tup is not 2d anymore
                ind_tup = ind_bool.T.nonzero()[::-1]  # R varies first
                # warnings.warn("ind is not 2d anymore!")

            else:
                ind_bool = ind_bool & cropi
                ind_tup = ind_bool.T.nonzero()[::-1]
        else:
            ind_bool &= cropi

    # ------------
    # tuple to return

    if returnas is bool:
        out = ind_bool
    elif returnas is int:
        out = ind_bool.nonzero()[0]
    elif returnas is tuple:
        out = ind_tup
    elif returnas == 'tuple-flat':
        # make sure R is varying first
        out = (ind_tup[0].T.ravel(), ind_tup[1].T.ravel())
    elif returnas is np.ndarray:
        out = ind_tup[0] + ind_tup[1]*nR
    elif returnas == 'array-flat':
        # make sure R is varying first
        out = (ind_tup[0] + ind_tup[1]*nR).T.ravel()
    else:
        out = ind_bool

    return out


# #############################################################################
# #############################################################################
#                           mesh elements
# #############################################################################


def _select_mesh(
    coll=None,
    key=None,
    ind=None,
    elements=None,
    returnas=None,
    return_ind_as=None,
    return_neighbours=None,
):
    """ ind is a tuple for rect """

    # ------------
    # check inputs

    wm = coll._which_mesh
    mtype = coll.dobj[wm][key]['type']

    (
        elements, returnas,
        return_ind_as, return_neighbours,
    ) = _generic_mesh._select_check(
        elements=elements,
        returnas=returnas,
        return_ind_as=return_ind_as,
        return_neighbours=return_neighbours,
    )

    # ------------
    # prepare

    if mtype in ['rect', 'tri']:
        kR, kZ = coll.dobj[wm][key][elements]
        R = coll.ddata[kR]['data']
        Z = coll.ddata[kZ]['data']
    else:
        kr = coll.dobj[wm][key][elements][0]
        rad = coll.ddata[kr]['data']

    # ------------
    # non-trivial case

    if returnas == 'ind':
        out = ind
    else:
        if mtype == 'rect':
            out = R[ind[0]], Z[ind[1]]
        elif mtype == 'tri':
            out = R[ind], Z[ind]
        else:
            out = rad[ind]

    # ------------
    # neighbours

    if return_neighbours is True:
        if mtype == 'rect':
            neigh = _select_mesh_neighbours_rect(
                coll=coll,
                key=key,
                ind=ind,
                elements=elements,
                returnas=returnas,
            )
        elif mtype == 'tri':
            neigh = _select_mesh_neighbours_tri(
                coll=coll,
                key=key,
                ind=ind,
                elements=elements,
                returnas=returnas,
                return_ind_as=return_ind_as,
            )
        else:
            neigh = _select_mesh_neighbours_1d(
                coll=coll,
                key=key,
                ind=ind,
                elements=elements,
                returnas=returnas,
                return_ind_as=return_ind_as,
            )

        return out, neigh
    else:
        return out


# TBF
def _select_mesh_neighbours_1d(
    coll=None,
    key=None,
    ind=None,
    elements=None,
    returnas=None,
    return_ind_as=None,
):
    """ ind is a bool

    if returnas = 'ind', ind is returned as a bool array
    (because the nb. of neighbours is not constant on a triangular mesh)

    """
    # ------------
    # neighbours

    elneig = 'cents' if elements == 'knots' else 'knots'
    kneig = coll.dobj[coll._which_mesh][key][f'{elneig}'][0]
    neig = coll.ddata[kneig]['data']
    nneig = neig.size

    # get tuple indices of neighbours
    ind = ind.nonzero()[0]
    shape = tuple(np.r_[ind.shape, 2].astype(int))
    ineig = np.zeros(shape, dtype=int)

    rsh = tuple(
        [2 if ii == len(shape)-1 else 1 for ii in range(len(shape))]
    )

    if elements == 'cents':
        ineig[...] = ind[..., None] + np.r_[0, 1].reshape(rsh)

    elif elements == 'knots':
        ineig[...] = ind[..., None] + np.r_[-1, 0].reshape(rsh)
        ineig[(ineig < 0) | (ineig >= nneig)] = -1

    # return neighbours in desired format
    if returnas == 'ind':
        neig_out = ineig
    else:
        neig_out = neig[ineig]
        neig_out[ineig == -1] = np.nan

    return neig_out



def _select_mesh_neighbours_rect(
    coll=None,
    key=None,
    ind=None,
    elements=None,
    returnas=None,
):
    """ ind is a tuple for rect """

    # ------------
    # neighbours

    elneig = 'cents' if elements == 'knots' else 'knots'
    kRneig, kZneig = coll.dobj[coll._which_mesh][key][f'{elneig}']
    Rneig = coll.ddata[kRneig]['data']
    Zneig = coll.ddata[kZneig]['data']
    nRneig = Rneig.size
    nZneig = Zneig.size

    # get tuple indices of neighbours
    shape = tuple(np.r_[ind[0].shape, 4].astype(int))
    neig = (
        np.zeros(shape, dtype=int),
        np.zeros(shape, dtype=int),
    )
    rsh = tuple(
        [4 if ii == len(shape)-1 else 1 for ii in range(len(shape))]
    )

    if elements == 'cents':
        neig[0][...] = ind[0][..., None] + np.r_[0, 1, 1, 0].reshape(rsh)
        neig[1][...] = ind[1][..., None] + np.r_[0, 0, 1, 1].reshape(rsh)
    elif elements == 'knots':
        neig[0][...] = ind[0][..., None] + np.r_[-1, 0, 0, -1].reshape(rsh)
        neig[1][...] = ind[1][..., None] + np.r_[-1, -1, 0, 0].reshape(rsh)
        neig[0][(neig[0] < 0) | (neig[0] >= nRneig)] = -1
        neig[1][(neig[1] < 0) | (neig[1] >= nZneig)] = -1

    # return neighbours in desired format
    if returnas == 'ind':
        neig_out = neig
    else:
        neig_out = np.array([Rneig[neig[0]], Zneig[neig[1]]])
        neig_out[:, (neig[0] == -1) | (neig[1] == -1)] = np.nan

    return neig_out


def _select_mesh_neighbours_tri(
    coll=None,
    key=None,
    ind=None,
    elements=None,
    returnas=None,
    return_ind_as=None,
):
    """ ind is a bool

    if returnas = 'ind', ind is returned as a bool array
    (because the nb. of neighbours is not constant on a triangular mesh)

    """
    # ------------
    # neighbours

    nind = ind.sum()
    kind = coll.dobj[coll._which_mesh][key]['ind']

    if returnas == 'data':
        elneig = 'cents' if elements == 'knots' else 'knots'
        kneig = coll.dobj[coll._which_mesh][key][elneig]
        Rneig = coll.ddata[kneig[0]]['data']
        Zneig = coll.ddata[kneig[1]]['data']

    if elements == 'cents':
        neig = coll.ddata[kind]['data'][ind, :]
        if returnas == 'ind':
            if return_ind_as is bool:
                kknots = coll.dobj[coll._which_mesh][key]['knots']
                nneig = coll.dref[f'{kknots}-ind']['size']
                neig_temp = np.zeros((nind, nneig), dtype=bool)
                for ii in range(nind):
                    neig_temp[ii, neig[ii, :]] = True
                neig = neig_temp
        else:
            neig = np.array([Rneig[neig], Zneig[neig]])
    else:
        ind_int = ind.nonzero()[0]
        neig = np.array([
            np.any(coll.ddata[kind]['data'] == ii, axis=1)
            for ii in ind_int
        ])
        c0 = returnas == 'ind' and return_ind_as is int
        if c0 or returnas == 'data':
            nmax = np.sum(neig, axis=1)
            if returnas == 'ind':
                neig_temp = -np.ones((nind, nmax.max()), dtype=int)
                for ii in range(nind):
                    neig_temp[ii, :nmax[ii]] = neig[ii, :].nonzero()[0]
            else:
                neig_temp = np.full((2, nind, nmax.max()), np.nan)
                for ii in range(nind):
                    neig_temp[0, ii, :nmax[ii]] = Rneig[neig[ii, :]]
                    neig_temp[1, ii, :nmax[ii]] = Zneig[neig[ii, :]]
            neig = neig_temp

    return neig


# TBF
# def _select_mesh_neighbours_polar(
    # coll=None,
    # key=None,
    # ind=None,
    # elements=None,
    # returnas=None,
    # return_ind_as=None,
# ):
    # """ ind is a bool

    # if returnas = 'ind', ind is returned as a bool array
    # (because the nb. of neighbours is not constant on a triangular mesh)

    # """

    # elneig = 'cents' if elements == 'knots' else 'knots'
    # kneig = coll.dobj[coll._which_mesh][key][f'{elneig}']
    # rneig = coll.ddata[kneig[0]]['data']
    # nrneig = rneig.size


    # # ----------------
    # # radius + angle

    # if len(kneig) == 2:
        # aneig = coll.ddata[kneig[1]]['data']
        # naneig = aneig.size

        # # prepare indices
        # shape = tuple(np.r_[ind[0].shape, 2])
        # neig = (
            # np.zeros((nrneig, 2), dtype=bool),
            # np.zeros((naneig, 2), dtype=bool),
        # )

        # # get indices of neighbours
        # if elements == 'cents':
            # neig[0][...] = ind[0][..., None] + np.r_[0, 1, 1, 0].reshape(rsh)
            # neig[1][...] = ind[1][..., None] + np.r_[0, 0, 1, 1].reshape(rsh)
        # elif elements == 'knots':
            # neig[0][...] = ind[0][..., None] + np.r_[-1, 0, 0, -1].reshape(rsh)
            # neig[1][...] = ind[1][..., None] + np.r_[-1, -1, 0, 0].reshape(rsh)
            # neig[0][(neig[0] < 0) | (neig[0] >= nRneig)] = -1
            # neig[1][(neig[1] < 0) | (neig[1] >= nZneig)] = -1


    # # ----------------
    # # radius only

    # else:
        # # prepare indices
        # neig = np.zeros((nrneig, 2), dtype=bool)

        # # get indices of neighbours
        # if elements == 'cents':
            # neig[0][...] = ind[0][..., None] + np.r_[0, 1, 1, 0].reshape(rsh)
            # neig[1][...] = ind[1][..., None] + np.r_[0, 0, 1, 1].reshape(rsh)

    # # return neighbours in desired format
    # if returnas == 'ind':
        # neig_out = neig
    # else:
        # if len(kneig) == 2:
            # neig_out = np.array([rneig[neig[0]], zneig[neig[1]]])
            # neig_out[:, (neig[0] == -1) | (neig[1] == -1)] = np.nan
        # else:
            # neig_out = rneig[neig]

    # return neig_out


# #############################################################################
# #############################################################################
#                           bsplines
# #############################################################################


def _select_bsplines(
    coll=None,
    key=None,
    ind=None,
    returnas=None,
    return_cents=None,
    return_knots=None,
    crop=None,
):
    """ ind is a tuple """

    # ------------
    # check inputs

    wm = coll._which_mesh
    wbs = coll._which_bsplines

    # check key
    keym, keybs, _ = _generic_mesh._get_key_mesh_vs_bplines(
        coll=coll,
        key=key,
        which=wbs,
    )

    nd = coll.dobj[wm][keym]['nd']
    mtype = coll.dobj[wm][keym]['type']

    # ----
    # ind

    if mtype == 'rect':
        returnasind = tuple
    else:
        returnasind = bool

    ind = _select_ind(
        coll=coll,
        key=key,
        ind=ind,
        elements=None,
        returnas=returnasind,
        crop=crop,
    )

    # ------------
    # knots, cents

    if mtype == 'rect':
        kRk, kZk = coll.dobj[wm][keym]['knots']
        kRc, kZc = coll.dobj[wm][keym]['cents']

        out = _mesh2DRect_bsplines_knotscents(
            returnas=returnas,
            return_knots=return_knots,
            return_cents=return_cents,
            ind=ind,
            deg=coll.dobj[wbs][key]['deg'],
            Rknots=coll.ddata[kRk]['data'],
            Zknots=coll.ddata[kZk]['data'],
            Rcents=coll.ddata[kRc]['data'],
            Zcents=coll.ddata[kZc]['data'],
        )

    elif mtype == 'tri':
        clas = coll.dobj[wbs][key]['class']
        out = clas._get_knotscents_per_bs(
            returnas=returnas,
            return_knots=return_knots,
            return_cents=return_cents,
            ind=ind,
        )

    else:
        kk = coll.dobj[wm][keym]['knots'][0]
        kc = coll.dobj[wm][keym]['cents'][0]

        out = _mesh1d_bsplines_knotscents(
            returnas=returnas,
            return_knots=return_knots,
            return_cents=return_cents,
            ind=ind,
            deg=coll.dobj[wbs][key]['deg'],
            knots=coll.ddata[kk]['data'],
            cents=coll.ddata[kc]['data'],
        )

    # ------------
    # return

    if return_cents is True and return_knots is True:
        return ind, out[0], out[1]
    elif return_cents is True or return_knots is True:
        return ind, out
    else:
        return ind


def _mesh1d_bsplines_knotscents(
    returnas=None,
    return_knots=None,
    return_cents=None,
    ind=None,
    deg=None,
    knots=None,
    cents=None,
):

    # -------------
    # check inputs

    return_knots = ds._generic_check._check_var(
        return_knots, 'return_knots',
        types=bool,
        default=True,
    )
    return_cents = ds._generic_check._check_var(
        return_cents, 'return_cents',
        types=bool,
        default=True,
    )
    if return_knots is False and return_cents is False:
        return

    # -------------
    # compute

    if return_knots is True:

        knots_per_bs = _utils_bsplines._get_knots_per_bs(
            knots, deg=deg, returnas=returnas,
        )
        if ind is not None:
            knots_per_bs = knots_per_bs[:, ind]

    if return_cents is True:

        cents_per_bs = _utils_bsplines._get_cents_per_bs(
            cents, deg=deg, returnas=returnas,
        )
        if ind is not None:
            cents_per_bs = cents_per_bs[:, ind]

    # -------------
    # return

    if return_knots is True and return_cents is True:
        out = (knots_per_bs, cents_per_bs)
    elif return_knots is True:
        out = knots_per_bs
    else:
        out = cents_per_bs
    return out


def _mesh2DRect_bsplines_knotscents(
    returnas=None,
    return_knots=None,
    return_cents=None,
    ind=None,
    deg=None,
    Rknots=None,
    Zknots=None,
    Rcents=None,
    Zcents=None,
):

    # -------------
    # check inputs

    return_knots = ds._generic_check._check_var(
        return_knots, 'return_knots',
        types=bool,
        default=True,
    )
    return_cents = ds._generic_check._check_var(
        return_cents, 'return_cents',
        types=bool,
        default=True,
    )
    if return_knots is False and return_cents is False:
        return

    # -------------
    # compute

    if return_knots is True:

        knots_per_bs_R = _utils_bsplines._get_knots_per_bs(
            Rknots, deg=deg, returnas=returnas,
        )
        knots_per_bs_Z = _utils_bsplines._get_knots_per_bs(
            Zknots, deg=deg, returnas=returnas,
        )
        if ind is not None:
            knots_per_bs_R = knots_per_bs_R[:, ind[0]]
            knots_per_bs_Z = knots_per_bs_Z[:, ind[1]]

        nknots = knots_per_bs_R.shape[0]
        knots_per_bs_R = np.tile(knots_per_bs_R, (nknots, 1))
        knots_per_bs_Z = np.repeat(knots_per_bs_Z, nknots, axis=0)

    if return_cents is True:

        cents_per_bs_R = _utils_bsplines._get_cents_per_bs(
            Rcents, deg=deg, returnas=returnas,
        )
        cents_per_bs_Z = _utils_bsplines._get_cents_per_bs(
            Zcents, deg=deg, returnas=returnas,
        )
        if ind is not None:
            cents_per_bs_R = cents_per_bs_R[:, ind[0]]
            cents_per_bs_Z = cents_per_bs_Z[:, ind[1]]

        ncents = cents_per_bs_R.shape[0]
        cents_per_bs_R = np.tile(cents_per_bs_R, (ncents, 1))
        cents_per_bs_Z = np.repeat(cents_per_bs_Z, ncents, axis=0)

    # -------------
    # return

    if return_knots is True and return_cents is True:
        out = (
            (knots_per_bs_R, knots_per_bs_Z), (cents_per_bs_R, cents_per_bs_Z)
        )
    elif return_knots is True:
        out = (knots_per_bs_R, knots_per_bs_Z)
    else:
        out = (cents_per_bs_R, cents_per_bs_Z)
    return out