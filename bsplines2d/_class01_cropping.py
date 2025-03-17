# -*- coding: utf-8 -*-


# Built-in


# Common
import numpy as np
from matplotlib.path import Path
import datastock as ds


# specific
from . import _class01_select as _select


# #############################################################################
# #############################################################################
#                          crop rect mesh
# #############################################################################


def crop(
    coll=None,
    key=None,
    crop=None,
    thresh_in=None,
    remove_isolated=None,
):
    """ Crop a rect mesh

    Parameters
    ----------
    key:        str
        key of the rect mesh to be cropped
    crop:      np.ndarray
        Can be either:
            - bool: a boolean mask array
            - float: a closed 2d polygon used for cropping
    threshin:   int
        minimum nb. of corners for a mesh element to be included
    remove_isolated: bool
        flag indicating whether to remove isolated mesh elements

    Return
    ------
    crop:       np.ndarray
        bool mask
    key:        str
        key of the rect mesh to be cropped
    thresh_in:  int
        minimum nb. of corners for a mesh element to be included

    """

    # ------------
    # check inputs
    # ------------

    key, mtype, cropbool, thresh_in, remove_isolated = _crop_check(
        coll=coll, key=key, crop=crop, thresh_in=thresh_in,
        remove_isolated=remove_isolated,
    )

    # -----------
    # if crop is a poly => compute as bool
    # -----------

    if not cropbool:

        # --------------
        # rect

        if mtype == 'rect':

            (Rc, Zc), (Rk, Zk) = coll.select_mesh_elements(
                key=key,
                elements='cents',
                return_neighbours=True,
                returnas='data',
            )

            nR, nZ = Rc.shape
            npts = Rk.shape[-1] + 1

            pts = np.concatenate(
                (
                    np.concatenate((Rc[:, :, None], Rk), axis=-1)[..., None],
                    np.concatenate((Zc[:, :, None], Zk), axis=-1)[..., None],
                ),
                axis=-1,
            ).reshape((npts*nR*nZ, 2))

            isin = Path(crop.T).contains_points(pts).reshape((nR, nZ, npts))
            crop = np.sum(isin, axis=-1) >= thresh_in

            # Remove isolated pixelsi
            if remove_isolated is True:
                # All pixels should have at least one neighbour in R and one in Z
                # This constraint is useful for discrete gradient evaluation (D1N2)
                crop0 = crop
                while True:

                    # neighR
                    neighR = np.copy(crop0)
                    neighR[0, :] &= neighR[1, :]
                    neighR[-1, :] &= neighR[-2, :]
                    neighR[1:-1, :] &= (neighR[:-2, :] | neighR[2:, :])

                    # neighZ
                    neighZ = np.copy(crop0)
                    neighZ[:, 0] &= neighZ[:, 1]
                    neighZ[:, -1] &= neighZ[:, -2]
                    neighZ[:, 1:-1] &= (neighZ[:, :-2] | neighZ[:, 2:])

                    # overall
                    crop = neighR & neighZ

                    # stop or continue
                    if np.all(crop[crop0]):
                        break
                    else:
                        crop0 = crop

        # --------------
        # tri

        else:

            Rc, Zc = coll.select_mesh_elements(
                key=key,
                elements='cents',
                return_neighbours=False,
                returnas='data',
            )

            pts = np.array([Rc, Zc]).T
            crop = Path(crop.T).contains_points(pts)

    return crop, key, thresh_in


# #############################################################################
# #############################################################################
#                           check
# #############################################################################


def _crop_check(
    coll=None,
    key=None,
    crop=None,
    thresh_in=None,
    remove_isolated=None,
):

    # -----------
    # key
    # -----------

    wm = coll._which_mesh
    lkm = list(coll.dobj[wm].keys())
    key = ds._generic_check._check_var(
        key, 'key',
        default=None,
        types=str,
        allowed=lkm,
    )
    meshtype = coll.dobj[wm][key]['type']

    if meshtype not in ['rect', 'tri']:
        raise NotImplementedError()

    # shape
    shape = coll.dobj[wm][key]['shape_c']

    # -----------
    # crop
    # -----------

    # basic
    c0 = (
        isinstance(crop, np.ndarray)
        and crop.ndim == 2
        and np.all(np.isfinite(crop))
        and (
            (
                crop.shape[0] == 2
                and np.allclose(crop[:, 0], crop[:, -1])
                and crop.dtype in [int, float]
            )
            or (
                crop.shape == shape
                and crop.dtype == bool
            )
        )
    )
    if not c0:
        msg = (
            "Arg crop must be either:\n"
            f"\t- array of bool: mask of shape {shape}\n"
            f"\t- array of floats: (2, npts) closed (R, Z) polygon\n"
            f"Provided:\n{crop}"
        )
        raise Exception(msg)

    cropbool = crop.dtype == np.bool_

    # thresh_in and maxth
    if thresh_in is None:
        thresh_in = 3
    maxth = 5 if coll.dobj[wm][key]['type'] == 'rect' else 4

    c0 = (
        isinstance(thresh_in, (int, np.integer))
        and (1 <= thresh_in <= maxth)
    )
    if not c0:
        msg = (
            f"Arg thresh_in must be a int in 1 <= thresh_in <= {maxth}\n"
            f"Provided: {thresh_in}"
        )
        raise Exception(msg)

    # ----------------
    # remove_isolated
    # ----------------

    remove_isolated = ds._generic_check._check_var(
        remove_isolated, 'remove_isolated',
        default=True,
        types=bool,
    )

    return key, meshtype, cropbool, thresh_in, remove_isolated


# #############################################################################
# #############################################################################
#                           crop bsplines
# #############################################################################


def add_cropbs_from_crop(coll=None, keybs=None, keym=None):

    # ----------------
    # get

    kcropbs = False
    wm = coll._which_mesh
    wbs = coll._which_bsplines

    if coll.dobj[wm][keym]['crop'] is not False:
        kcropm = coll.dobj[wm][keym]['crop']
        cropbs = _get_cropbs_from_crop(
            coll=coll,
            crop=coll.ddata[kcropm]['data'],
            keybs=keybs,
        )
        kcropbs = f'{keybs}_crop'
        kcroppedbs = f'{keybs}_nbs_crop'

    elif coll.dobj[wm][keym]['type'] == 'tri':
        kcropbs = f'{keybs}_crop'

    # ----------------
    # optional crop

    if kcropbs is not False:

        # add cropped flat reference
        coll.add_ref(
            key=kcroppedbs,
            size=int(cropbs.sum()),
        )

        coll.add_data(
            key=kcropbs,
            data=cropbs,
            ref=coll._dobj[wbs][keybs]['ref'],
            dim='bool',
            quant='bool',
        )
        coll._dobj[wbs][keybs]['crop'] = kcropbs


# #############################################################################
# #############################################################################
#                           Mesh2DRect - crop
# #############################################################################


def _get_cropbs_from_crop(coll=None, crop=None, keybs=None):

    if isinstance(crop, str) and crop in coll.ddata.keys():
        crop = coll.ddata[crop]['data']

    wm = coll._which_mesh
    wbs = coll._which_bsplines

    shref = coll.dobj[wm][coll.dobj[wbs][keybs][wm]]['shape_c']
    if crop.shape != shref:
        msg = "Arg crop seems to have the wrong shape!"
        raise Exception(msg)

    keym = coll.dobj[wbs][keybs][coll._which_mesh]
    kRk, kZk = coll.dobj[wm][keym]['knots']
    kRc, kZc = coll.dobj[wm][keym]['cents']

    cents_per_bs_R, cents_per_bs_Z = _select._mesh2DRect_bsplines_knotscents(
        returnas='ind',
        return_knots=False,
        return_cents=True,
        ind=None,
        deg=coll.dobj[wbs][keybs]['deg'],
        Rknots=coll.ddata[kRk]['data'],
        Zknots=coll.ddata[kZk]['data'],
        Rcents=coll.ddata[kRc]['data'],
        Zcents=coll.ddata[kZc]['data'],
    )

    shapebs = coll.dobj[wbs][keybs]['shape']
    cropbs = np.array([
        [
            np.all(crop[cents_per_bs_R[:, ii], cents_per_bs_Z[:, jj]])
            for jj in range(shapebs[1])
        ]
        for ii in range(shapebs[0])
    ], dtype=bool)

    return cropbs