# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 18:15:16 2025

@author: dvezinet
"""


import numpy as np
import datastock as ds


# #############################################################
# #############################################################
#                  Main
# #############################################################


def main(
    coll=None,
    # 3 componants
    key_XR=None,
    key_YZ=None,
    key_Zphi=None,
    # linear vs toroidal
    geometry=None,
    # starting points
    pts_X=None,
    pts_Y=None,
    pts_Z=None,
):

    # ------------------
    # check inputs
    # ------------------

    (
        dkey,
        dpts,
        geometry,
    ) = _check(
        coll=coll,
        # 3 componants
        key_XR=key_XR,
        key_YZ=key_YZ,
        key_Zphi=key_Zphi,
        # linear vs toroidal
        geometry=geometry,
        # starting points
        pts_X=pts_X,
        pts_Y=pts_Y,
        pts_Z=pts_Z,
    )

    # ------------------
    # compute
    # ------------------

    # -----------------------
    # loop on starting points

    if geometry == 'toroidal':

        for ii, ind in enumerate(zip(*dpts['iok'].nonzero())):

            _compute_toroidal(
                coll=coll,
                dkey=dkey,
                ptsX=dpts['pts_X'][ind],
                ptsY=dpts['pts_Y'][ind],
                ptsZ=dpts['pts_Z'][ind],
            )

    # ------------------
    # format output
    # ------------------

    dout = {}

    return dout


# #############################################################
# #############################################################
#                  Check
# #############################################################


def _check(
    coll=None,
    # 3 componants
    key_XR=None,
    key_YZ=None,
    key_Zphi=None,
    # linear vs toroidal
    geometry=None,
    # starting points
    pts_X=None,
    pts_Y=None,
    pts_Z=None,
):

    # ------------------
    # key coordinates
    # ------------------

    dkey = _check_keys_components(
        coll=coll,
        # 3 componants
        key_XR=key_XR,
        key_YZ=key_YZ,
        key_Zphi=key_Zphi,
    )

    # ------------------
    # geometry
    # ------------------

    geometry = ds._generic_check._check_var(
        geometry, 'geometry',
        types=str,
        default='toroidal',
        allowed=['toroidal', 'linear'],
    )

    # not implemented
    if geometry == 'linear':
        msg = "Linear field line tracing not implemented yet!"
        raise NotImplementedError(msg)

    # ------------------
    # starting points
    # ------------------

    dpts = _check_pts(
        pts_X=pts_X,
        pts_Y=pts_Y,
        pts_Z=pts_Z,
    )

    return (
        dkey,
        dpts,
        geometry,
    )


# #############################################################
# #############################################################
#                  Check components
# #############################################################


def _check_keys_components(
    coll=None,
    # 3 componants
    key_XR=None,
    key_YZ=None,
    key_Zphi=None,
):

    # --------------
    # prepare
    # --------------

    wm = coll._which_mesh
    wbs = coll._which_bsplines

    # check valid 2d bsplines exist
    lok_bs = [
        k0 for k0, v0 in coll.dobj.get(wbs, {}).items()
        if coll.dobj[wm][v0[wm]]['nd'] == '2d'
    ]
    if len(lok_bs) == 0:
        msg = (
            "Line tracing for vector field: no valid bsplines2d!\n\n"
            + coll.show(wbs, verb=False, returnas=str)
        )
        raise Exception(msg)

    # check valid data defined on 2d bspline exist
    lok_data = [
        k0 for k0, v0 in coll.ddata.items()
        if v0.get(wbs) is not None
        and any([bs in lok_bs for bs in v0[wbs]])
    ]
    if len(lok_data) == 0:
        msg = (
            "Line tracing for vector field:"
            " no valid data defined on 2d bsplines\n"
            f"\t- lok_bs = {lok_bs}\n"
        )
        raise Exception(msg)

    # -------------
    # initialize
    # -------------

    dkey = {'key_XR': key_XR, 'key_YZ': key_YZ, 'key_Zphi': key_Zphi}
    dbs = {}
    dref = {}

    # -----------------
    # check component by component
    # -----------------

    for k0, v0 in dkey.items():

        c0 = v0 in lok_data and len(coll.ddata[v0][wbs]) == 1
        if not c0:
            msg = (
                "Line tracing for vector field wrong component '{k0}'!\n"
                f"\t- Not dependent on a (single) valid 2d bsplines!\n"
                f"\t- depends on {v0[wbs]}\n"
                f"\t- Avilable valid 2d bsplines: {lok_bs}\n"
            )
            raise Exception(msg)

        dbs[k0] = coll.ddata[v0][wbs][0]
        dref[k0] = coll.ddata[v0]['ref']

    # -----------------
    # check uniformity
    # -----------------

    # bsplines
    lbs = list(set(dbs.values()))
    if len(lbs) != 1:
        lstr = [f"\t- {k0}: {v0}" for k0, v0 in dbs.items()]
        msg = (
            "Line tracing for vector field: non-uniform bsplines 2d!\n"
            + "\n".join(lstr)
        )
        raise Exception(msg)

    # ref
    lref = list(set(dref.values()))
    if len(lref) != 1:
        lstr = [f"\t- {k0}: {v0}" for k0, v0 in dref.items()]
        msg = (
            "Line tracing for vector field: non-uniform ref!\n"
            + "\n".join(lstr)
        )
        raise Exception(msg)

    # -------------------------
    # common bsplines and mesh

    dkey.update({
        'key_mesh': coll.dobj[wbs][lbs[0]][wm],
        'key_bs': lbs[0],
        'ref': lref[0],
    })

    return dkey


# #############################################################
# #############################################################
#                  Check pts
# #############################################################


def _check_pts(
    pts_X=None,
    pts_Y=None,
    pts_Z=None,
):

    # -------------
    # initialize
    # -------------

    dpts = {'pts_X': pts_X, 'pts_Y': pts_Y, 'pts_Z': pts_Z}

    # -------------
    # type and finite
    # -------------

    for k0, v0 in dpts.items():

        # convert
        if np.isscalar(v0):
            dpts[k0] = np.r_[v0]
        elif isinstance(v0, (list, tuple)):
            dpts[k0] = np.array(v0)

        # check
        if not isinstance(dpts[k0], np.ndarray):
            msg = (
                f"Line tracing for vector field, arg '{k0}' wrong type!\n"
                "\t- expected: np.ndarray\n"
                f"\t- Provided: {type(dpts[k0])}\n"
            )
            raise Exception(msg)

        # should have at least one finite value
        if not np.any(np.isfinite(dpts[k0])):
            msg = (
            )
            raise Exception(msg)

    # -------------
    # uniformity
    # -------------

    # shape
    dshape = {k0: v0.shape for k0, v0 in dpts.items()}

    if len(set(dshape.values())) != 1:
        lstr = [f"\t- {k0}: {v0}" for k0, v0 in dshape.items()]
        msg = (
            "Line tracing for vector field, "
            "pts coordinates have different shape!\n"
            + "\n".join(lstr)
        )
        raise Exception(msg)

    # iok
    iok = np.all(
        [np.isfinite(v0) for v0 in dpts.values()],
        axis=0,
    )

    if not np.any(iok):
        msg = (
            "Not a single common finite value found in pts_X, pts_Y, ptsZ!\n"
        )
        raise Exception(msg)

    dpts['iok'] = iok

    return dpts


# #############################################################
# #############################################################
#                  Toroidal
# #############################################################


def _compute_toroidal(
    coll=None,
    dkey=None,
    ptsX=None,
    ptsY=None,
    ptsZ=None,
):

    # ---------------
    #
    # ---------------

    # find dx, dy, dz /
    #   (dx, dy, dz) cross (FX, FY, FZ) = 0
    #   dl small


    return