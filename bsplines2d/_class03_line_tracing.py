# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 18:15:16 2025

@author: dvezinet
"""


import copy


import numpy as np
import scipy.integrate as scpinteg
import datastock as ds


# #############################################################
# #############################################################
#                  DEFAULT
# #############################################################


_DSOLVER = {
    'method': 'RK45',
    'first_step': None,
    'min_step': 0.,
    'max_step': np.inf,
    'rtol': 1e-6,
    'atol': 1e-6,
    'dense_output': False,
    'vectorized': False,   # only for ‘Radau’ and ‘BDF’
    'fac': None,           # only for ‘Radau’, ‘BDF’ and ‘LSODA’
}


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
    # domain
    domain=None,
    # linear vs toroidal
    geometry=None,
    # starting points
    pts_X=None,
    pts_Y=None,
    pts_Z=None,
    # length and resolution of lines
    res=None,
    length=None,
    npts=None,
    dl=None,
    # direction
    direct=None,
    # solver options
    dsolver=None,
):

    # ------------------
    # check inputs
    # ------------------

    (
        dkey,
        dpts,
        geometry,
        dl,
        direct,
        dsolver,
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
        # length and resolution of lines
        res=res,
        length=length,
        npts=npts,
        dl=dl,
        # direction
        direct=direct,
        # solver
        dsolver=dsolver,
    )

    # ------------------
    # prepare output
    # ------------------

    dout, update_sli_dom, shape_key_nobs = _prepare_output(
        coll=coll,
        dkey=dkey,
        dpts=dpts,
        dl=dl,
        direct=direct,
    )

    # ------------------
    # get derivative functions
    # ------------------

    func_fwd, func_bck = _get_func(
        coll=coll,
        dkey=dkey,
        geometry=geometry,
    )

    # ------------------
    # compute
    # ------------------

    # -----------------------
    # loop on starting points

    for i0, ipts in enumerate(zip(*dpts['iok'].nonzero())):

        # ------------------
        # loop on field refs

        for i1, iref in enumerate(np.ndindex(shape_key_nobs)):

            # update domain
            (
                sli_fwd, sli_bck,
                domaini,
                sli_ev,
            ) = update_sli_dom(ipts, iref)

            # -------
            # Forward

            if 'fwd' in direct:

                sol = scpinteg.solve_ivp(
                    func_fwd,
                    [dl[0], dl[-1]],
                    [pts_X[ipts], pts_Y[ipts], pts_Z[ipts]],
                    args=(domaini,),
                    t_eval=dl,
                    **dsolver,
                )

                dout['nfev_fwd']['data'][sli_ev] = sol.nfev
                dout['status_fwd']['data'][sli_ev] = sol.status

                if sol.success is True:
                    dout['pts_X']['data'][sli_fwd] = sol.y[0, :]
                    dout['pts_Y']['data'][sli_fwd] = sol.y[1, :]
                    dout['pts_Z']['data'][sli_fwd] = sol.y[2, :]

            # --------
            # Backward

            if 'bck' in direct:

                sol = scpinteg.solve_ivp(
                    func_bck,
                    [dl[0], dl[-1]],
                    [pts_X[ipts], pts_Y[ipts], pts_Z[ipts]],
                    args=(domaini,),
                    t_eval=dl,
                    **dsolver,
                )

                dout['nfev_bck']['data'][sli_ev] = sol.nfev
                dout['status_bck']['data'][sli_ev] = sol.status

                if sol.success is True:
                    dout['pts_X']['data'][sli_bck] = sol.y[0, :]
                    dout['pts_Y']['data'][sli_bck] = sol.y[1, :]
                    dout['pts_Z']['data'][sli_bck] = sol.y[2, :]

    # ------------------
    # return
    # ------------------

    if geometry == 'toroidal':
        _add_dphi(
            dout,
            direct=direct,
            npts=dl.size,
        )

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
    # length and resolution of lines
    res=None,
    length=None,
    npts=None,
    dl=None,
    # direction
    direct=None,
    # solver
    dsolver=None,
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

    # ------------------
    # starting points
    # ------------------

    dpts = _check_pts(
        pts_X=pts_X,
        pts_Y=pts_Y,
        pts_Z=pts_Z,
    )

    # ------------------
    # resolution and length
    # ------------------

    dl = _check_res_length(
        res=res,
        length=length,
        npts=npts,
        dl=dl,
    )

    # ------------------
    # direction
    # ------------------

    if isinstance(direct, str):
        direct = [direct]

    lok = ['fwd', 'bck']
    direct = ds._generic_check._check_var_iter(
        direct, 'direct',
        types=(list, tuple),
        types_iter=str,
        default=lok,
        allowed=lok,
    )

    # ------------------
    # solver
    # ------------------

    dsolver = _check_dsolver(dsolver)

    return (
        dkey,
        dpts,
        geometry,
        dl,
        direct,
        dsolver,
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
        'ref_bs': coll.dobj[wbs][lbs[0]]['ref'],
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
#                  Check resolution, length
# #############################################################


def _check_res_length(
    res=None,
    length=None,
    npts=None,
    dl=None,
):

    # ------------------
    # res vs length vs npts
    # ------------------

    lc = [
        res is not None and length is not None,
        res is not None and npts is not None,
        npts is not None and length is not None,
        dl is not None,
    ]
    if np.sum(lc) != 1:
        lnames = [
            ('res', 'length'), ('res', 'npts'),
            ('npts', 'length'), ('dl', '')
        ]
        lstr = [
            f"\t- '{k0}' and '{k1}':\t{lc[ii]}"
            for ii, (k0, k1) in enumerate(lnames)
        ]
        msg = (
            "Field line tracing: provide one only of the following input:\n"
            + "\n".join(lstr)
        )
        raise Exception(msg)

    # ------------------
    # check each
    # ------------------

    # res
    if res is not None:
        res = float(ds._generic_check._check_var(
            res, 'res',
            types=(int, float),
            sign='>0.',
        ))

    # npts
    if npts is not None:
        npts = int(ds._generic_check._check_var(
            npts, 'npts',
            types=(int, float),
            sign='>0.',
        ))

    # length
    if length is not None:
        length = float(ds._generic_check._check_var(
            length, 'length',
            types=(int, float),
            sign='>0.',
        ))

    # dl
    if dl is not None:
        dl = ds._generic_check._check_flat1darray(
            dl, 'dl',
            dtypes=float,
            unique=True,
        )

    # ------------------
    # case by case
    # ------------------

    if lc[0]:
        dl = np.linspace(0, length, int(np.ceil(length / res)))

    elif lc[1]:
        dl = np.linspace(0, res * (npts - 1), npts)

    elif lc[2]:
        dl = np.linspace(0, length, npts)

    return dl


# #############################################################
# #############################################################
#              check dsolver
# #############################################################


def _check_dsolver(dsolver=None):

    # ---------------
    # Default
    # ---------------

    ddef = copy.deepcopy(_DSOLVER)

    if dsolver is None:
        dsolver = ddef

    # ---------------
    # check
    # ---------------

    if not isinstance(dsolver, dict):
        _err_dsolver(dsolver, ddef)

    for k0, v0 in ddef.items():
        if dsolver.get(k0) is None:
            dsolver[k0] = v0

    if any([k0 not in ddef.keys() for k0 in dsolver.keys()]):
        _err_dsolver(dsolver, ddef)

    return dsolver


def _err_dsolver(dsolver, ddef):
    lstr = [f"\t- {k0}: {v0} (default)" for k0, v0 in ddef.items()]
    msg = (
        "Arg dsolver must be a dict with some of the following:\n"
        + "\n".join(lstr)
        + "\nIs fed to scipy.integrate.solve_ivp()"
    )

    raise Exception(msg)


# #############################################################
# #############################################################
#              prepare output
# #############################################################


def _prepare_output(
    coll=None,
    dkey=None,
    dpts=None,
    dl=None,
    direct=None,
):

    # -------------
    # npts_full
    # -------------

    # npts_full vs direct
    npts = dl.size
    if len(direct) == 1:
        npts_full = npts
    else:
        npts_full = 2*npts - 1

    # ref_npts
    ref_npts = True

    # -------------
    # ref and shape
    # -------------

    # ref
    ref_key = dkey['ref']
    ref_bs = dkey['ref_bs']
    ref_key_nobs = tuple([rr for rr in ref_key if rr not in ref_bs])

    ibs = ref_key.index(ref_bs[0])
    ref_out = (
        ref_key[:ibs]
        + (ref_npts,) + tuple([None for ss in dpts['pts_X'].shape])
        + ref_key[ibs+2:]
    )

    ref_ev = (
        ref_key[:ibs]
        + tuple([None for ss in dpts['pts_X'].shape])
        + ref_key[ibs+2:]
    )

    # shape
    shape_key = coll.ddata[dkey['key_XR']]['shape']
    shape_key_nobs = tuple([
        shape_key[ii] for ii, rr in enumerate(ref_key)
        if rr not in ref_bs
    ])

    shape_out = (
        shape_key[:ibs]
        + (npts_full,) + dpts['pts_X'].shape
        + shape_key[ibs+2:]
    )

    shape_ev = (
        shape_key[:ibs]
        + dpts['pts_X'].shape
        + shape_key[ibs+2:]
    )

    # nan-filled
    nan = np.full(shape_out, np.nan)
    nanev = np.full(shape_ev, np.nan)

    # ------------------
    # slices and domain
    # ------------------

    update_sli_dom = _get_update_sli_dom(
        ref_out=ref_out,
        ref_key_nobs=ref_key_nobs,
        direct=direct,
        npts=npts,
    )

    # -------------
    # prepare dl_full
    # -------------

    if len(direct) == 1:
        if 'bck' in direct:
            dl_full = -dl
        else:
            dl_full = dl

    else:
        dl_full = np.r_[-dl[:0:-1], dl]
    assert dl_full.size == npts_full

    # -------------
    # dout
    # -------------

    dout = {
        'pts_X': {
            'data': np.copy(nan),
            'units': 'm',
            'ref': ref_out,
        },
        'pts_Y': {
            'data': np.copy(nan),
            'units': 'm',
            'ref': ref_out,
        },
        'pts_Z': {
            'data': np.copy(nan),
            'units': 'm',
            'ref': ref_out,
        },
        'dl': {
            'data': dl_full,
            'ref': (True,),
            'units': 'm',
            'dim': 'length',
        },
    }

    # forward solver info
    if 'fwd' in direct:
        dout.update({
            'nfev_fwd': {
                'data': np.copy(nanev),
                'units': None,
                'ref': ref_ev,
            },
            'status_fwd': {
                'data': np.copy(nanev),
                'units': None,
                'ref': ref_ev,
            },
        })

    # backward solver info
    if 'bck' in direct:
        dout.update({
            'nfev_bck': {
                'data': np.copy(nanev),
                'units': None,
                'ref': ref_ev,
            },
            'status_bck': {
                'data': np.copy(nanev),
                'units': None,
                'ref': ref_ev,
            },
        })

    return dout, update_sli_dom, shape_key_nobs


def _get_update_sli_dom(
    ref_out=None,
    ref_key_nobs=None,
    direct=None,
    npts=None,
):

    # -------------------
    # prepare ifwd, ibck
    # -------------------

    if len(direct) == 2:
        ifwd = slice(npts-1, None, 1)
        ibck = slice(npts-1, None, -1)

    else:
        if 'fwd' in direct:
            ifwd = slice(None)
            ibck = None

        else:
            ibck = slice(None)
            ifwd = None

    # --------------
    # initialize
    # --------------

    # -------
    # slices

    sli_fwd = [ifwd if rr is True else 0 for rr in ref_out]
    sli_bck = [ibck if rr is True else 0 for rr in ref_out]

    sli_ev = [0 for rr in ref_out[:-1]]

    # ----------
    # indr, indp

    indr = [ii for ii, rr in enumerate(ref_out) if rr in ref_key_nobs]
    indp = [ii for ii, rr in enumerate(ref_out) if rr is None]

    ref_ev = [rr for rr in ref_out if rr is not True]
    indr_ev = [ii for ii, rr in enumerate(ref_ev) if rr in ref_key_nobs]
    indp_ev = [ii for ii, rr in enumerate(ref_ev) if rr is None]

    # -------
    # domain

    domain = {k0: {'ind': None} for k0 in ref_key_nobs}

    # -------------
    # update func
    # -------------

    def func(
        ipts, iref,
        sli_fwd=np.array(sli_fwd),
        sli_bck=np.array(sli_bck),
        sli_ev=np.array(sli_ev),
        indr=np.array(indr, dtype=int),
        indp=np.array(indp, dtype=int),
        indr_ev=np.array(indr_ev, dtype=int),
        indp_ev=np.array(indp_ev, dtype=int),
        domain=domain,
        ref_key_nobs=ref_key_nobs,
    ):

        # --------
        # ev

        sli_ev[indr_ev] = iref
        sli_ev[indp_ev] = ipts

        # --------
        # pts

        if sli_fwd is not None:
            sli_fwd[indr] = iref
            sli_fwd[indp] = ipts

        if sli_bck is not None:
            sli_bck[indr] = iref
            sli_bck[indp] = ipts

        # --------
        # domain

        for ii, k0 in enumerate(ref_key_nobs):
            domain[k0]['ind'] = iref[ii]

        return (
            tuple(sli_fwd), tuple(sli_bck),
            domain,
            tuple(sli_ev),
        )

    return func


# #############################################################
# #############################################################
#              Get derivative functions
# #############################################################


def _get_func(
    coll=None,
    dkey=None,
    geometry=None,
):

    if geometry == 'toroidal':

        func_fwd = _get_func_toroidal(
            coll=coll,
            dkey=dkey,
            sign=1.,
        )

        func_bck = _get_func_toroidal(
            coll=coll,
            dkey=dkey,
            sign=-1.,
        )

    else:

        func_fwd = _get_func_linear(
            coll=coll,
            dkey=dkey,
            sign=1.,
        )

        func_bck = _get_func_linear(
            coll=coll,
            dkey=dkey,
            sign=-1.,
        )

    return func_fwd, func_bck


def _get_func_toroidal(coll=None, dkey=None, sign=1.):

    def func(ll, pxyz, domain, coll=coll, sign=sign, dkey=dkey):

        R = np.hypot(pxyz[0], pxyz[1])

        dout = coll.interpolate(
            keys=[dkey['key_XR'], dkey['key_YZ'], dkey['key_Zphi']],
            ref_key=dkey['key_bs'],
            x0=R,
            x1=pxyz[2],
            grid=False,
            domain=domain,
            details=False,
            store=False,
        )

        # extract and normalized
        BR = dout[dkey['key_XR']]['data']
        BZ = dout[dkey['key_YZ']]['data']
        Bphi = dout[dkey['key_Zphi']]['data']
        BN_sign_inv = sign / np.sqrt(BR**2 + BZ**2 + Bphi**2)

        cosphi = pxyz[0] / R
        sinphi = pxyz[1] / R

        BXn = (BR * cosphi - Bphi * sinphi) * BN_sign_inv
        BYn = (BR * sinphi + Bphi * cosphi) * BN_sign_inv
        BZn = BZ * BN_sign_inv

        return [BXn[0], BYn[0], BZn[0]]

    return func


def _get_func_linear(coll=None, dkey=None, sign=1.):

    def func(ll, pxyz, domain, coll=coll, sign=sign, dkey=dkey):

        dout = coll.interpolate(
            keys=[dkey['key_XR'], dkey['key_YZ'], dkey['key_Zphi']],
            ref_key=dkey['key_bs'],
            x0=pxyz[0],
            x1=pxyz[1],
            grid=False,
            domain=domain,
            details=False,
            store=False,
        )

        # extract and normalized
        BX = dout[dkey['key_XR']]['data']
        BY = dout[dkey['key_YZ']]['data']
        BZ = dout[dkey['key_Zphi']]['data']
        BN_sign_inv = sign / np.sqrt(BX**2 + BY**2 + BZ**2)

        BXn = BX * BN_sign_inv
        BYn = BY * BN_sign_inv
        BZn = BZ * BN_sign_inv

        return [BXn[0], BYn[0], BZn[0]]

    return func


# #############################################################
# #############################################################
#              add dphi
# #############################################################


def _add_dphi(
    dout=None,
    direct=None,
    npts=None,
):

    # ------------
    # axis
    # ------------

    axis = dout['pts_X']['ref'].index(True)

    # ------------
    # phi
    # ------------

    phi = np.arctan2(dout['pts_Y']['data'], dout['pts_X']['data'])

    # ------------
    # get nturns
    # ------------

    if len(direct) == 2:
        i0 = npts - 1

        sli_bck = tuple([
            slice(i0, None, -1) if ii == axis else slice(None)
            for ii in range(phi.ndim)
        ])
        sli_fwd = tuple([
            slice(i0, None, 1) if ii == axis else slice(None)
            for ii in range(phi.ndim)
        ])
        nturns_fwd = _add_turns(phi[sli_fwd], axis=axis)
        nturns_bck = _add_turns(phi[sli_bck], axis=axis)

        sli_rev = tuple([
            slice(None, 0, -1) if ii == axis else slice(None)
            for ii in range(phi.ndim)
        ])
        nturns = np.concatenate(
            (nturns_bck[sli_rev], nturns_fwd),
            axis=axis,
        )

    else:
        i0 = 0
        nturns = _add_turns(phi, axis=axis)

    # -----------------
    # derive dphi
    # -----------------

    sli0 = tuple([
        slice(i0, i0+1) if ii == axis else slice(None)
        for ii in range(phi.ndim)
    ])
    dphi = (phi - phi[sli0]) + nturns * 2*np.pi

    # ------------
    # store
    # ------------

    dout.update({
        'dphi': {
            'data': dphi,
            'units': 'rad',
            'ref': dout['pts_X']['ref'],
            'dim': 'angle',
        }
    })

    return


def _add_turns(phi, axis=None):

    dphi = np.diff(phi, axis=axis)

    sli0 = tuple([
        slice(0, 1) if ii == axis else slice(None)
        for ii in range(dphi.ndim)
    ])
    sign = np.sign(dphi[sli0])

    ind = np.concatenate(
        (np.zeros(sign.shape), dphi * sign < 0),
        axis=axis,
    )

    return sign * np.cumsum(ind, axis=axis)