# -*- coding: utf-8 -*-


# Built-in
import itertools as itt
import warnings

# Common
import numpy as np
from matplotlib.path import Path
import datastock as ds


# tofu
# from . import _class02_checks as _checks
from . import _utils_bsplines


# ################################################################
# ################################################################
#                                       Main
# ################################################################


def interpolate(
    coll=None,
    # interpolation base
    keys=None,
    ref_key=None,
    # interpolation pts
    x0=None,
    x1=None,
    grid=None,
    res=None,
    mode=None,
    submesh=None,
    # domain limitation
    domain=None,
    # common ref
    ref_com=None,
    ref_vector_strategy=None,
    # bsplines-specific
    details=None,
    indbs_tf=None,
    # rect-specific
    crop=None,
    # parameters
    deg=None,
    deriv=None,
    val_out=None,
    log_log=None,
    nan0=None,
    # store vs return
    returnas=None,
    return_params=None,
    store=None,
    store_keys=None,
    inplace=None,
    # debug or unit tests
    debug=None,
):
    """ Interpolate at desired points on desired data

    Interpolate quantities (keys) on coordinates (ref_keys)
    All provided keys should share the same refs
    They should have dimension 2 or less

    ref_keys should be a list of monotonous data keys to be used as coordinates
    It should have one element per dimension in refs

    The interpolation points are provided as np.ndarrays in each dimension
    They should all have the same shape except if grid = True

    deg is the degree of the interpolation

    It is an interpolation, not a smoothing, so it will pass through all points
    Uses scpinterp.InterpolatedUnivariateSpline() for 1d

    """

    # -------------
    # check inputs

    if debug is None:
        debug = False

    # keys
    (
        isbs, keys, ref_key,
        daxis, dunits, units_ref,
        details, ktemp,
    ) = _check_keys(
        coll=coll,
        keys=keys,
        ref_key=ref_key,
        only1d=False,
        details=details,
    )

    # -----------
    # isbs

    x0_str = None
    crop_path = None
    if isbs:

        # -----------------
        # prepare bsplines

        keybs, keym, mtype, ref_key = _prepare_bsplines(
            coll=coll,
            ref_key=ref_key,
        )

        if mtype == 'tri' and ref_com is not None:
            msg = (
                "Arg ref_com cannot be used with a triangular mesh!"
            )
            raise NotImplementedError(msg)

        # ------------------------------
        # check bsplines-specific params

        (
            val_out, crop, cropbs, indbs_tf, nbs, submesh,
        ) = _check_params_bsplines(
            coll=coll,
            details=details,
            val_out=val_out,
            crop=crop,
            mtype=mtype,
            keybs=keybs,
            keym=keym,
            indbs_tf=indbs_tf,
            submesh=submesh,
        )

        if details is True:
            ref_com, domain = None, None

        if mtype == 'rect' and crop is not None:
            doutline = coll.get_mesh_outline(keym)
            crop_path = Path(
                np.array([doutline['x0']['data'], doutline['x1']['data']]).T
            )

        # # no x0 + no res + deg = 1 bsplines => x0 = knots
        wm = coll._which_mesh
        wbs = coll._which_bsplines
        deg = coll.dobj[wbs][keybs]['deg']
        if x0 is None and res is None and deg == 1:
            subbs = coll.dobj[wm][keym].get('subbs')
            if (
                submesh is True
                and subbs is not None
                and coll.dobj[wbs][subbs]['deg'] == 1
            ):
                bsx0 = subbs
            elif submesh is False or subbs is None:
                bsx0 = keybs

            kx01 = coll.dobj[bsx0]['apex']
            x0 = coll.ddata[kx01[0]]['data']
            if len(kx01) == 2:
                x1 = coll.ddata[kx01[1]]['data']

        # ---------------
        # no x0 => mesh

        if x0 is None:
            if submesh is True:
                km = coll.dobj[wm][keym]['submesh']
            else:
                km = keym

            out = coll.get_sample_mesh(
                key=km,
                res=res,
                grid=True,
                mode=mode,
                x0=None,
                x1=None,
                Dx0=None,
                Dx1=None,
                imshow=False,
            )

            x0 = out['x0']['data']
            x1 = out.get('x1', {}).get('data')

        # ---------------
        # submesh

        if submesh is True:

            # submesh
            kd0 = coll.dobj[wm][keym]['subkey']
            kbs0 = coll.dobj[wm][keym]['subbs']

            # interpolate
            dout_temp, dparams_temp = coll.interpolate(
                # interpolation base
                keys=kd0,
                ref_key=kbs0,
                # interpolation pts
                x0=x0,
                x1=x1,
                grid=grid,
                submesh=False,
                # domain
                domain=domain,
                # common ref
                ref_com=None,
                # bsplines-specific
                details=False,
                indbs_tf=None,
                # rect-specific
                crop=crop,
                return_params=True,
                store=False,
            )

            # ------------------------------------------
            # handle ref_com, add ref and data if needed

            ref_com = _submesh_ref_com(
                coll=coll,
                kd0=kd0,
                keys=keys,
                ref_com=ref_com,
                # coordinates
                x0=x0,
            )

            # temporary storing for ref handling
            lr_add, ld_add = _submesh_addtemp(
                coll=coll,
                kd0=kd0,
                # interp resut
                dout_temp=dout_temp,
                # coordinates
                x0=x0,
            )

            # --------
            # substitue x0, x1

            x0_str = isinstance(x0, str)
            # if ref_com is None:
            # x0 = dout_temp[kd0[0]]['data']
            # if len(kd0) > 1:
            # x1 = dout_temp[kd0[1]]['data']
            # else:
            x0 = dout_temp[kd0[0]]['key']
            if len(kd0) > 1:
                x1 = dout_temp[kd0[1]]['key']

            if len(kd0) == 1:
                x1 = None

        else:
            lr_add, ld_add = None, None

    # ---------------
    # prepare x0, x1

    # params
    (
        deg, deriv,
        kx0, kx1, x0, x1, refx, dref_com,
        ddata, dout, dsh_other, sli_c, sli_x, sli_v,
        log_log, nan0, grid, ndim, xunique,
        returnas, return_params, store, inplace,
        domain,
    ) = ds._class1_interpolate._check(
        coll=coll,
        # interpolation base
        keys=keys,
        ref_key=ref_key,      # ddata keys
        # interpolation pts
        x0=x0,
        x1=x1,
        # useful for shapes
        daxis=daxis,
        dunits=dunits,
        domain=domain,
        # common ref
        ref_com=ref_com,
        # ref_vector_strategy=ref_vector_strategy,
        # parameters
        grid=grid,
        deg=None,
        deriv=deriv,
        log_log=None,
        nan0=nan0,
        # return vs store
        returnas=returnas,
        return_params=return_params,
        store=store,
        inplace=inplace,
        x0_str=x0_str,
    )

    # ---------------
    # interpolate

    if isbs:

        _interp(
            # resources
            coll=coll,
            keybs=keybs,
            keys=keys,
            # interpolation points
            x0=x0,
            x1=x1,
            # options
            val_out=val_out,
            deriv=deriv,
            indbs_tf=indbs_tf,
            # cropping
            crop=crop,
            cropbs=cropbs,
            crop_path=crop_path,
            # details
            details=details,
            # others
            dout=dout,
            mtype=mtype,
            ddata=ddata,
            daxis=daxis,
            sli_c=sli_c,
            sli_x=sli_x,
            sli_v=sli_v,
            dsh_other=dsh_other,
            dref_com=dref_com,
            lr_add=lr_add,
            ld_add=ld_add,
            nan0=nan0,
        )

    else:

        ds._class1_interpolate._interp(
            coll=coll,
            keys=keys,
            ref_key=ref_key,
            x0=x0,
            x1=x1,
            ndim=ndim,
            log_log=log_log,
            dout=dout,
            ddata=ddata,
            dsh_other=dsh_other,
            daxis=daxis,
            deg=deg,
            deriv=deriv,
            dref_com=dref_com,
            sli_x=sli_x,
            sli_c=sli_c,
            sli_v=sli_v,
            nan0=nan0,
        )

    # ------------------------------
    # cleanup if details

    if details is True:
        coll.remove_data(ktemp)

    # ------------------------------
    # adjust data and ref if xunique

    if xunique:
        # try:
        ds._class1_interpolate._xunique(dout, domain=domain)
        # except Exception as err:
        #     msg = (
        #         err.args[0]
        #         + "\n\n"
        #         f"keys = {keys}\n"
        #         f"ref_key = {ref_key}\n"
        #         f"x0 = {x0}\n"
        #         f"x1 = {x1}\n"
        #     )
        #     err.args = (msg,)
        #     raise err

    # ----------
    # store

    if store is True:
        coll2 = ds._class1_interpolate._store(
            coll=coll,
            dout=dout,
            inplace=inplace,
            store_keys=store_keys,
        )

    # -------
    # return

    if returnas is object:
        return coll2

    elif returnas is dict:

        if return_params is True:
            dparam = {
                'keys': keys,
                'keybs': keybs,
                'ref_key': ref_key,
                'deg': deg,
                'deriv': deriv,
                'log_log': log_log,
                'kx0': kx0,
                'kx1': kx1,
                'x0': x0,
                'x1': x1,
                'grid': grid,
                'refx': refx,
                'ref_com': ref_com,
                'dref_com': dref_com,
                'daxis': daxis,
                'dsh_other': dsh_other,
                'domain': domain,
                'details': details,
                'crop': crop,
                'cropbs': cropbs,
                'indbs_tf': indbs_tf,
                'submesh': submesh,
                'subbs': kbs0 if submesh else None,
            }

            # debug
            if debug is True:
                dparam['x0.shape'] = x0.shape

            return dout, dparam
        else:
            return dout


# ####################################
# ####################################
#       check
# ####################################


def _check_keys(
    coll=None,
    keys=None,
    ref_key=None,
    only1d=None,
    details=None,
):

    # -------------
    # details

    details = ds._generic_check._check_var(
        details, 'details',
        types=bool,
        default=False,
    )

    # -------------
    # only1d

    only1d = ds._generic_check._check_var(
        only1d, 'only1d',
        types=bool,
        default=True,
    )

    maxd = 1 if only1d else 2

    # ---------------
    # details is True

    wbs = coll._which_bsplines
    lbs = list(coll.dobj[wbs].keys())

    ktemp = None
    if details is True:
        ref_key = ds._generic_check._check_var(
            ref_key, 'ref_key',
            types=str,
            allowed=lbs,
        )

        # add temporary data
        ktemp = f'{ref_key}_details'
        coll.add_data(
            key=ktemp,
            data=np.ones(coll.dobj[wbs][ref_key]['shape'], dtype=float),
            ref=ref_key,
        )
        keys = ktemp

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
                f"Arg ref_key shall not more than {maxd} ref!\n"
                "And can only contain one single bsplines!\n"
                f"Provided: {ref_key}"
            )
            raise Exception(msg)

        # bs vs refs
        if lc[1]:
            ref_key = ref_key[0]
            lok_bs = [
                k0 for k0, v0 in coll.ddata.items()
                if v0[wbs] is not None and ref_key in v0[wbs]
            ]
            lok_nobs = []

        else:

            # check vs valid vectors
            for ii, rr in enumerate(ref_key):

                if rr in lref:
                    kwd = {'ref': rr}
                else:
                    kwd = {'key0': rr}

                hasref, hasvect, ref, ref_key[ii] = coll.get_ref_vector(
                    **kwd
                )[:4]

                if not (hasref and hasvect):
                    msg = (
                        f"Provided ref_key[{ii}] not a valid ref or ref vector!\n"
                        f"Provided: {rr}"
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
                if all([
                    coll.ddata[rr]['ref'][0] in v0['ref']
                    for rr in ref_key
                ])
            ]
            lok_bs = []

        if keys is None:
            keys = lok_nobs + lok_bs

    else:
        # binning only for non-bsplines or 1d bsplines
        lok_nobs = [
            k0 for k0, v0 in coll.ddata.items()
            if k0 not in dbs.keys()
        ]
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
        types=(list, tuple),
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
        lbs = [
            [
                bs for bs in coll.ddata[k0][wbs]
                if len(coll.dobj[wbs][bs]['ref']) <= maxd
            ]
            for k0 in keys
            if any([
                len(coll.dobj[wbs][bs]['ref']) <= maxd
                for bs in coll.ddata[k0][wbs]
            ])
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
            k0: [
                coll.ddata[k0]['ref'].index(rr)
                for rr in coll.dobj[wbs][ref_key]['ref']
            ]
            for k0 in keys
        }

        # units_ref
        wbs = coll._which_bsplines
        units_ref = coll.ddata[coll.dobj[wbs][ref_key]['apex'][0]]['units']

    # ref_key
    else:
        if ref_key is None:
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

        # units_ref
        units_ref = [coll.ddata[rr]['units'] for rr in ref_key]

    # dunits
    dunits = {k0: coll.ddata[k0]['units'] for k0 in keys}

    return isbs, keys, ref_key, daxis, dunits, units_ref, details, ktemp


def _check_params_bsplines(
    coll=None,
    details=None,
    val_out=None,
    crop=None,
    mtype=None,
    keybs=None,
    keym=None,
    indbs_tf=None,
    submesh=None,
):

    # -------------
    # val_out

    val_out = ds._generic_check._check_var(
        val_out, 'val_out',
        default=np.nan,
        allowed=[0., np.nan, False],
    )

    # ---------------
    # cropping ?

    wbs = coll._which_bsplines
    cropbs = coll.dobj[wbs][keybs]['crop']
    if cropbs is None:
        cropbs = False
    if cropbs is not False:
        cropbs = coll.ddata[cropbs]['data']
        nbs = cropbs.sum()
    else:
        nbs = np.prod(coll.dobj[wbs][keybs]['shape'])

    # -------------
    # crop

    lok = [False]
    if mtype == 'rect' and cropbs is not False:
        lok.append(True)
    else:
        crop = False

    crop = ds._generic_check._check_var(
        crop, 'crop',
        types=bool,
        default=lok[-1],
        allowed=lok,
    )

    # ---------------
    # indbs_tf

    if details is True:

        if mtype == 'rect':
            returnas = 'tuple-flat'
            # returnas = 'array-flat'
        else:
            returnas = int

        # compute validated indbs array with appropriate form
        indbs_tf = coll.select_ind(
            key=keybs,
            returnas=returnas,
            ind=indbs_tf,
            crop=crop,
        )

        if isinstance(indbs_tf, tuple):
            nbs = indbs_tf[0].size
        else:
            nbs = indbs_tf.size

    # ---------
    # submesh

    submesh = ds._generic_check._check_var(
        submesh, 'submesh',
        types=bool,
        default=False,
    )
    if coll.dobj[coll._which_mesh][keym]['subkey'] is None:
        submesh = False

    return val_out, crop, cropbs, indbs_tf, nbs, submesh


# ################################################################
# ################################################################
#               Handle ref_com for submesh
# ################################################################


def _submesh_ref_com(
    coll=None,
    kd0=None,
    keys=None,
    ref_com=None,
    # coordinates
    x0=None,
):

    # ---------------------
    # find possible matches

    ref0 = coll.ddata[kd0[0]]['ref']
    lrcom = [
        (rr, ref0.index(rr))
        for ii, rr in enumerate(ref0)
        if rr in list(itt.chain.from_iterable([
            coll.ddata[kk]['ref'] for kk in keys
        ]))
        and ii in [0, len(ref0) - 1]
    ]

    # ----------
    # unused options

    if ref_com is None:
        if len(lrcom) > 0:
            msg = (
                f"\nPossible common ref for data {keys} and subkey '{kd0}':\n"
                + "\n".join([f"\t- {rr}" for rr in lrcom])
                + "\nIf you wish to use one, specify with ref_com=..."
            )
            warnings.warn(msg)

    else:
        # --------------
        # if ref_com

        lok = [rr[0] for rr in lrcom]
        ref_com = ds._generic_check._check_var(
            ref_com, 'ref_com',
            types=str,
            allowed=lok,
        )

        # -----------
        # safety check vs assumptions

        if isinstance(x0, str) and ref_com is coll.ddata[x0]['ref']:
            msg = (
                "For submesh=True, x0 '{x0}' should have no common ref!"
                f"\n\t- detected: {ref_com}"
            )
            raise NotImplementedError(msg)

    return ref_com


def _submesh_addtemp(
    coll=None,
    kd0=None,
    # interp resut
    dout_temp=None,
    # coordinates
    x0=None,
):

    # ----------------
    # add ref and data

    if isinstance(x0, str):

        assert not any([rr is None for rr in dout_temp[kd0[0]]['ref']])
        dr_add, lr_add = None, None

    else:

        # dref
        ii, dr_add, lref = 0, {}, []
        for jj, rr in enumerate(dout_temp[kd0[0]]['ref']):

            # ref name
            if rr is None:
                rr = f"r{len(coll.dref) + ii:03.0f}"
                lref.append(rr)
                ii += 1
            else:
                lref.append(None)

            # dr_add
            if rr not in coll.dref.keys():
                dr_add[rr] = {
                    'key': rr,
                    'size': dout_temp[kd0[0]]['data'].shape[jj],
                }

        # combine
        for kk in dout_temp.keys():
            dout_temp[kk]['ref'] = tuple([
                rr if rr is not None else dout_temp[kd0[0]]['ref'][jj]
                for jj, rr in enumerate(lref)
            ])

        lr_add = [vv['key'] for vv in dr_add.values()]

    # --------------
    # populate dd_add

    dd_add = dout_temp
    for ii, kk in enumerate(dd_add.keys()):
        dd_add[kk]['key'] = f'dddd{len(coll.ddata) + ii:03.0f}'
    ld_add = [vv['key'] for vv in dd_add.values()]

    # --------
    # add

    # ref
    if dr_add is not None:
        for vv in dr_add.values():
            coll.add_ref(**vv)

    # data
    for vv in dd_add.values():
        coll.add_data(**vv)

    return lr_add, ld_add


# ################################################################
# ################################################################
#                       Prepare bsplines
# ################################################################


def _prepare_bsplines(
    coll=None,
    ref_key=None,
):

    # ----------------
    # keys and types

    wbs = coll._which_bsplines
    keybs = ref_key
    keym = coll.dobj[wbs][keybs]['mesh']
    ref_key = coll.dobj[wbs][keybs]['apex']
    mtype = coll.dobj[coll._which_mesh][keym]['type']

    # -------------
    # azone

    # azone = ds._generic_check._check_var(
    # azone, 'azone',
    # types=bool,
    # default=True,
    # )

    # -----------
    # indbs

    return keybs, keym, mtype, ref_key


# #############################################################################
# #############################################################################
#                          core interpolation
# #############################################################################


def _interp(
    coll=None,
    keybs=None,
    keys=None,
    x0=None,
    x1=None,
    val_out=None,
    deriv=None,
    indbs_tf=None,
    crop=None,
    cropbs=None,
    crop_path=None,
    details=None,
    dout=None,
    mtype=None,
    ddata=None,
    daxis=None,
    sli_c=None,
    sli_x=None,
    sli_v=None,
    dsh_other=None,
    dref_com=None,
    lr_add=None,
    ld_add=None,
    nan0=None,
):

    # ----------
    # prepare

    # indokx
    indokx0 = np.isfinite(x0)
    if x1 is not None:
        indokx0 &= np.isfinite(x1)

    # loop on keys
    derr = {}
    wbs = coll._which_bsplines
    clas = coll.dobj[wbs][keybs]['class']
    for ii, k0 in enumerate(keys):

        try:

            if details is True:
                dout[k0]['data'] = clas.ev_details(
                    x0=x0,
                    x1=x1,
                    # options
                    val_out=val_out,
                    deriv=deriv,
                    # indices
                    indbs_tf=indbs_tf,
                    # rect-specific
                    crop=crop,
                    cropbs=cropbs,
                )

            else:

                # --------------------
                # prepare slicing func

                if mtype == 'tri':
                    shape_c = ddata[k0].shape

                    (
                        shape_v, axis_v, ind_c, ind_v, shape_o,
                    ) = _utils_bsplines._get_shapes_ind(
                        axis=daxis[k0],
                        shape_c=shape_c,
                        shape_x=x0.shape,
                    )

                    sli_c = _utils_bsplines._get_slice_cx(
                        axis=daxis[k0],
                        shape=shape_c,
                        ind_cv=ind_c,
                        reverse=True,
                    )

                    sli_v = _utils_bsplines._get_slice_cx(
                        axis=axis_v,
                        shape=shape_v,
                        ind_cv=ind_v,
                        reverse=True,
                    )
                else:
                    axis_v = None

                sli_o = _utils_bsplines._get_slice_out(
                    axis=daxis[k0],
                    shape_c=ddata[k0].shape,
                )

                # --------------------
                # actual interpolation

                dout[k0]['data'][...] = clas(
                    x0=x0,
                    x1=x1,
                    # coefs
                    coefs=ddata[k0],
                    axis=daxis[k0],
                    # options
                    val_out=val_out,
                    deriv=deriv,
                    # rect-specific
                    crop=crop,
                    cropbs=cropbs,
                    crop_path=crop_path,
                    # slicing
                    sli_c=sli_c,
                    sli_x=sli_x,
                    sli_v=sli_v,
                    sli_o=sli_o,
                    indokx0=indokx0,
                    shape_o=dsh_other[k0],
                    shape_v=dout[k0]['data'].shape,
                    dref_com=dref_com[k0],
                    axis_v=axis_v,
                )

                # nan0
                if nan0 is True:
                    dout[k0]['data'][dout[k0]['data'] == 0.] = np.nan

        except Exception as err:
            raise err
            derr[k0] = str(err)

    # ----------------------------
    # raise warning if any failure

    if len(derr) > 0:
        lstr = [f"\t- '{k0}': {v0}" for k0, v0 in derr.items()]
        msg = (
            "The following keys could not be interpolated:\n"
            + "\n".join(lstr)
        )
        warnings.warn(msg)

    # ----------------------------
    # clean-up temporary storage

    if lr_add is not None:
        for rr in lr_add:
            coll.remove_ref(rr, propagate=False)

    if ld_add is not None:
        for dd in ld_add:
            if dd in coll.ddata.keys():
                coll.remove_data(dd, propagate=False)

    # ----------
    # adjust ref

    if lr_add is not None:
        for k0, v0 in dout.items():
            dout[k0]['ref'] = tuple([
                None if rr in lr_add else rr
                for rr in v0['ref']
            ])

    if details:
        wbs = coll._which_bsplines
        for k0, v0 in dout.items():
            dout[k0]['ref'] = tuple(np.r_[
                v0['ref'], coll.dobj[wbs][keybs]['ref_bs'],
            ])

    return
