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
#                           Mesh2DRect - interp
# #############################################################################


def _interp2d_check_RZ(
    R=None,
    Z=None,
    grid=None,
):
    # R and Z provided
    if not isinstance(R, np.ndarray):
        try:
            R = np.atleast_1d(R).astype(float)
        except Exception as err:
            msg = "R must be convertible to np.arrays of floats"
            raise Exception(msg)
    if not isinstance(Z, np.ndarray):
        try:
            Z = np.atleast_1d(Z).astype(float)
        except Exception as err:
            msg = "Z must be convertible to np.arrays of floats"
            raise Exception(msg)

    # grid
    grid = ds._generic_check._check_var(
        grid, 'grid',
        default=R.shape != Z.shape,
        types=bool,
    )

    if grid is True and (R.ndim > 1 or Z.ndim > 1):
        msg = "If grid=True, R and Z must be 1d!"
        raise Exception(msg)
    elif grid is False and R.shape != Z.shape:
        msg = "If grid=False, R and Z must have the same shape!"
        raise Exception(msg)

    if grid is True:
        R = np.tile(R, Z.size)
        Z = np.repeat(Z, R.size)
    return R, Z, grid


def _interp2d_check(
    # ressources
    coll=None,
    # interpolation base, 1d or 2d
    key=None,
    # external coefs (optional)
    coefs=None,
    # interpolation points
    R=None,
    Z=None,
    radius=None,
    angle=None,
    grid=None,
    radius_vs_time=None,
    azone=None,
    # time: t or indt
    t=None,
    indt=None,
    indt_strict=None,
    # bsplines
    indbs=None,
    # parameters
    details=None,
    reshape=None,
    res=None,
    crop=None,
    nan0=None,
    val_out=None,
    imshow=None,
    return_params=None,
    store=None,
    inplace=None,
):

    # -------------
    # keys

    dk = {
        k0: v0['bsplines']
        for k0, v0 in coll.ddata.items()
        if v0.get('bsplines') not in ['', None]
        and 'crop' not in k0
    }
    dk.update({kk: kk for kk in coll.dobj['bsplines'].keys()})
    if key is None and len(dk) == 1:
        key = list(dk.keys())[0]
    if key not in dk.keys():
        msg = (
            "Arg key must the key to a data referenced on a bsplines set\n"
            f"\t- available: {list(dk.keys())}\n"
            f"\t- provided: {key}\n"
        )
        raise Exception(msg)

    keybs = dk[key]
    keym = coll.dobj['bsplines'][keybs]['mesh']
    mtype = coll.dobj[coll._which_mesh][keym]['type']

    # -------------
    # details

    details = ds._generic_check._check_var(
        details, 'details',
        types=bool,
        default=False,
    )

    # -------------
    # crop

    crop = ds._generic_check._check_var(
        crop, 'crop',
        types=bool,
        default=mtype == 'rect',
        allowed=[False, True] if mtype == 'rect' else [False],
    )

    # -------------
    # nan0

    nan0 = ds._generic_check._check_var(
        nan0, 'nan0',
        types=bool,
        default=True,
    )

    # -------------
    # val_out

    val_out = ds._generic_check._check_var(
        val_out, 'val_out',
        default=np.nan,
        allowed=[False, np.nan, 0.]
    )

    # -----------
    # time

    # refbs and hastime
    refbs = coll.dobj['bsplines'][keybs]['ref']

    if mtype == 'polar':
        radius2d = coll.dobj[coll._which_mesh][keym]['radius2d']
        # take out key if bsplines
        lk = [kk for kk in [key, radius2d] if kk != keybs]

        # 
        hastime, reft, keyt, t, dind = coll.get_time_common(
            keys=lk,
            t=t,
            indt=indt,
            ind_strict=indt_strict,
        )

        rad2d_hastime = radius2d in dind.keys()
        if key == keybs:
            kind = radius2d
        else:
            kind = key
            hastime = key in dind.keys()

        if hastime:

            indt = dind[kind].get('ind')

            # Special case: all times match
            if indt is not None:
                rtk = coll.get_time(key)[2]
                if indt.size == coll.dref[rtk]['size']:
                    if np.allclose(indt, np.arange(0, coll.dref[rtk]['size'])):
                        indt = None

            if indt is not None:
                indtu = np.unique(indt)
                indtr = np.array([indt == iu for iu in indtu])
            else:
                indtu, indtr = None, None
        else:
            indt, indtu, indtr = None, None, None

    elif key != keybs:
        # hastime, t, indit
        hastime, hasvect, reft, keyt, t, dind = coll.get_time(
            key=key,
            t=t,
            indt=indt,
            ind_strict=indt_strict,
        )
        if dind is None:
            indt, indtu, indtr = None, None, None
        else:
            indt, indtu, indtr = dind['ind'], dind['indu'], dind['indr']
    else:
        hastime, hasvect = False, False
        reft, keyt, indt, indtu, indtr = None, None, None, None, None

    # -----------
    # indbs

    if details is True:
        if mtype == 'rect':
            returnas = 'tuple-flat'
        elif mtype == 'tri':
            returnas = int
        elif mtype == 'polar':
            if len(coll.dobj['bsplines'][keybs]['shape']) == 2:
                returnas = 'array-flat'
            else:
                returnas = int

        # compute validated indbs array with appropriate form
        indbs_tf = coll.select_ind(
            key=keybs,
            returnas=returnas,
            ind=indbs,
        )
    else:
        indbs_tf = None

    # -----
    # store

    store = ds._generic_check._check_var(
        store, 'store',
        types=bool,
        default=False,
        allowed=[False] if details else [True, False],
    )

    # -----------
    # coordinates

    # (R, Z) vs (radius, angle)
    lc = [
        R is None and Z is None and radius is None,
        R is not None and Z is not None and store is False,
        (R is None and Z is None)
        and (radius is not None and mtype == 'polar') and store is False
    ]

    if not any(lc):
        msg = (
            "Please provide either:\n"
            "\t- R and Z: for any mesh type\n"
            "\t- radius (and angle): for polar mesh\n"
            "\t- None: for rect / tri mesh\n"
        )
        raise Exception(msg)

    # R, Z
    if lc[0]:

        if store is True:
            if imshow:
                msg = "Storing only works if imshow = False"
                raise Exception(msg)

        # no spec => sample mesh
        R, Z = coll.get_sample_mesh(
            key=keym,
            res=res,
            mode='abs',
            grid=True,
            R=R,
            Z=Z,
            imshow=imshow,
        )
        lc[1] = True

    if lc[1]:

        # check R, Z
        R, Z, grid = _interp2d_check_RZ(R=R, Z=Z, grid=grid)

        # special case if polar mesh => (radius, angle) from (R, Z)
        if mtype == 'polar':

            rad2d_indt = dind[radius2d].get('ind') if rad2d_hastime else None

            # compute radius2d at relevant times
            radius, _, _ = coll.interpolate_profile2d(
                # coordinates
                R=R,
                Z=Z,
                grid=False,
                # quantities
                key=radius2d,
                details=False,
                # time
                indt=rad2d_indt,
            )

            # compute angle2d at relevant times
            angle2d = coll.dobj[coll._which_mesh][keym]['angle2d']
            if angle2d is not None:
                angle, _, _ = coll.interpolate_profile2d(
                    # coordinates
                    R=R,
                    Z=Z,
                    grid=False,
                    # quantities
                    key=angle2d,
                    details=False,
                    # time
                    indt=rad2d_indt,
                )

            # simplify if not time-dependent
            if radius2d not in dind:
                assert radius.ndim == R.ndim
                if angle2d is not None:
                    assert angle.ndim == R.ndim

            # check consistency
            if rad2d_hastime != (radius.ndim == R.ndim + 1):
                msg = f"Inconsistency! {radius.shape}"
                raise Exception(msg)

            radius_vs_time = rad2d_hastime

        else:
            radius_vs_time = False

    else:
        radius_vs_time = ds._generic_check._check_var(
            radius_vs_time, 'radius_vs_time',
            types=bool,
            default=False,
        )

    # -------------
    # radius, angle

    if mtype == 'polar':

        # check same shape
        if not isinstance(radius, np.ndarray):
            radius = np.atleast_1d(radius)

        # angle vs angle2d
        angle2d = coll.dobj[coll._which_mesh][keym]['angle2d']
        if angle2d is not None and angle is None:
            msg = (
                f"Arg angle must be provided for bsplines {keybs}"
            )
            raise Exception(msg)

        # angle vs radius
        if angle is not None:
            if not isinstance(angle, np.ndarray):
                angle = np.atleast_1d(angle)

            if radius.shape != angle.shape:
                msg = (
                    "Args radius and angle must be np.ndarrays of same shape!\n"
                    f"\t- radius.shape: {radius.shape}\n"
                    f"\t- angle.shape: {angle.shape}\n"
                )
                raise Exception(msg)

    # -------------
    # coefs

    shapebs = coll.dobj['bsplines'][keybs]['shape']
    # 3 possible coefs shapes:
    #   - None (if details = True or key provided)
    #   - scalar or shapebs if details = False and key = keybs

    if details is True:
        coefs = None
        axis = None
        assert key == keybs, (key, keybs)

    else:
        if coefs is None:
            if key == keybs:
                if mtype == 'polar' and rad2d_hastime:
                    if rad2d_indt is None:
                        r2dnt = coll.dref[dind[radius2d]['ref']]['size']
                    else:
                        r2dnt = rad2d_indt.size
                    coefs = np.ones(tuple(np.r_[r2dnt, shapebs]), dtype=float)
                    axis = [ii + 1 for ii in range(len(shapebs))]
                else:
                    coefs = np.ones(shapebs, dtype=float)
                    axis = [ii for ii in range(len(shapebs))]
            else:
                coefs = coll.ddata[key]['data']
                axis = [
                    ii for ii in range(len(coll.ddata[key]['ref']))
                    if coll.ddata[key]['ref'][ii] in refbs
                ]

        elif key != keybs:
            msg = f"Arg coefs can only be provided if key = keybs!\n\t- key: {key}"
            raise Exception(msg)

        elif np.isscalar(coefs):
            coefs = np.full(shapebs, coefs)
            axis = [ii for ii in range(len(shapebs))]

        # consistency
        nshbs = len(shapebs)
        c0 = (
            coefs.shape[-nshbs:] == shapebs
            and (
                (hastime and coefs.ndim == nshbs + 1)           # pre-shaped
                or (not hastime and coefs.ndim == nshbs + 1)    # pre-shaped
                or (not hastime and coefs.ndim == nshbs)        # [None, ...]
            )
            and len(axis) == len(refbs)
        )
        if not c0:
            msg = (
                f"Inconsistency of '{key}' shape:\n"
                f"\t- shape: {coefs.shape}\n"
                f"\t- shapebs: {shapebs}\n"
                f"\t- hastime: {hastime}\n"
                f"\t- radius_vs_time: {radius_vs_time}\n"
                f"\t- refbs: {refbs}\n"
                f"\t- axis: {axis}"
            )
            raise Exception(msg)

        # Make sure coefs is time dependent
        if hastime:
            if indt is not None and (mtype == 'polar' or indtu is None):
                if coefs.shape[0] != indt.size:
                    # in case coefs is already provided with indt
                    coefs = coefs[indt, ...]

            if radius_vs_time and coefs.shape[0] != radius.shape[0]:
                msg = (
                    "Inconstistent coefs vs radius!\n"
                    f"\t- coefs.shape = {coefs.shape}\n"
                    f"\t- radius.shape = {radius.shape}\n"
                )
                raise Exception(msg)
        else:
            if radius_vs_time is True:
                sh = tuple([radius.shape[0]] + [1]*len(shapebs))
                coefs = np.tile(coefs, sh)
                axis = [aa + 1 for aa in axis]
            elif coefs.ndim == nshbs:
                coefs = coefs[None, ...]
                axis = [aa + 1 for aa in axis]

    # -------------
    # azone

    azone = ds._generic_check._check_var(
        azone, 'azone',
        types=bool,
        default=True,
    )

    # -------------
    # return_params

    return_params = ds._generic_check._check_var(
        return_params, 'return_params',
        types=bool,
        default=False,
    )

    # -------
    # inplace

    inplace = ds._generic_check._check_var(
        inplace, 'inplace',
        types=bool,
        default=store,
    )

    return (
        key, keybs,
        R, Z,
        radius, angle,
        coefs,
        axis,
        hastime,
        reft, keyt,
        shapebs,
        radius_vs_time,
        azone,
        indbs, indbs_tf,
        t, indt, indtu, indtr,
        details, crop,
        nan0, val_out,
        return_params,
        store, inplace,
    )


def interp2d(
    # ressources
    coll=None,
    # interpolation base, 1d or 2d
    key=None,
    # external coefs (instead of key, optional)
    coefs=None,
    # interpolation points
    R=None,
    Z=None,
    radius=None,
    angle=None,
    grid=None,
    radius_vs_time=None,        # if radius is provided, in case radius vs time 
    azone=None,                 # if angle2d is interpolated, exclusion zone
    # time: t or indt
    t=None,
    indt=None,
    indt_strict=None,
    # bsplines
    indbs=None,
    # parameters
    details=None,
    reshape=None,
    res=None,
    crop=None,
    nan0=None,
    val_out=None,
    imshow=None,
    return_params=None,
    store=None,
    inplace=None,
):

    # ---------------
    # check inputs

    (
        key, keybs,
        R, Z,
        radius, angle,
        coefs,
        axis,
        hastime,
        reft, keyt,
        shapebs,
        radius_vs_time,
        azone,
        indbs, indbs_tf,
        t, indt, indtu, indtr,
        details, crop,
        nan0, val_out,
        return_params,
        store, inplace,
    ) = _interp2d_check(
        # ressources
        coll=coll,
        # interpolation base, 1d or 2d
        key=key,
        # external coefs (optional)
        coefs=coefs,
        # interpolation points
        R=R,
        Z=Z,
        radius=radius,
        angle=angle,
        grid=grid,
        radius_vs_time=radius_vs_time,
        azone=azone,
        # time: t or indt
        t=t,
        indt=indt,
        indt_strict=indt_strict,
        # bsplines
        indbs=indbs,
        # parameters
        details=details,
        res=res,
        crop=crop,
        nan0=nan0,
        val_out=val_out,
        imshow=imshow,
        return_params=return_params,
        store=store,
        inplace=inplace,
    )
    keym = coll.dobj['bsplines'][keybs]['mesh']
    meshtype = coll.dobj['mesh'][keym]['type']

    # ----------------------
    # which function to call

    if details is True:
        fname = 'func_details'
    elif details is False:
        fname = 'func_sum'
    else:
        raise Exception("Unknown details!")

    # ---------------
    # cropping ?

    cropbs = coll.dobj['bsplines'][keybs]['crop']
    if cropbs not in [None, False]:
        cropbs = coll.ddata[cropbs]['data']
        nbs = cropbs.sum()
    else:
        nbs = np.prod(coll.dobj['bsplines'][keybs]['shape'])

    if indbs_tf is not None:
        if isinstance(indbs_tf, tuple):
            nbs = indbs_tf[0].size
        else:
            nbs = indbs_tf.size

    # -----------
    # Interpolate

    if meshtype in ['rect', 'tri']:

        # manage time
        if indtu is not None:
            val0 = np.full(tuple(np.r_[indt.size, R.shape]), np.nan)
            coefs = coefs[indtu, ...]

        val = coll.dobj['bsplines'][keybs][fname](
            R=R,
            Z=Z,
            coefs=coefs,
            axis=axis,
            crop=crop,
            cropbs=cropbs,
            indbs_tf=indbs_tf,
            val_out=val_out,
        )

        # manage time
        if indtu is not None:
            for ii, iu in enumerate(indtu):
                val0[indtr[ii], ...] = val[ii, ...]
            val = val0

        shape_pts = R.shape

    elif meshtype == 'polar':

        val = coll.dobj['bsplines'][keybs][fname](
            radius=radius,
            angle=angle,
            coefs=coefs,
            axis=axis,
            indbs_tf=indbs_tf,
            radius_vs_time=radius_vs_time,
            val_out=val_out,
        )

        shape_pts = radius.shape

    # ---------------
    # post-treatment for angle2d only (discontinuity)

    if azone is True:
        lkm0 = [
            k0 for k0, v0 in coll.dobj[coll._which_mesh].items()
            if key == v0.get('angle2d')
        ]
        if len(lkm0) > 0:
            ind = angle2d_inzone(
                coll=coll,
                keym0=lkm0[0],
                keya2d=key,
                R=R,
                Z=Z,
                t=t,
                indt=indt,
            )
            assert val.shape == ind.shape
            val[ind] = np.pi

    # ---------------
    # post-treatment

    if nan0 is True:
        val[val == 0] = np.nan

    if not hastime and not radius_vs_time:
        c0 = (
            (
                details is False
                and val.shape == tuple(np.r_[1, shape_pts])
            )
            or (
                details is True
                and val.shape == tuple(np.r_[shape_pts, nbs])
            )
        )
        if not c0:
            import pdb; pdb.set_trace()     # DB
            pass
        if details is False:
            val = val[0, ...]
            reft = None

    # ------
    # store

    # ref
    ct = (
        (hastime or radius_vs_time)
        and (
            reft in [None, False]
            or (
                reft not in [None, False]
                and indt is not None
                and not (
                    indt.size == coll.dref[reft]['size']
                    or np.allclose(indt, np.arange(0, coll.dref[reft]['size']))
                )
            )
        )
    )
    if ct:
        reft = f'{key}-nt'

    # store
    if store is True:
        Ru = np.unique(R)
        Zu = np.unique(Z)
        nR, nZ = Ru.size, Zu.size

        knR, knZ = f'{key}-nR', f'{key}-nZ'
        kR, kZ = f'{key}-R', f'{key}-Z'

        if inplace is True:
            coll2 = coll
        else:
            coll2 = ds.DataStock()

        # add ref nR, nZ
        coll2.add_ref(key=knR, size=nR)
        coll2.add_ref(key=knZ, size=nZ)

        # add data Ru, Zu
        coll2.add_data(
            key=kR,
            data=Ru,
            ref=knR,
            dim='distance',
            name='R',
            units='m',
        )
        coll2.add_data(
            key=kZ,
            data=Zu,
            ref=knZ,
            dim='distance',
            name='Z',
            units='m',
        )

        # ref
        if hastime or radius_vs_time:
            if ct:
                coll2.add_ref(key=reft, size=t.size)
                coll2.add_data(
                    key=f'{key}-t',
                    data=t,
                    ref=reft,
                    dim='time',
                    units='s',
                )
            else:
                if reft not in coll2.dref.keys():
                    coll2.add_ref(key=reft, size=coll.dref[reft]['size'])
                if keyt is not None and keyt not in coll2.ddata.keys():
                    coll2.add_data(
                        key=keyt,
                        data=coll.ddata[keyt]['data'],
                        dim='time',
                    )

            ref = (reft, knR, knZ)
        else:
            ref = (knR, knZ)

        coll2.add_data(
            data=val,
            key=f'{key}-map',
            ref=ref,
            dim=coll.ddata[key]['dim'],
            quant=coll.ddata[key]['quant'],
            name=coll.ddata[key]['name'],
            units=coll.ddata[key]['units'],
        )

    else:

        ref = []
        c0 = (
            reft not in [None, False]
            and (hastime or radius_vs_time)
            and not (
                meshtype == 'polar'
                and R is None
                and coefs is None
                and radius_vs_time is False
            )
        )
        if c0:
            ref.append(reft)

        if meshtype in ['rect', 'tri']:
            for ii in range(R.ndim):
                ref.append(None)
            if grid is True:
                for ii in range(Z.ndim):
                    ref.append(None)
        else:
            for ii in range(radius.ndim):
                if radius_vs_time and ii == 0:
                    continue
                ref.append(None)
            if grid is True and angle is not None:
                for ii in range(angle.ndim):
                    ref.append(None)

        if details is True:
            refbs = coll.dobj['bsplines'][keybs]['ref-bs'][0]
            if crop is True:
                refbs = f"{refbs}-crop"
            ref.append(refbs)
        ref = tuple(ref)

        if ref[0] == 'emiss-nt':
            import pdb; pdb.set_trace()     # DB

        if len(ref) != val.ndim:
            msg = (
                "Mismatching ref vs val.shape:\n"
                f"\t- key = {key}\n"
                f"\t- keybs = {keybs}\n"
                f"\t- val.shape = {val.shape}\n"
                f"\t- ref = {ref}\n"
                f"\t- reft = {reft}\n"
                f"\t- hastime = {hastime}\n"
                f"\t- radius_vs_time = {radius_vs_time}\n"
                f"\t- details = {details}\n"
                f"\t- indbs_tf = {indbs_tf}\n"
                f"\t- key = {key}\n"
                f"\t- meshtype = {meshtype}\n"
                f"\t- grid = {grid}\n"
            )
            if coefs is not None:
                msg += f"\t- coefs.shape = {coefs.shape}\n"
            if R is not None:
                msg += (
                    f"\t- R.shape = {R.shape}\n"
                    f"\t- Z.shape = {Z.shape}\n"
                )
            if meshtype == 'polar':
                msg += f"\t- radius.shape = {radius.shape}\n"
                if angle is not None:
                    msg += f"\t- angle.shape = {angle.shape}\n"
            raise Exception(msg)

    # ------
    # return

    if store and inplace is False:
        return coll2
    else:
        if return_params is True:
            return val, t, ref, dparams
        else:
            return val, t, ref


# #############################################################################
# #############################################################################
#                   2d points to 1d quantity interpolation
# #############################################################################


def interp2dto1d(
    coll=None,
    key=None,
    R=None,
    Z=None,
    indbs=None,
    indt=None,
    grid=None,
    details=None,
    reshape=None,
    res=None,
    crop=None,
    nan0=None,
    imshow=None,
    return_params=None,
):

    # ---------------
    # check inputs

    # TBD
    pass

    # ---------------
    # post-treatment

    if nan0 is True:
        val[val == 0] = np.nan

    # ------
    # return

    if return_params is True:
        return val, dparams
    else:
        return val


# #############################################################################
# #############################################################################
#                           Mesh2DRect - operators
# #############################################################################


def get_bsplines_operator(
    coll,
    key=None,
    operator=None,
    geometry=None,
    crop=None,
    store=None,
    returnas=None,
    # specific to deg = 0
    centered=None,
    # to return gradR, gradZ, for D1N2 deg 0, for tomotok
    returnas_element=None,
):

    # check inputs
    lk = list(coll.dobj.get('bsplines', {}).keys())
    key = ds._generic_check._check_var(
        key, 'key',
        types=str,
        allowed=lk,
    )

    store = ds._generic_check._check_var(
        store, 'store',
        default=True,
        types=bool,
    )

    returnas = ds._generic_check._check_var(
        returnas, 'returnas',
        default=store is False,
        types=bool,
    )

    crop = ds._generic_check._check_var(
        crop, 'crop',
        default=True,
        types=bool,
    )

    # cropbs
    cropbs = coll.dobj['bsplines'][key]['crop']
    keycropped = coll.dobj['bsplines'][key]['ref-bs'][0]
    if cropbs not in [None, False] and crop is True:
        cropbs_flat = coll.ddata[cropbs]['data'].ravel(order='F')
        if coll.dobj['bsplines'][key]['deg'] == 0:
            cropbs = coll.ddata[cropbs]['data']
        keycropped = f"{keycropped}-crop"
    else:
        cropbs = False
        cropbs_flat = False

    # compute and return
    (
        opmat, operator, geometry, dim,
    ) = coll.dobj['bsplines'][key]['class'].get_operator(
        operator=operator,
        geometry=geometry,
        cropbs_flat=cropbs_flat,
        # specific to deg=0
        cropbs=cropbs,
        centered=centered,
        # to return gradR, gradZ, for D1N2 deg 0, for tomotok
        returnas_element=returnas_element,
    )

    # cropping
    if operator == 'D1':
        ref = (keycropped, keycropped)
    elif operator == 'D0N1':
        ref = (keycropped,)
    elif 'N2' in operator:
        ref = (keycropped, keycropped)

    return opmat, operator, geometry, dim, ref, crop, store, returnas, key


# #################################################################
# #################################################################
#               Contour computation
# #################################################################


def _get_contours(
    RR=None,
    ZZ=None,
    val=None,
    levels=None,
    largest=None,
    uniform=None,
):
    """ Return R, Z coordinates of contours (time-dependent)

    For contourpy algorithm, the dimensions shoud be (ny, nx), from meshgrid

    RR = (nz, nr)
    ZZ = (nz, nr)
    val = (nt, nz, nr)
    levels = (nlevels,)

    cR = (nt, nlevels, nmax) array of R coordinates
    cZ = (nt, nlevels, nmax) array of Z coordinates

    The contour coordinates are uniformzied to always have the same nb of pts

    """

    # -------------
    # check inputs

    if largest is None:
        largest = False

    if uniform is None:
        uniform = True

    # val.shape = (nt, nR, nZ)
    lc = [
        val.shape == RR.shape,
        val.ndim == RR.ndim + 1 and val.shape[1:] == RR.shape,
    ]
    if lc[0]:
        val = val[None, ...]
    elif lc[1]:
        pass
    else:
        msg = "Incompatible val.shape!"
        raise Exception(msg)

    nt, nR, nZ = val.shape

    # ------------------------
    # Compute list of contours

    # compute contours at rknots
    # see https://github.com/matplotlib/matplotlib/blob/main/src/_contour.h

    contR = [[] for ii in range(nt)]
    contZ = [[] for ii in range(nt)]
    for ii in range(nt):

        # define map
        contgen = contour_generator(
            x=RR,
            y=ZZ,
            z=val[ii, ...],
            name='serial',
            corner_mask=None,
            line_type='Separate',
            fill_type=None,
            chunk_size=None,
            chunk_count=None,
            total_chunk_count=None,
            quad_as_tri=True,       # for sub-mesh precision
            # z_interp=<ZInterp.Linear: 1>,
            thread_count=0,
        )

        for jj in range(len(levels)):

            # compute concatenated contour
            no_cont = False
            cj = contgen.lines(levels[jj])

            c0 = (
                isinstance(cj, list)
                and all([
                    isinstance(cjj, np.ndarray)
                    and cjj.ndim == 2
                    and cjj.shape[1] == 2
                    for cjj in cj
                ])
            )
            if not c0:
                msg = f"Wrong output from contourpy!\n{cj}"
                raise Exception(msg)

            if len(cj) > 0:
                cj = [
                    cc[np.all(np.isfinite(cc), axis=1), :]
                    for cc in cj
                    if np.sum(np.all(np.isfinite(cc), axis=1)) >= 3
                ]

                if len(cj) == 0:
                    no_cont = True
                elif len(cj) == 1:
                    cj = cj[0]
                elif len(cj) > 1:
                    if largest:
                        nj = [
                            0.5*np.abs(np.sum(
                                (cc[1:, 0] + cc[:-1, 0])
                                *(cc[1:, 1] - cc[:-1, 1])
                            ))
                            for cc in cj
                        ]
                        cj = cj[np.argmax(nj)]
                    else:
                        ij = np.cumsum([cc.shape[0] for cc in cj])
                        cj = np.concatenate(cj, axis=0)
                        cj = np.insert(cj, ij, np.nan, axis=0)

                elif np.sum(np.all(~np.isnan(cc), axis=1)) < 3:
                    no_cont = True
            else:
                no_cont = True

            if no_cont is True:
                cj = np.full((3, 2), np.nan)

            contR[ii].append(cj[:, 0])
            contZ[ii].append(cj[:, 1])

    # ------------------------------------------------
    # Interpolate / concatenate to uniformize as array

    if uniform:
        ln = [[pp.size for pp in cc] for cc in contR]
        nmax = np.max(ln)
        cR = np.full((nt, len(levels), nmax), np.nan)
        cZ = np.full((nt, len(levels), nmax), np.nan)

        for ii in range(nt):
            for jj in range(len(levels)):
                cR[ii, jj, :] = np.interp(
                    np.linspace(0, ln[ii][jj], nmax),
                    np.arange(0, ln[ii][jj]),
                    contR[ii][jj],
                )
                cZ[ii, jj, :] = np.interp(
                    np.linspace(0, ln[ii][jj], nmax),
                    np.arange(0, ln[ii][jj]),
                    contZ[ii][jj],
                )

        return cR, cZ
    else:
        return contR, contZ


# #############################################################################
# #############################################################################
#                   Polygon simplification
# #############################################################################


def _simplify_polygon(pR=None, pZ=None, res=None):
    """ Use convex hull with a constraint on the maximum discrepancy """

    # ----------
    # preliminary 1: check there is non redundant point

    dp = np.sqrt((pR[1:] - pR[:-1])**2 + (pZ[1:] - pZ[:-1])**2)
    ind = (dp > 1.e-6).nonzero()[0]
    pR = pR[ind]
    pZ = pZ[ind]

    # check new poly is closed
    if (pR[0] != pR[-1]) or (pZ[0] != pZ[-1]):
        pR = np.append(pR, pR[0])
        pZ = np.append(pZ, pZ[0])

    # check it is counter-clockwise
    clock = np.nansum((pR[1:] - pR[:-1]) * (pZ[1:] + pZ[:-1]))
    if clock > 0:
        pR = pR[::-1]
        pZ = pZ[::-1]

    # threshold = diagonal of resolution + 10%
    thresh = res * np.sqrt(2) * 1.1

    # ----------
    # preliminary 2: get convex hull and copy

    poly = np.array([pR, pZ]).T
    iconv = ConvexHull(poly, incremental=False).vertices

    # close convex hull to iterate on edges
    pR_conv = np.append(pR[iconv], pR[iconv[0]])
    pZ_conv = np.append(pZ[iconv], pZ[iconv[0]])

    # copy to create new polygon that will serve as buffer
    pR_bis, pZ_bis = np.copy(pR), np.copy(pZ)

    # -------------------------
    # loop on convex hull edges

    for ii in range(pR_conv.size - 1):

        pR1, pR2 = pR_conv[ii], pR_conv[ii+1]
        pZ1, pZ2 = pZ_conv[ii], pZ_conv[ii+1]
        i0 = np.argmin(np.hypot(pR_bis - pR1, pZ_bis - pZ1))

        # make sure it starts from p1
        pR_bis = np.append(pR_bis[i0:], pR_bis[:i0])
        pZ_bis = np.append(pZ_bis[i0:], pZ_bis[:i0])

        # get indices of closest points to p1, p2
        i1 = np.argmin(np.hypot(pR_bis - pR1, pZ_bis - pZ1))
        i2 = np.argmin(np.hypot(pR_bis - pR2, pZ_bis - pZ2))

        # get corresponding indices of poly points to be included
        if i2 == i1 + 1:
            itemp = [i1, i2]

        else:
            # several points in-between
            # => check they are all within distance before exclusing them

            # get unit vector of segment
            norm12 = np.hypot(pR2 - pR1, pZ2 - pZ1)
            u12R = (pR2 - pR1) / norm12
            u12Z = (pZ2 - pZ1) / norm12

            # get points standing between p1 nd p2
            lpR = pR_bis[i1 + 1:i2]
            lpZ = pZ_bis[i1 + 1:i2]

            # indices of points standing too far from edge (use cross-product)
            iout = np.abs(u12R*(lpZ - pZ1) - u12Z*(lpR - pR1)) > thresh

            # if any pts too far => include all pts
            if np.any(iout):
                itemp = np.arange(i1, i2 + 1)
            else:
                itemp = [i1, i2]

        # build pts_in
        pR_in = pR_bis[itemp]
        pZ_in = pZ_bis[itemp]

        # concatenate to add to new polygon
        pR_bis = np.append(pR_in, pR_bis[i2 + 1:])
        pZ_bis = np.append(pZ_in, pZ_bis[i2 + 1:])

    # check new poly is closed
    if (pR_bis[0] != pR_bis[-1]) or (pZ_bis[0] != pZ_bis[-1]):
        pR_bis = np.append(pR_bis, pR_bis[0])
        pZ_bis = np.append(pZ_bis, pZ_bis[0])

    return pR_bis, pZ_bis


# #############################################################################
# #############################################################################
#                   radius2d special points handling
# #############################################################################


def radius2d_special_points(
    coll=None,
    key=None,
    keym0=None,
    res=None,
):

    keybs = coll.ddata[key]['bsplines']
    keym = coll.dobj['bsplines'][keybs]['mesh']
    mtype = coll.dobj[coll._which_mesh][keym]['type']
    assert mtype in ['rect', 'tri']

    # get map sampling
    RR, ZZ = coll.get_sample_mesh(
        key=keym,
        res=res,
        grid=True,
    )

    # get map
    val, t, _ = coll.interpolate_profile2d(
        key=key,
        R=RR,
        Z=ZZ,
        grid=False,
        imshow=True,        # for contour
    )

    # get min max values
    rmin = np.nanmin(val)
    rmax = np.nanmax(val)

    # get contour of 0
    cR, cZ = _get_contours(
        RR=RR,
        ZZ=ZZ,
        val=val,
        levels=[rmin + 0.05*(rmax-rmin)],
    )

    # dref
    ref_O = f'{keym0}-pts-O-n'
    dref = {
        ref_O: {'size': 1},
    }

    # get barycenter 
    if val.ndim == 3:
        assert cR.shape[1] == 1
        ax_R = np.nanmean(cR[:, 0, :], axis=-1)[:, None]
        ax_Z = np.nanmean(cZ[:, 0, :], axis=-1)[:, None]
        reft = coll.ddata[key]['ref'][0]
        ref = (reft, ref_O)
    else:
        ax_R = np.r_[np.nanmean(cR)]
        ax_Z = np.r_[np.nanmean(cZ)]
        ref = (ref_O,)

    kR = f'{keym0}-pts-O-R'
    kZ = f'{keym0}-pts-O-Z'
    ddata = {
        kR: {
            'ref': ref,
            'data': ax_R,
            'dim': 'distance',
            'quant': 'R',
            'name': 'O-points_R',
            'units': 'm',
        },
        kZ: {
            'ref': ref,
            'data': ax_Z,
            'dim': 'distance',
            'quant': 'Z',
            'name': 'O-points_Z',
            'units': 'm',
        },
    }

    return dref, ddata, kR, kZ


# #############################################################################
# #############################################################################
#                   angle2d discontinuity handling
# #############################################################################


def angle2d_zone(
    coll=None,
    key=None,
    keyrad2d=None,
    key_ptsO=None,
    res=None,
    keym0=None,
):

    keybs = coll.ddata[key]['bsplines']
    keym = coll.dobj['bsplines'][keybs]['mesh']
    mtype = coll.dobj[coll._which_mesh][keym]['type']
    assert mtype in ['rect', 'tri']

    # --------------
    # prepare

    hastime, hasvect, reft, keyt = coll.get_time(key=key)[:4]
    if hastime:
        nt = coll.dref[reft]['size']
    else:
        msg = (
            "Non time-dependent angle2d not implemented yet\n"
            "=> ping @Didou09 on Github to open an issue"
        )
        raise NotImplementedError(msg)

    if res is None:
        res = _get_sample_mesh_res(
            coll=coll,
            keym=keym,
            mtype=mtype,
        )

    # get map sampling
    RR, ZZ = coll.get_sample_mesh(
        key=keym,
        res=res/2.,
        grid=True,
        imshow=True,    # for contour
    )

    # get map
    val, t, _ = coll.interpolate_profile2d(
        key=key,
        R=RR,
        Z=ZZ,
        grid=False,
        azone=False,
    )
    val[np.isnan(val)] = 0.
    amin = np.nanmin(val)
    amax = np.nanmax(val)

    # get contours of absolute value
    cRmin, cZmin = _get_contours(
        RR=RR,
        ZZ=ZZ,
        val=val,
        levels=[amin + 0.10*(amax - amin)],
        largest=True,
        uniform=True,
    )
    cRmax, cZmax = _get_contours(
        RR=RR,
        ZZ=ZZ,
        val=val,
        levels=[amax - 0.10*(amax - amin)],
        largest=True,
        uniform=True,
    )

    cRmin, cZmin = cRmin[:, 0, :], cZmin[:, 0, :]
    cRmax, cZmax = cRmax[:, 0, :], cZmax[:, 0, :]

    rmin = np.full(cRmin.shape, np.nan)
    rmax = np.full(cRmax.shape, np.nan)

    # get points inside contour 
    for ii in range(nt):
        rmin[ii, :], _, _ = coll.interpolate_profile2d(
            key=keyrad2d,
            R=cRmin[ii, :],
            Z=cZmin[ii, :],
            grid=False,
            indt=ii,
        )
        rmax[ii, :], _, _ = coll.interpolate_profile2d(
            key=keyrad2d,
            R=cRmax[ii, :],
            Z=cZmax[ii, :],
            grid=False,
            indt=ii,
        )

    # get magnetic axis
    kR, kZ = key_ptsO
    axR = coll.ddata[kR]['data']
    axZ = coll.ddata[kZ]['data']
    assert coll.ddata[kR]['ref'][0] == coll.ddata[key]['ref'][0]

    start_min = np.nanargmin(rmin, axis=-1)
    start_max = np.nanargmin(rmax, axis=-1)

    # re-order from start_min, start_max
    lpR, lpZ = [], []
    for ii in range(rmin.shape[0]):
        imin = np.r_[
            np.arange(start_min[ii], rmin.shape[1]),
            np.arange(0, start_min[ii]),
        ]

        cRmin[ii] = cRmin[ii, imin]
        cZmin[ii] = cZmin[ii, imin]
        rmin[ii] = rmin[ii, imin]
        # check it is counter-clockwise
        clock = np.nansum(
            (cRmin[ii, 1:] - cRmin[ii, :-1])
            *(cZmin[ii, 1:] + cZmin[ii, :-1])
        )
        if clock > 0:
            cRmin[ii, :] = cRmin[ii, ::-1]
            cZmin[ii, :] = cZmin[ii, ::-1]
            rmin[ii, :] = rmin[ii, ::-1]

        imax = np.r_[
            np.arange(start_max[ii], rmax.shape[1]),
            np.arange(0, start_max[ii])
        ]
        cRmax[ii] = cRmax[ii, imax]
        cZmax[ii] = cZmax[ii, imax]
        rmax[ii] = rmax[ii, imax]
        # check it is clockwise
        clock = np.nansum(
            (cRmax[ii, 1:] - cRmax[ii, :-1])
            *(cZmax[ii, 1:] + cZmax[ii, :-1])
        )
        if clock < 0:
            cRmax[ii, :] = cRmax[ii, ::-1]
            cZmax[ii, :] = cZmax[ii, ::-1]
            rmax[ii, :] = rmax[ii, ::-1]

        # i0
        dr = np.diff(rmin[ii, :])
        i0 = (np.isnan(dr) | (dr < 0)).nonzero()[0][0]
        # rmin[ii, i0-1:] = np.nan
        dr = np.diff(rmax[ii, :])
        i1 = (np.isnan(dr) | (dr < 0)).nonzero()[0][0]
        # rmax[ii, i1-1:] = np.nan

        # polygon
        pR = np.r_[axR[ii], cRmin[ii, :i0-1], cRmax[ii, :i1-1][::-1]]
        pZ = np.r_[axZ[ii], cZmin[ii, :i0-1], cZmax[ii, :i1-1][::-1]]

        pR, pZ = _simplify_polygon(pR=pR, pZ=pZ, res=res)

        lpR.append(pR)
        lpZ.append(pZ)

    # Ajust sizes
    nb = np.array([pR.size for pR in lpR])

    # 
    nmax = np.max(nb)
    pR = np.full((nt, nmax), np.nan)
    pZ = np.full((nt, nmax), np.nan)

    for ii in range(nt):
        pR[ii, :] = np.interp(
            np.linspace(0, nb[ii], nmax),
            np.arange(0, nb[ii]),
            lpR[ii],
        )
        pZ[ii, :] = np.interp(
            np.linspace(0, nb[ii], nmax),
            np.arange(0, nb[ii]),
            lpZ[ii],
        )

    # ----------------
    # prepare output dict

    # ref
    kref = f'{keym0}-azone-npt'
    dref = {
        kref: {'size': nmax}
    }

    # data
    kR = f'{keym0}-azone-R'
    kZ = f'{keym0}-azone-Z'
    ddata = {
        kR: {
            'data': pR,
            'ref': (reft, kref),
            'units': 'm',
            'dim': 'distance',
            'quant': 'R',
            'name': None,
        },
        kZ: {
            'data': pZ,
            'ref': (reft, kref),
            'units': 'm',
            'dim': 'distance',
            'quant': 'R',
            'name': None,
        },
    }

    return dref, ddata, kR, kZ


def angle2d_inzone(
    coll=None,
    keym0=None,
    keya2d=None,
    R=None,
    Z=None,
    t=None,
    indt=None,
):


    # ------------
    # prepare points

    if R.ndim == 1:
        shape0 = None
        pts = np.array([R, Z]).T
    else:
        shape0 = R.shape
        pts = np.array([R.ravel(), Z.ravel()]).T

    # ------------
    # prepare path

    kazR, kazZ = coll.dobj[coll._which_mesh][keym0]['azone']
    pR = coll.ddata[kazR]['data']
    pZ = coll.ddata[kazZ]['data']

    hastime, hasvect, reft, keyt, tnew, dind = coll.get_time(
        key=kazR,
        t=t,
        indt=indt,
    )

    # ------------
    # test points

    if hastime:
        if dind is None:
            nt = coll.dref[reft]['size']
            ind = np.zeros((nt, R.size), dtype=bool)
            for ii in range(nt):
                path = Path(np.array([pR[ii, :], pZ[ii, :]]).T)
                ind[ii, :] = path.contains_points(pts)
        else:
            import pdb; pdb.set_trace()     # DB
            raise NotImplementedError()
            # TBC / TBF
            nt = None
            ind = np.zeros((nt, R.size), dtype=bool)
            for ii in range(nt):
                path = Path(np.array([pR[ii, :], pZ[ii, :]]).T)
                ind[ii, :] = path.contains_points(pts)

    else:
        path = Path(np.array([pR, pZ]).T)
        ind = path.contains_points(pts)

    # -------------------------
    # fromat output and return

    if shape0 is not None:
        if hastime:
            ind = ind.reshape(tuple(np.r_[nt, shape0]))
        else:
            ind = ind.reshape(shape0)

    return ind