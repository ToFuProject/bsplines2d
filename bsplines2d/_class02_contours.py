# -*- coding: utf-8 -*-


# Built-in
import itertools as itt
import warnings


# Common
import numpy as np
from contourpy import contour_generator
import datastock as ds


# specific


# #################################################################
# #################################################################
#               Main
# #################################################################


def get_contours(
    coll=None,
    key=None,
    levels=None,
    res=None,
    npts=None,
    largest=None,
    ref_com=None,
    # return vs store
    returnas=None,
    return_dref=None,
    store=None,
    key_npts=None,
    key_lvls=None,
    key_cont0=None,
    key_cont1=None,
):

    # ----------
    # check

    (
        key, keybs, keym0, levels,
        store, returnas, return_dref,
        key_npts, key_lvls, key_cont0, key_cont1,
    ) = _check(
        coll=coll,
        key=key,
        levels=levels,
        # return vs store
        returnas=returnas,
        return_dref=return_dref,
        store=store,
        key_npts=key_npts,
        key_lvls=key_lvls,
        key_cont0=key_cont0,
        key_cont1=key_cont1,
    )

    # ---------------------
    # prepare (interpolate)

    dsamp = coll.get_sample_mesh(
        key=keym0,
        res=res,
        grid=True,
        store=False,
    )

    # temporary addition of sample
    if ref_com is not None:

        # ref
        lkr_temp = []
        for ii, ss in enumerate(dsamp['x0']['data'].shape):
            kri = f'{key}_refsamptemp{ii}'
            lkr_temp.append(kri)
            coll.add_ref(key=kri, size=ss)

        # data
        lkd_temp = [f'{key}_samptemp0', f'{key}_samptemp1']
        coll.add_data(key=lkd_temp[0], data=dsamp['x0']['data'], ref=lkr_temp)
        coll.add_data(key=lkd_temp[1], data=dsamp['x1']['data'], ref=lkr_temp)
        x0, x1 = lkd_temp
    else:
        x0 = dsamp['x0']['data']
        x1 = dsamp['x1']['data']

    # interpolate
    dinterp = coll.interpolate(
        keys=key,
        ref_key=keybs,
        x0=x0,
        x1=x1,
        grid=False,
        submesh=True,
        ref_com=ref_com,
        store=False,
    )[key]

    # remove temporary data
    if ref_com is not None:
        coll.remove_data(x0, propagate=False)
        coll.remove_data(x1, propagate=False)
        for rr in lkr_temp:
            coll.remove_ref(rr, propagate=False)

        axis = [
            dinterp['ref'].index(lkr_temp[0]),
            dinterp['ref'].index(lkr_temp[1]),
        ]

    else:
        # axis
        axis = [
            ii for ii, rr in enumerate(dinterp['ref'])
            if rr is None
        ]

    # ----------------
    # compute contours

    cont0, cont1 = _get_contours(
        key=key,
        xx0=dsamp['x0']['data'],
        xx1=dsamp['x1']['data'],
        val=dinterp['data'],
        axis=axis,
        levels=levels,
        npts=npts,
        largest=largest,
    )

    # ----------------
    # format output

    dout, dref = _format(
        coll=coll,
        key=key,
        keybs=keybs,
        dinterp=dinterp,
        cont0=cont0,
        cont1=cont1,
        axis=axis,
        keym0=keym0,
        key_npts=key_npts,
        key_lvls=key_lvls,
        key_cont0=key_cont0,
        key_cont1=key_cont1,
    )

    # ----------------
    # store and return

    if store is True:
        _store(
            coll=coll,
            dout=dout,
            dref=dref,
        )

    if returnas is True:
        if return_dref is True:
            return dout, dref
        else:
            return dout


# #################################################################
# #################################################################
#               check
# #################################################################


def _check(
    coll=None,
    key=None,
    levels=None,
    # store vs return
    returnas=None,
    return_dref=None,
    store=None,
    key_npts=None,
    key_lvls=None,
    key_cont0=None,
    key_cont1=None,
):

    # ------
    # key

    dp2d = coll.get_profiles2d()
    key = ds._generic_check._check_var(
        key, 'key',
        types=str,
        allowed=list(dp2d.keys()),
    )

    wm = coll._which_mesh
    wbs = coll._which_bsplines
    keybs = dp2d[key]
    keym0 = coll.dobj[wbs][keybs][wm]
    submesh = coll.dobj[wm][keym0]['submesh']
    if submesh is not None:
        keym0 = submesh

    # ----------
    # store

    store = ds._generic_check._check_var(
        store, 'store',
        types=bool,
        default=False,
    )

    # ----------
    # returnas

    returnas = ds._generic_check._check_var(
        returnas, 'returnas',
        types=bool,
        default=not store,
    )

    # ----------
    # return_dref

    return_dref = ds._generic_check._check_var(
        return_dref, 'return_dref',
        types=bool,
        default=False,
    )

    # -----------------
    # keys for storing

    # key_npts
    if store is True:
        lout = list(coll.dref.keys())
    else:
        lout = []
    key_npts = ds._generic_check._check_var(
        key_npts, 'key_npts',
        types=str,
        default=f'{key}_cont_npts',
        excluded=lout,
    )

    # key_lvls
    lout.append(key_npts)
    key_lvls = ds._generic_check._check_var(
        key_lvls, 'key_lvls',
        types=str,
        default=f'{key}_cont_nlvls',
        excluded=lout,
    )

    # key_cont0
    if store is True:
        lout = list(coll.ddata.keys())
    else:
        lout = []
    key_cont0 = ds._generic_check._check_var(
        key_cont0, 'key_cont0',
        types=str,
        default=f'{key}_cont0',
        excluded=lout,
    )

    # key_cont0
    lout.append(key_cont0)
    key_cont1 = ds._generic_check._check_var(
        key_cont1, 'key_cont1',
        types=str,
        default=f'{key}_cont1',
        excluded=lout,
    )

    return (
        key, keybs, keym0, levels,
        store, returnas, return_dref,
        key_npts, key_lvls, key_cont0, key_cont1,
    )


# #################################################################
# #################################################################
#               compute contours
# #################################################################


def _get_contours(
    key=None,
    xx0=None,
    xx1=None,
    val=None,
    axis=None,
    levels=None,
    npts=None,
    largest=None,
):
    """ Return x0, x1 coordinates of contours (time-dependent)

    For contourpy algorithm, the dimensions shoud be (ny, nx), from meshgrid

    xx0 = (n1, n0)
    xx1 = (n1, n0)
    val = (nt, n1, n0)
    levels = (nlevels,)

    c0 = (nt, nlevels, nmax) array of x0 coordinates
    c1 = (nt, nlevels, nmax) array of x1 coordinates

    The contour coordinates are uniformized to always have the same nb of pts

    DISCLAIMER: not entirely dimension agnostic yet!

    """

    # -------------
    # check inputs

    (
        axis, n0, n1, shape_other, shape_cont, npts, levels, largest,
    ) = _compute_check(
        xx0=xx0,
        xx1=xx1,
        val=val,
        axis=axis,
        levels=levels,
        npts=npts,
        largest=largest,
    )

    # --------
    # prepare

    # compute contours at rknots
    # see https://github.com/matplotlib/matplotlib/blob/main/src/_contour.h

    cont0 = np.full(shape_cont, np.nan)
    cont1 = np.full(shape_cont, np.nan)

    # prepare slices
    sli_v = [slice(None) if ii in axis else 0 for ii in range(val.ndim)]
    sli_c = [slice(None) if ii == axis[0] else 0 for ii in range(cont0.ndim)]

    # prepare indices
    ind_v = [ii for ii in range(val.ndim) if ii not in axis]
    ind_c = [ii for ii in range(cont0.ndim) if ii not in axis]

    # ------------------------
    # loop on slices

    maxnpts = 0
    maxlvl = levels[0]
    for ind in itt.product(*[range(aa) for aa in shape_other]):

        # slice
        for ii, jj in enumerate(ind):
            sli_v[ind_v[ii]] = jj
            sli_c[ind_c[ii]] = jj

        # define map
        contgen = contour_generator(
            x=xx0,
            y=xx1,
            z=val[tuple(sli_v)],
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

        # loop on levels
        for jj in range(len(levels)):

            # get contour
            no_cont, cj = _get_contours_lvls(
                contgen=contgen,
                level=levels[jj],
                largest=largest,
            )

            # if contour was found
            if no_cont is False:

                # maxnpts
                if cj.shape[0] > maxnpts:
                    maxnpts = cj.shape[0]
                    maxlvl = levels[jj]

                # slice
                sli_c[axis[1]] = jj

                # interpolate on desired nb of points
                cont0[tuple(sli_c)] = np.interp(
                    np.linspace(0, cj.shape[0], npts),
                    np.arange(0, cj.shape[0]),
                    cj[:, 0],
                )
                cont1[tuple(sli_c)] = np.interp(
                    np.linspace(0, cj.shape[0], npts),
                    np.arange(0, cj.shape[0]),
                    cj[:, 1],
                )

    # -----------------
    # Warning + return

    msg = (
        f"\nContour computing for '{key}':\n"
        f"\t- npts {npts} vs {maxnpts} (level {maxlvl})"
    )
    if maxnpts > npts:
        warnings.warn(msg)
    else:
        print(msg)

    return cont0, cont1


def _get_contours_lvls(
    contgen=None,
    level=None,
    largest=None,
):
    # compute concatenated contour
    no_cont = False
    cj = contgen.lines(level)

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

    # if one or several contours exist
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

        elif np.sum(np.all(~np.isnan(cj), axis=1)) < 3:
            no_cont = True
    else:
        no_cont = True

    return no_cont, cj


def _compute_check(
    xx0=None,
    xx1=None,
    val=None,
    axis=None,
    levels=None,
    npts=None,
    largest=None,
):

    # ---------------------
    # check axis and shapes

    if axis is None:
        axis = [
            ii for ii, ss in enumerate(val.shape)
            if ss in xx0.shape
        ]
    assert xx0.ndim == 2
    assert xx0.shape == xx1.shape
    assert len(axis) == 2

    # ---------------------
    # check axis and shapes

    npts = ds._generic_check._check_var(
        npts, 'npts',
        types=int,
        default=500,
        sign='>0',
    )

    # ----------
    # levels

    levels = ds._generic_check._check_flat1darray(
        levels, 'levels',
        dtype=float,
    )

    # ---------------------
    # check axis and shapes

    n0, n1 = xx0.shape
    shape_other = [
        ss for ii, ss in enumerate(val.shape)
        if ii not in axis
    ]
    shape_cont = tuple([
        npts if ii == axis[0]
        else (levels.size if ii == axis[1] else ss)
        for ii, ss in enumerate(val.shape)
    ])

    # ----------
    # largest

    largest = ds._generic_check._check_var(
        largest, 'largest',
        types=bool,
        default=True,
    )

    return axis, n0, n1, shape_other, shape_cont, npts, levels, largest


# #################################################################
# #################################################################
#               format
# #################################################################


def _format(
    coll=None,
    key=None,
    keybs=None,
    dinterp=None,
    cont0=None,
    cont1=None,
    axis=None,
    keym0=None,
    key_npts=None,
    key_lvls=None,
    key_cont0=None,
    key_cont1=None,
):

    npts = cont0.shape[axis[0]]
    nlvls = cont0.shape[axis[1]]

    # --------
    # ref

    ref = tuple([
        key_npts if ii == axis[0]
        else (key_lvls if ii == axis[1] else rr)
        for ii, rr in enumerate(dinterp['ref'])
    ])

    # -------
    # units

    wm = coll._which_mesh
    knots0, knots1 = coll.dobj[wm][keym0]['knots']

    lf = ['units', 'dim', 'quant', 'name']
    dd0 = {k0: coll.ddata[knots0][k0] for k0 in lf}
    dd1 = {k0: coll.ddata[knots1][k0] for k0 in lf}

    # -------------
    # populate dict

    dref = {
        'npts': {'key': key_npts, 'size': npts},
        'levels': {'key': key_lvls, 'size': nlvls},
    }

    dout = {
        'cont0': {
            'key': key_cont0,
            'data': cont0,
            'ref': ref,
            **dd0,
        },
        'cont1': {
            'key': key_cont1,
            'data': cont1,
            'ref': ref,
            **dd1,
        },
    }

    return dout, dref


# #################################################################
# #################################################################
#               store
# #################################################################


def _store(
    coll=None,
    dout=None,
    dref=None,
):

    # --------
    # add ref

    for k0, v0 in dref.items():
        coll.add_ref(**v0)

    # --------
    # add data

    for k0, v0 in dout.items():
        coll.add_data(**v0)


# #############################################################################
# #############################################################################
#                   Polygon simplification
# #############################################################################


# def _simplify_polygon(pR=None, pZ=None, res=None):
    # """ Use convex hull with a constraint on the maximum discrepancy """

    # # ----------
    # # preliminary 1: check there is non redundant point

    # dp = np.sqrt((pR[1:] - pR[:-1])**2 + (pZ[1:] - pZ[:-1])**2)
    # ind = (dp > 1.e-6).nonzero()[0]
    # pR = pR[ind]
    # pZ = pZ[ind]

    # # check new poly is closed
    # if (pR[0] != pR[-1]) or (pZ[0] != pZ[-1]):
        # pR = np.append(pR, pR[0])
        # pZ = np.append(pZ, pZ[0])

    # # check it is counter-clockwise
    # clock = np.nansum((pR[1:] - pR[:-1]) * (pZ[1:] + pZ[:-1]))
    # if clock > 0:
        # pR = pR[::-1]
        # pZ = pZ[::-1]

    # # threshold = diagonal of resolution + 10%
    # thresh = res * np.sqrt(2) * 1.1

    # # ----------
    # # preliminary 2: get convex hull and copy

    # poly = np.array([pR, pZ]).T
    # iconv = ConvexHull(poly, incremental=False).vertices

    # # close convex hull to iterate on edges
    # pR_conv = np.append(pR[iconv], pR[iconv[0]])
    # pZ_conv = np.append(pZ[iconv], pZ[iconv[0]])

    # # copy to create new polygon that will serve as buffer
    # pR_bis, pZ_bis = np.copy(pR), np.copy(pZ)

    # # -------------------------
    # # loop on convex hull edges

    # for ii in range(pR_conv.size - 1):

        # pR1, pR2 = pR_conv[ii], pR_conv[ii+1]
        # pZ1, pZ2 = pZ_conv[ii], pZ_conv[ii+1]
        # i0 = np.argmin(np.hypot(pR_bis - pR1, pZ_bis - pZ1))

        # # make sure it starts from p1
        # pR_bis = np.append(pR_bis[i0:], pR_bis[:i0])
        # pZ_bis = np.append(pZ_bis[i0:], pZ_bis[:i0])

        # # get indices of closest points to p1, p2
        # i1 = np.argmin(np.hypot(pR_bis - pR1, pZ_bis - pZ1))
        # i2 = np.argmin(np.hypot(pR_bis - pR2, pZ_bis - pZ2))

        # # get corresponding indices of poly points to be included
        # if i2 == i1 + 1:
            # itemp = [i1, i2]

        # else:
            # # several points in-between
            # # => check they are all within distance before exclusing them

            # # get unit vector of segment
            # norm12 = np.hypot(pR2 - pR1, pZ2 - pZ1)
            # u12R = (pR2 - pR1) / norm12
            # u12Z = (pZ2 - pZ1) / norm12

            # # get points standing between p1 nd p2
            # lpR = pR_bis[i1 + 1:i2]
            # lpZ = pZ_bis[i1 + 1:i2]

            # # indices of points standing too far from edge (use cross-product)
            # iout = np.abs(u12R*(lpZ - pZ1) - u12Z*(lpR - pR1)) > thresh

            # # if any pts too far => include all pts
            # if np.any(iout):
                # itemp = np.arange(i1, i2 + 1)
            # else:
                # itemp = [i1, i2]

        # # build pts_in
        # pR_in = pR_bis[itemp]
        # pZ_in = pZ_bis[itemp]

        # # concatenate to add to new polygon
        # pR_bis = np.append(pR_in, pR_bis[i2 + 1:])
        # pZ_bis = np.append(pZ_in, pZ_bis[i2 + 1:])

    # # check new poly is closed
    # if (pR_bis[0] != pR_bis[-1]) or (pZ_bis[0] != pZ_bis[-1]):
        # pR_bis = np.append(pR_bis, pR_bis[0])
        # pZ_bis = np.append(pZ_bis, pZ_bis[0])

    # return pR_bis, pZ_bis


# #############################################################################
# #############################################################################
#                   radius2d special points handling
# #############################################################################


# def radius2d_special_points(
    # coll=None,
    # key=None,
    # keym0=None,
    # res=None,
# ):

    # keybs = coll.ddata[key]['bsplines']
    # keym = coll.dobj['bsplines'][keybs]['mesh']
    # mtype = coll.dobj[coll._which_mesh][keym]['type']
    # assert mtype in ['rect', 'tri']

    # # get map sampling
    # RR, ZZ = coll.get_sample_mesh(
        # key=keym,
        # res=res,
        # grid=True,
    # )

    # # get map
    # val, t, _ = coll.interpolate_profile2d(
        # key=key,
        # R=RR,
        # Z=ZZ,
        # grid=False,
        # imshow=True,        # for contour
    # )

    # # get min max values
    # rmin = np.nanmin(val)
    # rmax = np.nanmax(val)

    # # get contour of 0
    # cR, cZ = _get_contours(
        # RR=RR,
        # ZZ=ZZ,
        # val=val,
        # levels=[rmin + 0.05*(rmax-rmin)],
    # )

    # # dref
    # ref_O = f'{keym0}-pts-O-n'
    # dref = {
        # ref_O: {'size': 1},
    # }

    # # get barycenter
    # if val.ndim == 3:
        # assert cR.shape[1] == 1
        # ax_R = np.nanmean(cR[:, 0, :], axis=-1)[:, None]
        # ax_Z = np.nanmean(cZ[:, 0, :], axis=-1)[:, None]
        # reft = coll.ddata[key]['ref'][0]
        # ref = (reft, ref_O)
    # else:
        # ax_R = np.r_[np.nanmean(cR)]
        # ax_Z = np.r_[np.nanmean(cZ)]
        # ref = (ref_O,)

    # kR = f'{keym0}-pts-O-R'
    # kZ = f'{keym0}-pts-O-Z'
    # ddata = {
        # kR: {
            # 'ref': ref,
            # 'data': ax_R,
            # 'dim': 'distance',
            # 'quant': 'R',
            # 'name': 'O-points_R',
            # 'units': 'm',
        # },
        # kZ: {
            # 'ref': ref,
            # 'data': ax_Z,
            # 'dim': 'distance',
            # 'quant': 'Z',
            # 'name': 'O-points_Z',
            # 'units': 'm',
        # },
    # }

    # return dref, ddata, kR, kZ


# #############################################################################
# #############################################################################
#                   angle2d discontinuity handling
# #############################################################################


# def angle2d_zone(
    # coll=None,
    # key=None,
    # keyrad2d=None,
    # key_ptsO=None,
    # res=None,
    # keym0=None,
# ):

    # keybs = coll.ddata[key]['bsplines']
    # keym = coll.dobj['bsplines'][keybs]['mesh']
    # mtype = coll.dobj[coll._which_mesh][keym]['type']
    # assert mtype in ['rect', 'tri']

    # # --------------
    # # prepare

    # hastime, hasvect, reft, keyt = coll.get_time(key=key)[:4]
    # if hastime:
        # nt = coll.dref[reft]['size']
    # else:
        # msg = (
            # "Non time-dependent angle2d not implemented yet\n"
            # "=> ping @Didou09 on Github to open an issue"
        # )
        # raise NotImplementedError(msg)

    # if res is None:
        # res = _get_sample_mesh_res(
            # coll=coll,
            # keym=keym,
            # mtype=mtype,
        # )

    # # get map sampling
    # RR, ZZ = coll.get_sample_mesh(
        # key=keym,
        # res=res/2.,
        # grid=True,
        # imshow=True,    # for contour
    # )

    # # get map
    # val, t, _ = coll.interpolate_profile2d(
        # key=key,
        # R=RR,
        # Z=ZZ,
        # grid=False,
        # azone=False,
    # )
    # val[np.isnan(val)] = 0.
    # amin = np.nanmin(val)
    # amax = np.nanmax(val)

    # # get contours of absolute value
    # cRmin, cZmin = _get_contours(
        # RR=RR,
        # ZZ=ZZ,
        # val=val,
        # levels=[amin + 0.10*(amax - amin)],
        # largest=True,
        # uniform=True,
    # )
    # cRmax, cZmax = _get_contours(
        # RR=RR,
        # ZZ=ZZ,
        # val=val,
        # levels=[amax - 0.10*(amax - amin)],
        # largest=True,
        # uniform=True,
    # )

    # cRmin, cZmin = cRmin[:, 0, :], cZmin[:, 0, :]
    # cRmax, cZmax = cRmax[:, 0, :], cZmax[:, 0, :]

    # rmin = np.full(cRmin.shape, np.nan)
    # rmax = np.full(cRmax.shape, np.nan)

    # # get points inside contour
    # for ii in range(nt):
        # rmin[ii, :], _, _ = coll.interpolate_profile2d(
            # key=keyrad2d,
            # R=cRmin[ii, :],
            # Z=cZmin[ii, :],
            # grid=False,
            # indt=ii,
        # )
        # rmax[ii, :], _, _ = coll.interpolate_profile2d(
            # key=keyrad2d,
            # R=cRmax[ii, :],
            # Z=cZmax[ii, :],
            # grid=False,
            # indt=ii,
        # )

    # # get magnetic axis
    # kR, kZ = key_ptsO
    # axR = coll.ddata[kR]['data']
    # axZ = coll.ddata[kZ]['data']
    # assert coll.ddata[kR]['ref'][0] == coll.ddata[key]['ref'][0]

    # start_min = np.nanargmin(rmin, axis=-1)
    # start_max = np.nanargmin(rmax, axis=-1)

    # # re-order from start_min, start_max
    # lpR, lpZ = [], []
    # for ii in range(rmin.shape[0]):
        # imin = np.r_[
            # np.arange(start_min[ii], rmin.shape[1]),
            # np.arange(0, start_min[ii]),
        # ]

        # cRmin[ii] = cRmin[ii, imin]
        # cZmin[ii] = cZmin[ii, imin]
        # rmin[ii] = rmin[ii, imin]
        # # check it is counter-clockwise
        # clock = np.nansum(
            # (cRmin[ii, 1:] - cRmin[ii, :-1])
            # *(cZmin[ii, 1:] + cZmin[ii, :-1])
        # )
        # if clock > 0:
            # cRmin[ii, :] = cRmin[ii, ::-1]
            # cZmin[ii, :] = cZmin[ii, ::-1]
            # rmin[ii, :] = rmin[ii, ::-1]

        # imax = np.r_[
            # np.arange(start_max[ii], rmax.shape[1]),
            # np.arange(0, start_max[ii])
        # ]
        # cRmax[ii] = cRmax[ii, imax]
        # cZmax[ii] = cZmax[ii, imax]
        # rmax[ii] = rmax[ii, imax]
        # # check it is clockwise
        # clock = np.nansum(
            # (cRmax[ii, 1:] - cRmax[ii, :-1])
            # *(cZmax[ii, 1:] + cZmax[ii, :-1])
        # )
        # if clock < 0:
            # cRmax[ii, :] = cRmax[ii, ::-1]
            # cZmax[ii, :] = cZmax[ii, ::-1]
            # rmax[ii, :] = rmax[ii, ::-1]

        # # i0
        # dr = np.diff(rmin[ii, :])
        # i0 = (np.isnan(dr) | (dr < 0)).nonzero()[0][0]
        # # rmin[ii, i0-1:] = np.nan
        # dr = np.diff(rmax[ii, :])
        # i1 = (np.isnan(dr) | (dr < 0)).nonzero()[0][0]
        # # rmax[ii, i1-1:] = np.nan

        # # polygon
        # pR = np.r_[axR[ii], cRmin[ii, :i0-1], cRmax[ii, :i1-1][::-1]]
        # pZ = np.r_[axZ[ii], cZmin[ii, :i0-1], cZmax[ii, :i1-1][::-1]]

        # pR, pZ = _simplify_polygon(pR=pR, pZ=pZ, res=res)

        # lpR.append(pR)
        # lpZ.append(pZ)

    # # Ajust sizes
    # nb = np.array([pR.size for pR in lpR])

    # #
    # nmax = np.max(nb)
    # pR = np.full((nt, nmax), np.nan)
    # pZ = np.full((nt, nmax), np.nan)

    # for ii in range(nt):
        # pR[ii, :] = np.interp(
            # np.linspace(0, nb[ii], nmax),
            # np.arange(0, nb[ii]),
            # lpR[ii],
        # )
        # pZ[ii, :] = np.interp(
            # np.linspace(0, nb[ii], nmax),
            # np.arange(0, nb[ii]),
            # lpZ[ii],
        # )

    # # ----------------
    # # prepare output dict

    # # ref
    # kref = f'{keym0}-azone-npt'
    # dref = {
        # kref: {'size': nmax}
    # }

    # # data
    # kR = f'{keym0}-azone-R'
    # kZ = f'{keym0}-azone-Z'
    # ddata = {
        # kR: {
            # 'data': pR,
            # 'ref': (reft, kref),
            # 'units': 'm',
            # 'dim': 'distance',
            # 'quant': 'R',
            # 'name': None,
        # },
        # kZ: {
            # 'data': pZ,
            # 'ref': (reft, kref),
            # 'units': 'm',
            # 'dim': 'distance',
            # 'quant': 'R',
            # 'name': None,
        # },
    # }

    # return dref, ddata, kR, kZ


# def angle2d_inzone(
    # coll=None,
    # keym0=None,
    # keya2d=None,
    # R=None,
    # Z=None,
    # t=None,
    # indt=None,
# ):


    # # ------------
    # # prepare points

    # if R.ndim == 1:
        # shape0 = None
        # pts = np.array([R, Z]).T
    # else:
        # shape0 = R.shape
        # pts = np.array([R.ravel(), Z.ravel()]).T

    # # ------------
    # # prepare path

    # kazR, kazZ = coll.dobj[coll._which_mesh][keym0]['azone']
    # pR = coll.ddata[kazR]['data']
    # pZ = coll.ddata[kazZ]['data']

    # hastime, hasvect, reft, keyt, tnew, dind = coll.get_time(
        # key=kazR,
        # t=t,
        # indt=indt,
    # )

    # # ------------
    # # test points

    # if hastime:
        # if dind is None:
            # nt = coll.dref[reft]['size']
            # ind = np.zeros((nt, R.size), dtype=bool)
            # for ii in range(nt):
                # path = Path(np.array([pR[ii, :], pZ[ii, :]]).T)
                # ind[ii, :] = path.contains_points(pts)
        # else:
            # import pdb; pdb.set_trace()     # DB
            # raise NotImplementedError()
            # # TBC / TBF
            # nt = None
            # ind = np.zeros((nt, R.size), dtype=bool)
            # for ii in range(nt):
                # path = Path(np.array([pR[ii, :], pZ[ii, :]]).T)
                # ind[ii, :] = path.contains_points(pts)

    # else:
        # path = Path(np.array([pR, pZ]).T)
        # ind = path.contains_points(pts)

    # # -------------------------
    # # fromat output and return

    # if shape0 is not None:
        # if hastime:
            # ind = ind.reshape(tuple(np.r_[nt, shape0]))
        # else:
            # ind = ind.reshape(shape0)

    # return ind