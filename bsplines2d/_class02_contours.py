# -*- coding: utf-8 -*-


# Built-in
import warnings


# Common
import numpy as np
from scipy.spatial import ConvexHull
from matplotlib.path import Path
from contourpy import contour_generator
import datastock as ds


# specific


# #################################################################
# #################################################################
#               Contour computation
# #################################################################


def _get_contours(
    xx0=None,
    xx1=None,
    val=None,
    levels=None,
    largest=None,
    uniform=None,
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

    if largest is None:
        largest = False

    if uniform is None:
        uniform = True

    # val.shape = (nt, nR, nZ)
    lc = [
        val.shape == xx0.shape,
        val.ndim == xx0.ndim + 1 and val.shape[1:] == xx0.shape,
    ]
    if lc[0]:
        val = val[None, ...]
    elif lc[1]:
        pass
    else:
        msg = "Incompatible val.shape!"
        raise Exception(msg)

    nt, n0, n1 = val.shape

    # ------------------------
    # Compute list of contours

    # compute contours at rknots
    # see https://github.com/matplotlib/matplotlib/blob/main/src/_contour.h

    cont0 = [[] for ii in range(nt)]
    cont1 = [[] for ii in range(nt)]
    for ii in range(nt):

        # define map
        contgen = contour_generator(
            x=xx0,
            y=xx1,
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

                elif np.sum(np.all(~np.isnan(cj), axis=1)) < 3:
                    no_cont = True
            else:
                no_cont = True

            if no_cont is True:
                cj = np.full((3, 2), np.nan)

            cont0[ii].append(cj[:, 0])
            cont1[ii].append(cj[:, 1])

    # ------------------------------------------------
    # Interpolate / concatenate to uniformize as array

    if uniform:
        ln = [[pp.size for pp in cc] for cc in cont0]
        nmax = np.max(ln)
        c0 = np.full((nt, len(levels), nmax), np.nan)
        c1 = np.full((nt, len(levels), nmax), np.nan)

        for ii in range(nt):
            for jj in range(len(levels)):
                c0[ii, jj, :] = np.interp(
                    np.linspace(0, ln[ii][jj], nmax),
                    np.arange(0, ln[ii][jj]),
                    cont0[ii][jj],
                )
                c1[ii, jj, :] = np.interp(
                    np.linspace(0, ln[ii][jj], nmax),
                    np.arange(0, ln[ii][jj]),
                    cont1[ii][jj],
                )

        return c0, c1
    else:
        return cont0, cont1


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
