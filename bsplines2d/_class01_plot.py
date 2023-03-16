# -*- coding: utf-8 -*-


# Built-in


# Common
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import datastock as ds

# specific
from . import _generic_check


# #############################################################################
# #############################################################################
#                     Main plotting
# #############################################################################


def plot_mesh(
    coll=None,
    key=None,
    ind_knot=None,
    ind_cent=None,
    crop=None,
    bck=None,
    nmax=None,
    color=None,
    dax=None,
    dmargin=None,
    fs=None,
    dleg=None,
    connect=None,
):
    """ Plot the desired mesh

    rect and tri meshes are constant
    polar meshes can vary in time

    """

    # --------------
    # check input

    (
     key, nd, mtype, ind_knot, ind_cent, crop, bck, color, dleg,
     return_neighbours,
     ) = _plot_mesh_check(
        coll=coll,
        key=key,
        ind_knot=ind_knot,
        ind_cent=ind_cent,
        crop=crop,
        bck=bck,
        color=color,
        dleg=dleg,
    )

    # ------------------------
    # call appropriate routine

    if nd == '1d':
        return plot_mesh_1d(
            coll=coll,
            key=key,
            ind_knot=ind_knot,
            ind_cent=ind_cent,
            #units=units,
            return_neighbours=return_neighbours,
            nmax=nmax,
            color=color,
            dax=dax,
            dmargin=dmargin,
            fs=fs,
            dleg=dleg,
            connect=connect,
        )

    elif mtype in ['rect', 'tri']:
        # time-fixed meshes
        return _plot_mesh_2d_recttri(
            coll=coll,
            key=key,
            mtype=mtype,
            ind_knot=ind_knot,
            ind_cent=ind_cent,
            return_neighbours=return_neighbours,
            crop=crop,
            bck=bck,
            color=color,
            dax=dax,
            fs=fs,
            dmargin=dmargin,
            dleg=dleg,
        )

    else:
        # possibly time-varying mesh
        raise NotImplementedError()
        # return _plot_mesh_2d_polar(
            # coll=coll,
            # key=key,
            # nmax=nmax,
            # color=color,
            # dax=dax,
            # fs=fs,
            # dmargin=dmargin,
            # dleg=dleg,
            # connect=connect,
        # )


# #############################################################################
# #############################################################################
#                           checks
# #############################################################################


def _plot_mesh_check(
    coll=None,
    key=None,
    ind_knot=None,
    ind_cent=None,
    crop=None,
    bck=None,
    color=None,
    dleg=None,
):

    # key
    key = ds._generic_check._check_var(
        key, 'key',
        default=None,
        types=str,
        allowed=list(coll.dobj.get(coll._which_mesh, {}).keys()),
    )

    nd = coll.dobj[coll._which_mesh][key]['nd']
    mtype = coll.dobj[coll._which_mesh][key]['type']

    # crop, bck
    crop = ds._generic_check._check_var(crop, 'crop', default=True, types=bool)
    bck = ds._generic_check._check_var(bck, 'bck', default=True, types=bool)

    if mtype in ['rect', 'tri']:
        return_neighbours = True
    else:
        return_neighbours = False

    # ind_knot
    if ind_knot is not None:
        ind_knot = coll.select_mesh_elements(
            key=key, ind=ind_knot, elements='knots',
            returnas='data', return_neighbours=return_neighbours, crop=crop,
        )
        if return_neighbours is False:
            ind_knot = [ind_knot]

    # ind_cent
    if ind_cent is not None:
        ind_cent = coll.select_mesh_elements(
            key=key, ind=ind_cent, elements='cents',
            returnas='data', return_neighbours=return_neighbours, crop=crop,
        )
        if return_neighbours is False:
            ind_cent = [ind_cent]

    # color
    if color is None:
        color = 'k'
    if not mcolors.is_color_like(color):
        msg = (
            "Arg color must be a valid matplotlib color identifier!\n"
            f"Provided: {color}"
        )
        raise Exception(msg)

    # dleg
    defdleg = {
        'bbox_to_anchor': (1.1, 1.),
        'loc': 'upper left',
        'frameon': True,
    }
    dleg = ds._generic_check._check_var(
        dleg, 'dleg',
        default=defdleg,
        types=(bool, dict),
    )

    return (
        key, nd, mtype, ind_knot, ind_cent, crop, bck, color, dleg,
        return_neighbours,
    )


# #############################################################################
# #############################################################################
#                           prepare
# #############################################################################


def _plot_mesh_prepare_1d(
    coll=None,
    key=None,
    **kwd,
):

    # --------
    # prepare

    kknots = coll.dobj[coll._which_mesh][key]['knots'][0]
    knots = coll.ddata[kknots]['data']

    xx = np.array([knots, knots, np.full(knots.shape, np.nan)]).T.ravel()
    yy = np.array([
        np.zeros(knots.shape),
        np.ones(knots.shape),
        np.ones(knots.shape),
    ]).T.ravel()

    return xx, yy


def _plot_mesh_prepare_2d_rect(
    coll=None,
    key=None,
    crop=None,
    bck=None,
):

    # --------
    # prepare

    Rk, Zk = coll.dobj[coll._which_mesh][key]['knots']
    R = coll.ddata[Rk]['data']
    Z = coll.ddata[Zk]['data']

    vert = np.array([
        np.repeat(R, 3),
        np.tile((Z[0], Z[-1], np.nan), R.size),
    ])
    hor = np.array([
        np.tile((R[0], R[-1], np.nan), Z.size),
        np.repeat(Z, 3),
    ])

    # --------
    # compute

    grid_bck = None
    if crop is False or coll.dobj[coll._which_mesh][key]['crop'] is False:
        grid = np.concatenate((vert, hor), axis=1)

    else:

        crop = coll.ddata[coll.dobj[coll._which_mesh][key]['crop']]['data']

        grid = []
        icropR = np.r_[range(R.size-1), R.size-2]
        jcropZ = np.r_[range(Z.size-1), Z.size-2]

        # vertical lines  TBC
        for ii, ic in enumerate(icropR):
            if np.any(crop[ic, :]):
                if ii in [0, R.size-1]:
                    cropi = crop[ic, :]
                else:
                    cropi = crop[ic, :] | crop[ic-1, :]
                lseg = []
                for jj, jc in enumerate(jcropZ):
                    if jj == 0 and cropi[jc]:
                        lseg.append(Z[jj])
                    elif jj == Z.size-1 and cropi[jc]:
                        lseg.append(Z[jj])
                    elif cropi[jc] and not cropi[jc-1]:
                        if len(lseg) > 0:
                            lseg.append(np.nan)
                        lseg.append(Z[jj])
                    elif (not cropi[jc]) and cropi[jc-1]:
                        lseg.append(Z[jc])
                grid.append(np.concatenate(
                    (
                        np.array([R[ii]*np.ones((len(lseg),)), lseg]),
                        np.full((2, 1), np.nan)
                    ),
                    axis=1,
                ))

        # horizontal lines
        for jj, jc in enumerate(jcropZ):
            if np.any(crop[:, jc]):
                if jj in [0, Z.size-1]:
                    cropj = crop[:, jc]
                else:
                    cropj = crop[:, jc] | crop[:, jc-1]
                lseg = []
                for ii, ic in enumerate(icropR):
                    if ii in [0, R.size-1] and cropj[ic]:
                        lseg.append(R[ii])
                    elif cropj[ic] and not cropj[ic-1]:
                        if len(lseg) > 0:
                            lseg.append(np.nan)
                        lseg.append(R[ii])
                    elif (not cropj[ic]) and cropj[ic-1]:
                        lseg.append(R[ic])
                grid.append(np.concatenate(
                    (
                        np.array([lseg, Z[jj]*np.ones((len(lseg),))]),
                        np.full((2, 1), np.nan)
                    ),
                    axis=1,
                ))

        grid = np.concatenate(tuple(grid), axis=1)

        if bck is True:
            grid_bck = np.concatenate((vert, hor), axis=1)

    return grid, grid_bck


def _plot_mesh_prepare_2d_tri(
    coll=None,
    key=None,
    crop=None,
    bck=None,
):

    # --------
    # prepare

    grid_bck = None

    kknots = coll.dobj[coll._which_mesh][key]['knots']
    R = coll.ddata[kknots[0]]['data']
    Z = coll.ddata[kknots[1]]['data']

    indtri = coll.ddata[coll.dobj[coll._which_mesh][key]['ind']]['data']

    # find unique segments from all triangles
    segs = np.unique(
        np.sort(np.concatenate(
            (indtri[:, 0:2], indtri[:, 1:], indtri[:, ::2]),
            axis=0,
        )),
        axis=0,
    )

    # build long segments if possible
    ind = np.ones((segs.shape[0],), dtype=bool)
    ind[0] = False
    lseg = [segs[0, :]]
    last = segs[0, :]
    while np.any(ind):
        ii = segs[ind, 0] == last[-1]
        if np.any(ii):
            ii = ind.nonzero()[0][ii]
            dR0 = R[last[1]] - R[last[0]]
            dZ0 = Z[last[1]] - Z[last[0]]
            dR = np.diff(R[segs[ii, :]], axis=1)[:, 0]
            dZ = np.diff(Z[segs[ii, :]], axis=1)[:, 0]
            norm0 = np.sqrt(dR0**2 + dZ0**2)
            norm = np.sqrt(dR**2 + dZ**2)
            sca = (dR0*dR + dZ0*dZ) / (norm0 * norm)
            iwin = ii[np.argmax(sca)]
            lseg.append([segs[iwin, 1]])

        else:
            lseg.append([-1])
            iwin = ind.nonzero()[0][0]
            lseg.append(segs[iwin, :])

        last = segs[iwin, :]
        ind[iwin] = False

    lseg = np.concatenate(lseg)
    grid = np.array([R[lseg], Z[lseg]])
    grid[0, lseg == -1] = np.nan

    return grid, grid_bck


# def _plot_mesh_prepare_polar_cont(
    # coll=None,
    # key=None,
    # k2d=None,
    # RR=None,
    # ZZ=None,
    # ind=None,
    # nn=None,
# ):

    # # ---------------------
    # # sample mesh if needed

    # # ---------------------
    # # get map of rr / angle

    # if callable(k2d):

        # # check RR
        # if RR is None:
            # msg = (
                # "radius2d / angle2d are callable => provide RR and ZZ!"
            # )
            # raise Exception(msg)

        # # compute map
        # rr = k2d(RR, ZZ)[None, ...]
        # assert rr.ndim == RR.ndim + 1
        # reft = None
        # nt = 1

        # if nn is None:
            # nn = 50

        # # create vector
        # rad = np.linspace(np.nanmin(rr), np.nanmax(rr), nn)

    # else:
        # kn = coll.dobj[coll._which_mesh][key]['knots'][ind]
        # rad = coll.ddata[kn]['data']
        # kb2 = coll.ddata[k2d]['bsplines']

        # if RR is None:
            # km2 = coll.dobj['bsplines'][kb2]['mesh']
            # RR, ZZ = coll.get_sample_mesh(
                # key=km2,
                # res=None,
                # grid=True,
                # mode=None,
                # R=None,
                # Z=None,
                # DR=None,
                # DZ=None,
                # imshow=True,
            # )

        # rr = coll.interpolate_profile2d(
            # key=k2d,
            # R=RR,
            # Z=ZZ,
            # grid=False,
            # return_params=False,
        # )[0]

        # refr2d = coll.ddata[k2d]['ref']
        # refbs = coll.dobj['bsplines'][kb2]['ref']
        # if refr2d == refbs:
            # reft = None
            # nt = 1
            # rr = rr[None, ...]
        # elif len(refr2d) == len(refbs) + 1 and refr2d[1:] == refbs:
            # reft = refr2d[0]
            # nt = coll.dref[reft]['size']

    # assert rr.shape[0] == nt

    # # ----------------
    # # Compute contours

    # contR, contZ = _contours._get_contours(
        # xx0=RR,
        # xx1=ZZ,
        # val=rr,
        # levels=rad,
    # )

    # # refrad
    # refrad = coll.dobj[coll._which_mesh][key]['knots'][ind]

    # return contR, contZ, rad, reft, refrad, RR, ZZ


# def _plot_mesh_prepare_polar(
    # coll=None,
    # key=None,
    # # Necessary for callable radius2d
    # RR=None,
    # ZZ=None,
# ):

    # # --------
    # # prepare

    # # create rectangular grid and compute radius at each point
    # k2d = coll.dobj[coll._which_mesh][key]['radius2d']
    # (
        # contRrad, contZrad,
        # rad, reft, refrad,
        # RR, ZZ,
    # ) = _plot_mesh_prepare_polar_cont(
        # coll=coll,
        # key=key,
        # k2d=k2d,
        # RR=RR,
        # ZZ=ZZ,
        # ind=0,
        # nn=None,        # nrad if k2d callable
    # )

    # # -----------
    # # contour of angle if angle not None

    # contRang, contZang, ang, refang = None, None, None, None
    # if len(coll.dobj[coll._which_mesh][key]['shape-c']) == 2:
        # # create rectangular grid and compute radius at each point
        # k2d = coll.dobj[coll._which_mesh][key]['angle2d']
        # (
            # contRang, contZang,
            # ang, _, refang,
            # _, _,
        # ) = _plot_mesh_prepare_polar_cont(
            # coll=coll,
            # key=key,
            # k2d=k2d,
            # RR=RR,
            # ZZ=ZZ,
            # ind=1,
            # nn=None,        # nang if k2d callable
        # )

    # return (
        # contRrad, contZrad, rad, refrad,
        # contRang, contZang, ang, refang,
        # reft,
    # )


# #############################################################################
# #############################################################################
#                           mesh type specific
# #############################################################################


def plot_mesh_1d(
    coll=None,
    key=None,
    ind_knot=None,
    ind_cent=None,
    return_neighbours=None,
    units=None,
    nmax=None,
    color=None,
    dax=None,
    dmargin=None,
    fs=None,
    dleg=None,
    connect=None,
):
    """ Plot the desired spectral mesh

    """

    # --------------
    #  Prepare data

    xx, yy = _plot_mesh_prepare_1d(
        coll=coll,
        key=key,
    )

    # if units not in [None, 'eV']:
        # xx, _, _, cat = _spectralunits.convert_spectral(
            # data_in=xx,
            # units_in='eV',
            # units_out=units,
        # )
        # xlab = cat + r" ($" + units + "$)"

    # else:
        # xlab = r'energy ($eV$)'
    xlab = None


    # --------------
    # plot - prepare

    if dax is None:

        if dmargin is None:
            dmargin = {
                'left': 0.1, 'right': 0.9,
                'bottom': 0.1, 'top': 0.9,
                'hspace': 0.1, 'wspace': 0.1,
            }

        fig = plt.figure(figsize=fs)
        gs = gridspec.GridSpec(ncols=1, nrows=1, **dmargin)
        ax0 = fig.add_subplot(gs[0, 0])
        ax0.set_xlabel(xlab)

        dax = {'spectral': ax0}

    dax = _generic_check._check_dax(dax=dax, main='spectral')

    # --------------
    # plot

    kax = 'spectral'
    if dax.get(kax) is not None:
        ax = dax[kax]['handle']

        ax.plot(
            xx,
            yy,
            ls='-',
            lw=0.5,
            color=color,
            alpha=0.5,
            label=key,
        )

        if ind_knot is not None:
            ax.plot(
                ind_knot[0][0],
                0.5,
                marker='o',
                ms=8,
                ls='None',
                color=color,
                label='knots',
            )

            # if return_neighbours:
                # ax.plot(
                #     ind_knot[1][0, :, :],
                #     marker='x',
                #     ms=4,
                #     ls='None',
                #     color=color,
                # )

        if ind_cent is not None:
            ax.plot(
                ind_cent[0][0],
                0.5,
                marker='x',
                ms=8,
                ls='None',
                color=color,
                label='cents',
            )

            # if return_neighbours:
                # ax.plot(
                #     ind_cent[1][0, :, :],
                #     marker='o',
                #     ms=4,
                #     ls='None',
                #     color=color,
                # )

    # --------------
    # dleg

    if dleg is not False:
        for kax in dax.keys():
            dax[kax]['handle'].legend(**dleg)

    return dax


def _plot_mesh_2d_recttri(
    coll=None,
    key=None,
    mtype=None,
    ind_knot=None,
    ind_cent=None,
    return_neighbours=None,
    crop=None,
    bck=None,
    color=None,
    dax=None,
    fs=None,
    dmargin=None,
    dleg=None,
):

    # --------------
    #  Prepare data

    if mtype == 'rect':
        grid, grid_bck = _plot_mesh_prepare_2d_rect(
            coll=coll,
            key=key,
            crop=crop,
            bck=bck,
        )
    else:
        grid, grid_bck = _plot_mesh_prepare_2d_tri(
            coll=coll,
            key=key,
            crop=crop,
            bck=bck,
        )

    # --------------
    # plot - prepare

    if dax is None:

        if dmargin is None:
            dmargin = {
                'left': 0.1, 'right': 0.9,
                'bottom': 0.1, 'top': 0.9,
                'hspace': 0.1, 'wspace': 0.1,
            }

        fig = plt.figure(figsize=fs)
        gs = gridspec.GridSpec(ncols=1, nrows=1, **dmargin)
        ax0 = fig.add_subplot(gs[0, 0], aspect='equal')
        ax0.set_xlabel(f'R (m)')
        ax0.set_ylabel(f'Z (m)')

        dax = {'cross': ax0}

    dax = _generic_check._check_dax(dax=dax, main='cross')

    # --------------
    # plot

    kax = 'cross'
    if dax.get(kax) is not None:
        ax = dax[kax]['handle']

        if grid_bck is not None and bck is True:
            ax.plot(
                grid_bck[0, :],
                grid_bck[1, :],
                ls='-',
                lw=0.5,
                color=color,
                alpha=0.5,
                label=key,
            )

        ax.plot(
            grid[0, :],
            grid[1, :],
            color=color,
            ls='-',
            lw=1.,
            label=key,
        )

        if ind_knot is not None:

            ax.plot(
                ind_knot[0][0],
                ind_knot[0][1],
                marker='o',
                ms=8,
                ls='None',
                color=color,
                label='knots',
            )

            if return_neighbours:
                ax.plot(
                    ind_knot[1][0, :, :],
                    ind_knot[1][1, :, :],
                    marker='x',
                    ms=4,
                    ls='None',
                    color=color,
                )

        if ind_cent is not None:

            ax.plot(
                ind_cent[0][0],
                ind_cent[0][1],
                marker='x',
                ms=8,
                ls='None',
                color=color,
                label='cents',
            )

            if return_neighbours:
                ax.plot(
                    ind_cent[1][0, :, :],
                    ind_cent[1][1, :, :],
                    marker='o',
                    ms=4,
                    ls='None',
                    color=color,
                )

    # --------------
    # dleg

    if dleg is not False:
        for kax in dax.keys():
            dax[kax]['handle'].legend(**dleg)

    return dax


# def _plot_mesh_2d_polar(
    # coll=None,
    # key=None,
    # npts=None,
    # nmax=None,
    # color=None,
    # dax=None,
    # fs=None,
    # dmargin=None,
    # dleg=None,
    # connect=None,
# ):

    # # --------------
    # #  Prepare data

    # if nmax is None:
        # nmax = 2

    # (
        # contRrad, contZrad, rad, refrad,
        # contRang, contZang, ang, refang,
        # reft,
    # ) = _plot_mesh_prepare_polar(
        # coll=coll,
        # key=key,
    # )
    # refptsr = 'ptsr'
    # nt, nr, nptsr = contRrad.shape
    # if contRang is not None:
        # refptsa = 'ptsa'
        # _, nang, nptsa = contRang.shape

    # # --------------------
    # # Instanciate Plasma2D

    # coll2 = coll.__class__()

    # # ref
    # coll2.add_ref(
        # key=reft,
        # size=nt,
    # )
    # reft = list(coll2.dref.keys())[0]
    # coll2.add_ref(
        # key=refrad,
        # size=nr,
    # )
    # coll2.add_ref(
        # key=refptsr,
        # size=nptsr,
    # )

    # if contRang is not None:
        # coll2.add_ref(
            # key=refang,
            # size=nang,
        # )
        # coll2.add_ref(
            # key=refptsa,
            # size=nptsa,
        # )

    # # data
    # coll2.add_data(
        # key='radius',
        # data=rad,
        # ref=(refrad,)
    # )
    # coll2.add_data(
        # key='contRrad',
        # data=contRrad,
        # ref=(reft, refrad, refptsr)
    # )
    # coll2.add_data(
        # key='contZrad',
        # data=contZrad,
        # ref=(reft, refrad, refptsr)
    # )

    # if contRang is not None:
        # coll2.add_data(
            # key='angle',
            # data=ang,
            # ref=(refang,)
        # )
        # coll2.add_data(
            # key='contRang',
            # data=contRang,
            # ref=(reft, refang, refptsa)
        # )
        # coll2.add_data(
            # key='contZang',
            # data=contZang,
            # ref=(reft, refang, refptsa)
        # )

    # # -----
    # # plot

    # if contRang is None:
        # return coll2.plot_as_mobile_lines(
            # keyX='contRrad',
            # keyY='contZrad',
            # key_time=reft,
            # key_chan='radius',
            # connect=connect,
        # )

    # else:

        # daxrad, dgrouprad = coll2.plot_as_mobile_lines(
            # keyX='contRrad',
            # keyY='contZrad',
            # key_time=reft,
            # key_chan='radius',
            # connect=False,
            # inplace=False,
        # )

        # daxang, dgroupang = coll2.plot_as_mobile_lines(
            # keyX='contRang',
            # keyY='contZang',
            # key_time=reft,
            # key_chan='angle',
            # connect=False,
            # inplace=False,
        # )

        # # connect
        # if connect is False:
            # return (daxrad, daxang), (dgrouprad, dgroupang)

        # else:
            # daxrad.setup_interactivity(
                # kinter='inter0', dgroup=dgrouprad, dinc=None,
            # )
            # daxrad.disconnect_old()
            # daxrad.connect()

            # daxang.setup_interactivity(
                # kinter='inter0', dgroup=dgroupang, dinc=None,
            # )
            # daxang.disconnect_old()
            # daxang.connect()

            # daxrad.show_commands()
            # return daxrad, daxang
