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
from . import _class01_checks as _checks


# #############################################################################
# #############################################################################
#                           Main
# #############################################################################


def plot_as_profile2d(
    # ressources
    coll=None,
    # inputs
    key=None,
    # parameters
    dres=None,
    # contours
    levels=None,
    # figure
    vmin=None,
    vmax=None,
    cmap=None,
    dax=None,
    dmargin=None,
    fs=None,
    dcolorbar=None,
    dleg=None,
    # interactivity
    dinc=None,
    connect=None,
):

    # --------------
    # check input

    if connect is None:
        connect = True

    (
        key, keybs, keym, mtype,
        subbs, submtype,
        cmap, dcolorbar, dleg,
        connect,
    ) = _check(
        coll=coll,
        key=key,
        # figure
        cmap=cmap,
        dcolorbar=dcolorbar,
        dleg=dleg,
        # interactivity
        connect=connect,
    )

    # --------------
    #  Prepare data

    (
        coll2, dkeys, interp,
        key_cont, axis, dout, dref, refZ, refU,
    ) = _prepare(
        coll=coll,
        key=key,
        keybs=keybs,
        keym=keym,
        dres=dres,
        mtype=mtype,
        # submesh
        subbs=subbs,
        submtype=submtype,
        # levels
        levels=levels,
    )

    # ---------------
    # call right function

    if levels is None:
        return coll2.plot_as_array(
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            dax=dax,
            dmargin=dmargin,
            fs=fs,
            dcolorbar=dcolorbar,
            dleg=dleg,
            interp=interp,
            connect=connect,
            **dkeys,
        )

    else:
        dax, dgroup = coll2.plot_as_array(
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            dax=dax,
            dmargin=dmargin,
            fs=fs,
            dcolorbar=dcolorbar,
            dleg=dleg,
            interp=interp,
            connect=False,
            **dkeys,
        )

        # ---------------
        # add contour

        for k0, v0 in dout.items():
            sh = dout[k0]['data'].shape
            shnan = [1 if ii == axis[0] else ss for ii, ss in enumerate(sh)]

            dout[k0]['data'] = np.append(
                v0['data'],
                np.full(tuple(shnan), np.nan),
                axis=axis[0],
            )

            sh = dout[k0]['data'].shape
            newpts = sh[axis[0]]*sh[axis[1]]
            sh = tuple(np.r_[sh[:axis[0]], newpts, sh[axis[1]+1:]].astype(int))
            newref = tuple(np.r_[
                v0['ref'][:axis[0]],
                [dref['npts']['key']],
                v0['ref'][axis[1]+1:],
            ])

            dout[k0]['data'] = dout[k0]['data'].swapaxes(axis[0], axis[1]).reshape(sh)
            dout[k0]['ref'] = newref

        dref['npts']['size'] = newpts
        dax.add_ref(**dref['npts'])

        for k0, v0 in dout.items():
            dax.add_data(**v0)

        levels = np.atleast_1d(levels)

        # ---------------
        # prepare contour

        ndim = len(dgroup)
        cont0 = dax.ddata[key_cont[0]]['data']
        cont1 = dax.ddata[key_cont[1]]['data']

        if ndim == 3:
            axisZ = None
            sli = [
                slice(None) if ii == axis[0] else 0
                for ii in range(ndim-1)
            ]

        elif ndim == 4:
            raise NotImplementedError()


        # -----------
        # add contour

        kax = 'matrix'
        if dax.dax.get(kax) is not None:
            ax = dax.dax[kax]['handle']

            if ndim == 2:

                for ii in range(len(levels)):
                    ax.plot(
                        cont0,
                        cont1,
                        ls='-',
                        lw=1.,
                        c='k',
                    )

            elif ndim == 3:

                l0, = ax.plot(
                    cont0[sli],
                    cont1[sli],
                    ls='-',
                    lw=1.,
                    c='k',
                )

                km = f'contours'
                dax.add_mobile(
                    key=km,
                    handle=l0,
                    # group_vis='Z',
                    # refs=[(refZ, ref_lvls), (refZ, ref_lvls)],
                    refs=[(refZ,), (refZ,)],
                    data=[key_cont[0], key_cont[1]],
                    dtype=['xdata', 'ydata'],
                    axes=kax,
                    ind=0,
                )

            else:
                raise NotImplementedError()

                # sli = [slice(None)]

                # l0, = ax.plot(
                # )

                # k0 = f'contour{ii}'
                # dax.add_mobile(
                # )

        # -----------
        # connect

        if connect is True:
            dax.setup_interactivity(kinter='inter0', dgroup=dgroup, dinc=dinc)
            dax.disconnect_old()
            dax.connect()

            dax.show_commands()
            return dax
        else:
            return dax, dgroup


# #############################################################################
# #############################################################################
#                       Check
# #############################################################################


def _check(
    coll=None,
    key=None,
    # figure
    cmap=None,
    dcolorbar=None,
    dleg=None,
    # interactivity
    connect=None,
):

    # ----------
    # keys

    # key
    dk = coll.get_profiles2d()
    key = ds._generic_check._check_var(
        key, 'key',
        types=str,
        allowed=list(dk.keys()),
    )

    wm = coll._which_mesh
    wbs = coll._which_bsplines

    keybs = dk[key]
    refbs = coll.dobj[wbs][keybs]['ref']

    keym = coll.dobj[wbs][keybs][wm]
    mtype = coll.dobj[wm][keym]['type']

    submesh = coll.dobj[wm][keym]['submesh']
    if submesh == '':
        submesh = None

    if submesh is None:
        subbs = keybs
        submtype = mtype
    else:
        subbs = coll.dobj[wm][keym]['subbs']
        subkey = coll.dobj[wm][keym]['subkey']
        submtype = coll.dobj[wm][keym]['type']

    # ----------
    # derived

    # ----------
    # figure

    # cmap
    if cmap is None:
        cmap = 'viridis'

    # dcolorbar
    defdcolorbar = {
        # 'location': 'right',
        'fraction': 0.15,
        'orientation': 'vertical',
    }
    dcolorbar = ds._generic_check._check_var(
        dcolorbar, 'dcolorbar',
        default=defdcolorbar,
        types=dict,
    )

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

    # ----------
    # interactivity

    connect = ds._generic_check._check_var(
        connect, 'connect',
        types=bool,
        default=True,
    )

    return (
        key, keybs, keym, mtype,
        subbs, submtype,
        cmap, dcolorbar, dleg,
        connect,
    )


# #############################################################################
# #############################################################################
#                       Prepare
# #############################################################################


def _prepare(
    coll=None,
    key=None,
    keybs=None,
    keym=None,
    coefs=None,
    indt=None,
    dres=None,
    mtype=None,
    # submesh
    subbs=None,
    submtype=None,
    # levels
    levels=None,
):

    # ------------
    # misc

    # deg and
    wm = coll._which_mesh
    wbs = coll._which_bsplines
    deg = coll.dobj[wbs][subbs]['deg']
    if deg == 0:
        interp = 'nearest'
    elif deg == 1:
        interp = 'bilinear'
    elif deg >= 2:
        interp = 'bicubic'

    # -----------------
    # get plotting mesh

    coll2, dbs = coll.interpolate_all_bsplines(
        key=key,
        dres=dres,
        submesh=True,
    )

    bs2d = [k0 for k0, v0 in dbs.items() if len(v0['ref']) == 2][0]
    rX, rY = dbs[bs2d]['ref']
    lr1d = [k0 for k0 in coll2.ddata[key]['ref'] if k0 not in [rX, rY]]
    ndim = coll2.ddata[key]['data'].ndim

    dkeys = {
        'key': key,
        'keyX': coll2.get_ref_vector(ref=rX)[3],
        'keyY': coll2.get_ref_vector(ref=rY)[3],
        'keyZ': None,
        'keyU': None,
    }

    if ndim >= 3:
        dkeys['keyZ'] = coll2.get_ref_vector(ref=lr1d[0])[3]
        # uniform = ds._plot_as_array._check_uniform_lin(
            # k0=keyZ, ddata=coll2.ddata,
        # )
        # if not uniform:
            # keyZ = None
        if ndim == 4:
            dkeys['keyU'] = coll2.get_ref_vector(ref=lr1d[1])[3]

    # -----------------
    # optional contours

    key_cont = None
    axis = None
    dout, dref = None, None
    refZ, refU = None, None
    if levels is not None:

        # get contours
        key_cont = ['cont0', 'cont1']
        dout, dref = coll.get_profile2d_contours(
            key=key,
            levels=levels,
            res=dres if isinstance(dres, (int, float)) else None,
            store=False,
            return_dref=True,
            key_cont0=key_cont[0],
            key_cont1=key_cont[1],
        )

        # axis
        ref = dout['cont0']['ref']
        axis = [
            ii for ii, rr in enumerate(ref)
            if rr not in coll.dref.keys()
        ]

        # refZ, refU
        if len(ref) == 3:
            refZ = [
                rr for ii, rr in enumerate(ref)
                if ii not in axis
            ][0]

    return (
        coll2, dkeys, interp,
        key_cont, axis, dout, dref, refZ, refU,
    )


# #############################################################################
# #############################################################################
#                       Utilities
# #############################################################################


def _plot_bsplines_get_dx01(coll=None, km=None):
    # Get minimum distances

    wm = coll._which_mesh
    mtype = coll.dobj[wm][km]['type']
    if mtype == 'rect':
        knots0, knots1 = coll.dobj[wm][km]['knots']
        knots0 = coll.ddata[knots0]['data']
        knots1 = coll.ddata[knots1]['data']
        dx0 = np.min(np.diff(knots0))
        dx1 = np.min(np.diff(knots1))

    elif mtype == 'tri':
        indtri = coll.ddata[coll.dobj['mesh'][km]['ind']]['data']
        kknots = coll.dobj['mesh'][km]['knots']
        knots0 = coll.ddata[kknots[0]]['data']
        knots1 = coll.ddata[kknots[1]]['data']
        x0 = knots0[indtri]
        x1 = knots1[indtri]
        dist = np.mean(np.array([
            np.sqrt((x0[:, 1] - x0[:, 0])**2 + (x1[:, 1] - x1[:, 0])**2),
            np.sqrt((x0[:, 2] - x0[:, 1])**2 + (x1[:, 2] - x1[:, 1])**2),
            np.sqrt((x0[:, 2] - x0[:, 0])**2 + (x1[:, 2] - x1[:, 0])**2),
        ]))
        dx0, dx1 = dist, dist

    x0minmax = [knots0.min(), knots0.max()]
    x1minmax = [knots1.min(), knots1.max()]
    return dx0, dx1, x0minmax, x1minmax


# #############################################################################
# #############################################################################
#                       Back-up
# #############################################################################


# def _plot_profile2d_polar_add_radial(
    # coll=None,
    # key=None,
    # keym=None,
    # keybs=None,
    # dax=None,
# ):

    # # key to radius
    # kr2d = coll.dobj[coll._which_mesh][keym]['radius2d']
    # kr = coll.dobj[coll._which_mesh][keym]['knots'][0]
    # rr = coll.ddata[kr]['data']
    # rad = np.linspace(rr[0], rr[-1], rr.size*20)

    # # get angle if any
    # clas = coll.dobj['bsplines'][keybs]['class']
    # if clas.knotsa is None:
        # angle = None
    # elif len(clas.shapebs) == 2:
        # ka = coll.dobj['bsplines'][keybs]['apex'][1]
        # angle = coll.ddata[ka]['data']
    # elif np.sum(clas.nbs_a_per_r > 1) == 1:
        # i0 = (clas.nbs_a_per_r > 1).nonzero()[0][0]
        # angle = coll.dobj['bsplines'][keybs]['class'].apex_per_bs_a[i0]
    # else:
        # pass

    # if angle is None:
        # radmap = rad
        # anglemap = angle
    # else:
        # radmap = np.repeat(rad[:, None], angle.size, axis=1)
        # anglemap = np.repeat(angle[None, :], rad.size, axis=0)

    # # reft
    # reft, keyt, _, dind = coll.get_time_common(keys=[key, kr2d])[1:]

    # # radial total profile
    # radial, t_radial, _ = coll.interpolate_profile2d(
        # key=key,
        # radius=radmap,
        # angle=anglemap,
        # grid=False,
        # t=keyt,
    # )

    # if reft is not None and radial.ndim == radmap.ndim:
        # radial = np.repeat(radial[None, ...], t_radial.size, axis=0)

    # # details for purely-radial cases
    # if clas.knotsa is None:
        # radial_details, t_radial, _ = coll.interpolate_profile2d(
            # key=keybs,
            # radius=rad,
            # angle=None,
            # grid=False,
            # details=True,
        # )

        # if reft is None:
            # radial_details = radial_details * coll.ddata[key]['data'][None, :]
            # refdet = ('nradius',)
        # else:
            # refdet = (reft, 'nradius')
            # if reft == coll.get_time(key)[2]:
                # radial_details = (
                    # radial_details[None, :, :]
                    # * coll.ddata[key]['data'][:, None, :]
                # )
            # elif key in dind.keys():
                # radial_details = (
                    # radial_details[None, :, :]
                    # * coll.ddata[key]['data'][dind[key]['ind'], None, :]
                # )

        # nbs = radial_details.shape[-1]

    # # add to dax
    # dax.add_ref(key='nradius', size=rad.size)
    # if angle is not None:
        # dax.add_ref(key='nangle', size=angle.size)

    # if reft is not None:
        # assert radial.ndim > 1 and radial.shape[0] > 1
        # # if angle is not None:
            # # ref = (reft, 'nangle', 'nradius')
        # # else:
        # ref = (reft, 'nradius')
    # else:
        # # if angle is not None:
            # # ref = ('nangle', 'nradius')
        # # else:
        # ref = 'nradius'

    # # add to ddata
    # kradius = 'radius'
    # dax.add_data(key=kradius, data=rad, ref='nradius')
    # if angle is None:
        # lk = ['radial']
        # dax.add_data(key=lk[0], data=radial, ref=ref)
        # lkdet = [f'radial-detail-{ii}' for ii in range(nbs)]
        # for ii in range(nbs):
            # dax.add_data(
                # key=lkdet[ii], data=radial_details[..., ii], ref=refdet,
            # )

    # else:
        # kangle = 'angle'
        # dax.add_data(key=kangle, data=angle, ref='nangle')
        # lkdet = None
        # lk = [f'radial-{ii}' for ii in range(angle.size)]
        # if reft is None:
            # for ii in range(angle.size):
                # dax.add_data(key=lk[ii], data=radial[:, ii], ref=ref)
        # else:
            # for ii in range(angle.size):
                # dax.add_data(key=lk[ii], data=radial[:, :, ii], ref=ref)

    # return kradius, lk, lkdet, reft





# def _plot_profile2d_polar_create_axes(
    # fs=None,
    # dmargin=None,
# ):

    # if fs is None:
        # fs = (15, 9)

    # if dmargin is None:
        # dmargin = {
            # 'left': 0.05, 'right': 0.95,
            # 'bottom': 0.05, 'top': 0.95,
            # 'hspace': 0.4, 'wspace': 0.3,
        # }

    # fig = plt.figure(figsize=fs)
    # gs = gridspec.GridSpec(ncols=6, nrows=6, **dmargin)

    # # axes for image
    # ax0 = fig.add_subplot(gs[:4, 2:4], aspect='auto')

    # # axes for vertical profile
    # ax1 = fig.add_subplot(gs[:4, 4], sharey=ax0)

    # # axes for horizontal profile
    # ax2 = fig.add_subplot(gs[4:, 2:4], sharex=ax0)

    # # axes for traces
    # ax3 = fig.add_subplot(gs[2:4, :2])

    # # axes for traces
    # ax7 = fig.add_subplot(gs[:2, :2], sharey=ax2)

    # # axes for text
    # ax4 = fig.add_subplot(gs[:3, 5], frameon=False)
    # ax5 = fig.add_subplot(gs[3:, 5], frameon=False)
    # ax6 = fig.add_subplot(gs[4:, :2], frameon=False)

    # # dax
    # dax = {
        # # data
        # 'matrix': {'handle': ax0, 'type': 'matrix'},
        # 'vertical': {'handle': ax1, 'type': 'misc'},
        # 'horizontal': {'handle': ax2, 'type': 'misc'},
        # 'traces': {'handle': ax3, 'type': 'misc'},
        # 'radial': {'handle': ax7, 'type': 'misc'},
        # # text
        # 'textX': {'handle': ax4, 'type': 'text'},
        # 'textY': {'handle': ax5, 'type': 'text'},
        # 'textZ': {'handle': ax6, 'type': 'text'},
    # }
    # return dax
