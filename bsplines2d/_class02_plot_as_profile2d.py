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
from . import _class01_compute as _compute


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
    res=None,
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
        coll2, dkeys, interp, lcol,
    ) = _prepare(
        coll=coll,
        key=key,
        keybs=keybs,
        subbs=subbs,
        keym=keym,
        res=res,
        mtype=mtype,
        # submesh
        subbs=subbs,
        submtype=submtype,
    )

    # ---------------
    # call right function

    if mtype in ['rect', 'tri']:
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

        if dax is None:
            dax = _plot_profile2d_polar_create_axes(
                fs=fs,
                dmargin=dmargin,
            )

        dax, dgroup = coll2.plot_as_array(
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            dax=dax,
            dmargin=dmargin,
            fs=fs,
            dcolorbar=dcolorbar,
            dleg=dleg,
            connect=False,
            interp=interp,
            label=True,
            **dkeys,
        )

        # ------------------
        # add radial profile to dax

        kradius, lkradial, lkdet, reft = _plot_profile2d_polar_add_radial(
            coll=coll,
            key=key,
            keym=keym,
            keybs=keybs,
            dax=dax,
        )


        assert (reft is not None) == ('Z' in dgroup.keys())
        if reft is not None and reft not in dgroup['Z']['ref']:
            dgroup['Z']['ref'].append(reft)
            dgroup['Z']['data'].append('index')

        # ------------------
        # add radial profile

        kax = 'radial'
        if dax.dax.get(kax) is not None:
            ax = dax.dax[kax]['handle']
            for ii in range(len(lkradial)):

                if reft is None:
                    l0, = ax.plot(
                        dax.ddata[kradius]['data'],
                        dax.ddata[lkradial[ii]]['data'],
                        c=lcol[ii],
                        ls='-',
                        lw=2,
                    )
                else:
                    l0, = ax.plot(
                        dax.ddata[kradius]['data'],
                        dax.ddata[lkradial[ii]]['data'][0, :],
                        c=lcol[ii],
                        ls='-',
                        lw=2,
                    )

                    kl = f"radial{ii}"
                    dax.add_mobile(
                        key=kl,
                        handle=l0,
                        refs=(reft,),
                        data=[lkradial[ii]],
                        dtype=['ydata'],
                        axes=kax,
                        ind=0,
                    )

            if lkdet is not None:
                for ii in range(len(lkdet)):
                    if reft is None:
                        l0, = ax.plot(
                            dax.ddata[kradius]['data'],
                            dax.ddata[lkdet[ii]]['data'],
                            ls='-',
                            lw=1,
                        )
                    else:
                        l0, = ax.plot(
                            dax.ddata[kradius]['data'],
                            dax.ddata[lkdet[ii]]['data'][0, :],
                            ls='-',
                            lw=1,
                        )

                        kl = f"radial_det{ii}"
                        dax.add_mobile(
                            key=kl,
                            handle=l0,
                            refs=(reft,),
                            data=[lkdet[ii]],
                            dtype=['ydata'],
                            axes=kax,
                            ind=0,
                        )

            ax.set_xlim(
                dax.ddata[kradius]['data'].min(),
                dax.ddata[kradius]['data'].max(),
            )

            if vmin is not None:
                ax.set_ylim(bottom=vmin)
            if vmax is not None:
                ax.set_ylim(top=vmax)

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
    res=None,
    mtype=None,
    # submesh
    subbs=None,
    submtype=None,
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

    # get dR, dZ
    dx0, dx1, x0minmax, x1minmax = _plot_bsplines_get_dx01(
        coll=coll,
        km=coll.dobj[wbs][subbs][wm],
    )

    if res is None:
        res_coef = 0.2
        res = [res_coef*dx0, res_coef*dx1]

    # compute
    coll2 = coll.interpolate_profile2d(
        key=key,
        res=res,
        details=False,
        # return vs store
        returnas=object,
        return_params=False,
        store=True,
        inplace=False,
    )


    keymap = [k0 for k0, v0 in coll2.ddata.items() if v0['data'].ndim > 1][0]
    ndim = coll2.ddata[keymap]['data'].ndim

    refmap = coll2.ddata[keymap]['ref']
    dkeys = {
        'key': keymap,
        'keyX': coll2.get_ref_vector(key=keymap, ref=refmap[-2])[3],
        'keyY': coll2.get_ref_vector(key=keymap, ref=refmap[-1])[3],
        'keyZ': None,
    }
    if ndim == 3:
        keyZ = coll2.get_ref_vector(key=keymap, ref=refmap[0])[3]
        import datastock as ds
        uniform = ds._plot_as_array._check_uniform_lin(
            k0=keyZ, ddata=coll2.ddata,
        )
        if not uniform:
            keyZ = None


    # radial of polar
    if mtype == 'polar':
        # lcol
        lcol = ['k', 'r', 'b', 'g', 'm', 'c', 'y']
    else:
        lcol = None

    return coll2, dkeys, interp, lcol


def _plot_bsplines_get_dx01(coll=None, km=None):
    # Get minimum distances

    wm = coll._which_mesh
    mtype = coll.dobj[wm][km]
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


def _plot_profile2d_polar_add_radial(
    coll=None,
    key=None,
    keym=None,
    keybs=None,
    dax=None,
):

    # key to radius
    kr2d = coll.dobj[coll._which_mesh][keym]['radius2d']
    kr = coll.dobj[coll._which_mesh][keym]['knots'][0]
    rr = coll.ddata[kr]['data']
    rad = np.linspace(rr[0], rr[-1], rr.size*20)

    # get angle if any
    clas = coll.dobj['bsplines'][keybs]['class']
    if clas.knotsa is None:
        angle = None
    elif len(clas.shapebs) == 2:
        ka = coll.dobj['bsplines'][keybs]['apex'][1]
        angle = coll.ddata[ka]['data']
    elif np.sum(clas.nbs_a_per_r > 1) == 1:
        i0 = (clas.nbs_a_per_r > 1).nonzero()[0][0]
        angle = coll.dobj['bsplines'][keybs]['class'].apex_per_bs_a[i0]
    else:
        pass

    if angle is None:
        radmap = rad
        anglemap = angle
    else:
        radmap = np.repeat(rad[:, None], angle.size, axis=1)
        anglemap = np.repeat(angle[None, :], rad.size, axis=0)

    # reft
    reft, keyt, _, dind = coll.get_time_common(keys=[key, kr2d])[1:]

    # radial total profile
    radial, t_radial, _ = coll.interpolate_profile2d(
        key=key,
        radius=radmap,
        angle=anglemap,
        grid=False,
        t=keyt,
    )

    if reft is not None and radial.ndim == radmap.ndim:
        radial = np.repeat(radial[None, ...], t_radial.size, axis=0)

    # details for purely-radial cases
    if clas.knotsa is None:
        radial_details, t_radial, _ = coll.interpolate_profile2d(
            key=keybs,
            radius=rad,
            angle=None,
            grid=False,
            details=True,
        )

        if reft is None:
            radial_details = radial_details * coll.ddata[key]['data'][None, :]
            refdet = ('nradius',)
        else:
            refdet = (reft, 'nradius')
            if reft == coll.get_time(key)[2]:
                radial_details = (
                    radial_details[None, :, :]
                    * coll.ddata[key]['data'][:, None, :]
                )
            elif key in dind.keys():
                radial_details = (
                    radial_details[None, :, :]
                    * coll.ddata[key]['data'][dind[key]['ind'], None, :]
                )

        nbs = radial_details.shape[-1]

    # add to dax
    dax.add_ref(key='nradius', size=rad.size)
    if angle is not None:
        dax.add_ref(key='nangle', size=angle.size)

    if reft is not None:
        assert radial.ndim > 1 and radial.shape[0] > 1
        # if angle is not None:
            # ref = (reft, 'nangle', 'nradius')
        # else:
        ref = (reft, 'nradius')
    else:
        # if angle is not None:
            # ref = ('nangle', 'nradius')
        # else:
        ref = 'nradius'

    # add to ddata
    kradius = 'radius'
    dax.add_data(key=kradius, data=rad, ref='nradius')
    if angle is None:
        lk = ['radial']
        dax.add_data(key=lk[0], data=radial, ref=ref)
        lkdet = [f'radial-detail-{ii}' for ii in range(nbs)]
        for ii in range(nbs):
            dax.add_data(
                key=lkdet[ii], data=radial_details[..., ii], ref=refdet,
            )

    else:
        kangle = 'angle'
        dax.add_data(key=kangle, data=angle, ref='nangle')
        lkdet = None
        lk = [f'radial-{ii}' for ii in range(angle.size)]
        if reft is None:
            for ii in range(angle.size):
                dax.add_data(key=lk[ii], data=radial[:, ii], ref=ref)
        else:
            for ii in range(angle.size):
                dax.add_data(key=lk[ii], data=radial[:, :, ii], ref=ref)

    return kradius, lk, lkdet, reft





def _plot_profile2d_polar_create_axes(
    fs=None,
    dmargin=None,
):

    if fs is None:
        fs = (15, 9)

    if dmargin is None:
        dmargin = {
            'left': 0.05, 'right': 0.95,
            'bottom': 0.05, 'top': 0.95,
            'hspace': 0.4, 'wspace': 0.3,
        }

    fig = plt.figure(figsize=fs)
    gs = gridspec.GridSpec(ncols=6, nrows=6, **dmargin)

    # axes for image
    ax0 = fig.add_subplot(gs[:4, 2:4], aspect='auto')

    # axes for vertical profile
    ax1 = fig.add_subplot(gs[:4, 4], sharey=ax0)

    # axes for horizontal profile
    ax2 = fig.add_subplot(gs[4:, 2:4], sharex=ax0)

    # axes for traces
    ax3 = fig.add_subplot(gs[2:4, :2])

    # axes for traces
    ax7 = fig.add_subplot(gs[:2, :2], sharey=ax2)

    # axes for text
    ax4 = fig.add_subplot(gs[:3, 5], frameon=False)
    ax5 = fig.add_subplot(gs[3:, 5], frameon=False)
    ax6 = fig.add_subplot(gs[4:, :2], frameon=False)

    # dax
    dax = {
        # data
        'matrix': {'handle': ax0, 'type': 'matrix'},
        'vertical': {'handle': ax1, 'type': 'misc'},
        'horizontal': {'handle': ax2, 'type': 'misc'},
        'traces': {'handle': ax3, 'type': 'misc'},
        'radial': {'handle': ax7, 'type': 'misc'},
        # text
        'textX': {'handle': ax4, 'type': 'text'},
        'textY': {'handle': ax5, 'type': 'text'},
        'textZ': {'handle': ax6, 'type': 'text'},
    }
    return dax
