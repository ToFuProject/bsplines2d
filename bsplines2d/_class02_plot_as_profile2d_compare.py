# -*- coding: utf-8 -*-


# Common
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datastock as ds


from . import _class02_interpolate_all as _interpolate_all
from . import _class02_plot_as_profile2d as _plot_as_profile2d


# ##############################################################
# ##############################################################
#                           Main
# ##############################################################


def plot_as_profile2d_compare(
    # ressources
    coll=None,
    # inputs
    keys=None,
    # parameters
    dres=None,
    # contours
    dlevels=None,
    ref_com=None,
    # details
    plot_details=None,
    # ref vectors
    dref_vectorZ=None,
    dref_vectorU=None,
    ref_vector_strategy=None,
    uniform=None,
    # interpolation
    val_out=None,
    nan0=None,
    # plotting
    vmin=None,
    vmax=None,
    cmap=None,
    vmin_err=None,
    vmax_err=None,
    cmap_err=None,
    # figure
    dax=None,
    dmargin=None,
    fs=None,
    dcolorbar=None,
    dleg=None,
    # interactivity
    dinc=None,
    connect=None,
    show_commands=None,
):

    # --------------
    # check input

    if connect is None:
        connect = True

    (
        dkeys,
        dlevels,
        cmap, cmap_err,
        dvminmax,
        dcolorbar, dleg,
        connect,
    ) = _plot_as_profile2d._check(
        coll=coll,
        keys=keys,
        dlevels=dlevels,
        # plotting
        cmap=cmap,
        cmap_err=cmap_err,
        # figure
        dcolorbar=dcolorbar,
        dleg=dleg,
        # interactivity
        connect=connect,
    )

    # ---------------
    # prepare dax

    if dax is None:
        dax = _create_axes(
            dkeys=dkeys,
            fs=fs,
            dmargin=dmargin,
        )

    # safety check
    if keys[0] == keys[1]:
        msg = (
            "Please provide 2 different data keys!\n"
            f"Provided: {keys}"
        )
        raise Exception(msg)

    # prepare dax2
    lk0 = ['prof0', 'vert', 'hor', 'traces', 'spectrum', 'radial']
    lk1 = ['prof1', 'vert', 'hor', 'traces', 'spectrum', 'radial']
    dax2 = {
        keys[0]: {k0: dax[k0] for k0 in lk0 if k0 in dax.keys()},
        keys[1]: {k0: dax[k0] for k0 in lk1 if k0 in dax.keys()},
    }

    # ----------------------
    # get unique 2d sampling

    keys, _, dres, submesh = _interpolate_all._check(
        coll=coll,
        keys=keys,
        # sampling
        dres=dres,
        submesh=True,
    )

    # list meshes 2d
    wm = coll._which_mesh
    lm = list(dres.keys())
    lm = [k0 for k0 in lm if coll.dobj[wm][k0]['nd'] == '2d']

    if dres[lm[0]]['x0'] is None:
        dunique = coll.get_sample_mesh(
            key=lm[0],
            res=dres[lm[0]]['res'],
            mode=dres[lm[0]]['mode'],
            grid=False,
            # store
            store=False,
            kx0='unique_x0',
            kx1='unique_x1',
        )
    else:
        kx01 = coll.dobj[wm][lm[0]]['knots']
        dunique = {
            'x0': {
                'key': 'x0_temp',
                'data':  dres[lm[0]]['x0'],
                'units': coll.ddata[kx01[0]]['units'],
                'name': coll.ddata[kx01[0]]['name'],
                'quant': coll.ddata[kx01[0]]['quant'],
                'dim': coll.ddata[kx01[0]]['dim'],
            },
            'x1': {
                'key': 'x1_temp',
                'data':  dres[lm[0]]['x1'],
                'units': coll.ddata[kx01[1]]['units'],
                'name': coll.ddata[kx01[1]]['name'],
                'quant': coll.ddata[kx01[1]]['quant'],
                'dim': coll.ddata[kx01[1]]['dim'],
            },
        }

    # ---------------
    # plot profiles2d

    collax, dgroup = coll.plot_as_profile2d(
        key=keys,
        dres=dres,
        dunique_mesh_2d=dunique,
        dlevels=dlevels,
        ref_com=ref_com,
        # details
        plot_details=plot_details,
        # ref vectors
        dref_vectorZ=dref_vectorZ,
        dref_vectorU=dref_vectorU,
        ref_vector_strategy=ref_vector_strategy,
        uniform=uniform,
        # interpolation
        val_out=val_out,
        nan0=nan0,
        # plotting
        dvminmax=dvminmax,
        cmap=cmap,
        # figure
        dax=dax2,
        dcolorbar=None,
        # interactivity
        connect=False,
    )

    # ----------
    # get error
    # ----------

    # data and vmin, vmax
    err = collax.ddata[keys[0]]['data'] - collax.ddata[keys[1]]['data']
    vmax = np.nanmax(np.abs(err))

    # deg and interp
    k0 = f"{keys[0]}_im"
    interp = collax.dobj['mobile'][k0]['handle'].get_interpolation()
    keyX = collax.dobj['axes']['prof0']['datax'][0]
    keyY = collax.dobj['axes']['prof0']['datay'][0]

    # add to collax
    collax.add_data(
        key='err',
        data=err,
        ref=collax.ddata[keys[0]]['ref'],
        units=collax.ddata[keys[0]]['units'],
    )

    # dax
    laxerr = ['err', 'vert_err', 'hor_err']

    # plot
    collax, dgroup2 = collax.plot_as_array(
        # keys
        key='err',
        keyX=keyX,
        keyY=keyY,
        uniform=uniform,
        # options
        dvminmax={'data': {'min': -vmax, 'max': vmax}},
        # vmin=-vmax,
        # vmax=vmax,
        cmap=plt.cm.seismic,
        dax={kax: dax[kax] for kax in laxerr},
        dmargin=dmargin,
        fs=fs,
        dcolorbar=dcolorbar,
        dleg=dleg,
        interp=interp,
        color_dict=None,
        # nmax=nmax,
        connect=False,
        inplace=True,
    )

    # -----------
    # connect

    if connect is True:
        collax.setup_interactivity(kinter='inter0', dgroup=dgroup, dinc=dinc)
        collax.disconnect_old()
        collax.connect()

        collax.show_commands(verb=show_commands)
        return collax
    else:
        return collax, dgroup


# ################################################################
# ################################################################
#                       Create axes
# ################################################################


def _create_axes(
    dkeys=None,
    fs=None,
    dmargin=None,
):

    # ------------------
    # check and prepare

    ndim_extra = len(list(dkeys.values())[0]['ref_other'])
    hassubmesh = any([v0['submesh'] is not None for v0 in dkeys.values()])
    if hassubmesh:
        ndim_extra += 1

    if fs is None:
        fs = (15, 9)

    if dmargin is None:
        dmargin = {
            'left': 0.04, 'right': 0.98,
            'bottom': 0.06, 'top': 0.92,
            'hspace': 0.4, 'wspace': 2.,
        }

    # ----------------
    # create axis grid

    # axes for images
    dgs = {}
    nrows = 3
    if ndim_extra == 0:
        ncols = 17
    else:
        ncols = 22

    gs = gridspec.GridSpec(ncols=ncols, nrows=nrows, **dmargin)

    dgs['prof0'] = gs[:2, -17:-13]
    dgs['prof1'] = gs[:2, -13:-9]
    dgs['vert'] = gs[:2, -9:-7]
    dgs['err'] = gs[:2, -6:-2]
    dgs['vert_err'] = gs[:2, -2:]
    dgs['hor'] = gs[2, -17:-13]
    dgs['hor_err'] = gs[2, -6:-2]

    for ii in range(0, ndim_extra):
        if ii == 0 and hassubmesh:
            kk = 'radial'
        else:
            kk = ['traces', 'spectrum'][ii - hassubmesh]
        dgs[kk] = gs[ii, :4]

    # ----------------
    # prepare figure and dax

    fig = plt.figure(figsize=fs)

    dax = {}

    # ------------
    # 2d profiles

    kax = 'prof0'
    if dgs.get(kax) is not None:
        ax = fig.add_subplot(dgs[kax])
        dax[kax] = {'handle': ax, 'type': 'matrix'}

    kax = 'prof1'
    if dgs.get(kax) is not None:
        ax = fig.add_subplot(
            dgs[kax],
            sharex=dax['prof0']['handle'],
            sharey=dax['prof0']['handle'],
        )
        dax[kax] = {'handle': ax, 'type': 'matrix'}

    kax = 'err'
    if dgs.get(kax) is not None:
        ax = fig.add_subplot(
            dgs[kax],
            sharex=dax['prof0']['handle'],
            sharey=dax['prof0']['handle'],
        )
        dax[kax] = {'handle': ax, 'type': 'matrix'}

        ax.set_title('difference', size=12, fontweight='bold')

    # --------------------
    # hor and vert slices

    kax = 'hor'
    if dgs.get(kax) is not None:
        ax = fig.add_subplot(
            dgs[kax],
            sharex=dax['prof0']['handle'],
        )
        dax[kax] = {'handle': ax, 'type': 'horizontal'}

    kax = 'hor_err'
    if dgs.get(kax) is not None:
        ax = fig.add_subplot(
            dgs[kax],
            sharex=dax['prof0']['handle'],
        )
        ax.axhline(0., c='k', lw=1., ls='--')
        dax[kax] = {'handle': ax, 'type': 'horizontal'}

    kax = 'vert'
    if dgs.get(kax) is not None:
        ax = fig.add_subplot(
            dgs[kax],
            sharey=dax['prof0']['handle'],
        )
        dax[kax] = {'handle': ax, 'type': 'vertical'}

    kax = 'vert_err'
    if dgs.get(kax) is not None:
        ax = fig.add_subplot(
            dgs[kax],
            sharey=dax['prof0']['handle'],
        )
        ax.axvline(0., c='k', lw=1., ls='--')
        dax[kax] = {'handle': ax, 'type': 'vertical'}

    # --------------------
    # extra dimensions

    kax = 'traces'
    if dgs.get(kax) is not None:
        ax = fig.add_subplot(
            dgs[kax],
            sharey=dax['hor']['handle'],
        )
        dax[kax] = {'handle': ax, 'type': 'tracesZ'}

    kax = 'spectrum'
    if dgs.get(kax) is not None:
        ax = fig.add_subplot(
            dgs[kax],
            sharey=dax['hor']['handle'],
        )
        dax[kax] = {'handle': ax, 'type': 'tracesU'}

    kax = 'radial'
    if dgs.get(kax) is not None:
        ax = fig.add_subplot(
            dgs[kax],
            sharey=dax['hor']['handle'],
        )
        dax[kax] = {'handle': ax}

    kax = 'radial_err'
    if dgs.get(kax) is not None:
        ax = fig.add_subplot(
            dgs[kax],
            sharex=dax['radial']['handle'],
            sharey=dax['hor']['handle'],
        )
        dax[kax] = {'handle': ax}

    return dax