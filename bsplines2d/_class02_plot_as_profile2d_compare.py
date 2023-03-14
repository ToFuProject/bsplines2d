# -*- coding: utf-8 -*-


# Common
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datastock as ds


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
):

    # --------------
    # check input

    if connect is None:
        connect = True

    (
        dkeys,
        dlevels,
        cmap, cmap_err, dcolorbar, dleg,
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

    lk0 = ['prof0', 'vert', 'hor', 'traces', 'spectrum', 'radial']
    lk1 = ['prof1', 'vert', 'hor', 'traces', 'spectrum', 'radial']
    dax2 = {
        keys[0]: {k0: dax[k0] for k0 in lk0 if k0 in dax.keys()},
        keys[1]: {k0: dax[k0] for k0 in lk1 if k0 in dax.keys()},
    }

    # ---------------
    # plot profiles2d

    coll2, dgroup = coll.plot_as_profile2d(
        key=keys,
        dres=dres,
        dlevels=dlevels,
        ref_com=ref_com,
        # details
        plot_details=plot_details,
        # plotting
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        # figure
        dax=dax2,
        dcolorbar=None,
        # interactivity
        connect=False,
    )

    # ----------
    # plot error

    # -----------
    # connect

    if connect is True:
        coll2.setup_interactivity(kinter='inter0', dgroup=dgroup, dinc=dinc)
        coll2.disconnect_old()
        coll2.connect()

        coll2.show_commands()
        return coll2
    else:
        return coll2, dgroup


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