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


def main(
    coll=None,
    key=None,
    # options for 2d mesh
    ind_knot=None,
    ind_cent=None,
    crop=None,
    bck=None,
    nmax=None,
    # plotting options
    color=None,
    dax=None,
    dmargin=None,
    fs=None,
    dleg=None,
):
    """ Plot the desired 3d mesh

    rect and tri meshes are constant
    polar meshes can vary in time

    """

    # --------------
    # check inputs
    # --------------

    (
        key, key_mesh2d, key_mesh1d_angle,
        color, dleg,
     ) = _check(
        coll=coll,
        key=key,
        color=color,
        dleg=dleg,
    )

    # ------------
    # plot
    # ------------

    dax = _plot(
        coll=coll,
        key=key,
        key_mesh2d=key_mesh2d,
        key_mesh1d_angle=key_mesh1d_angle,
        # plotting
        color=color,
        dleg=dleg,
    )

    return dax

# ###############################################
# ###############################################
#                   checks inputs
# ###############################################


def _check(
    coll=None,
    key=None,
    ind_knot=None,
    ind_cent=None,
    crop=None,
    bck=None,
    color=None,
    dleg=None,
):

    # ----------
    # key
    # ----------

    wm3d = coll._which_mesh3d
    lok = list(coll.dobj.get(wm3d, {}).keys())
    key = ds._generic_check._check_var(
        key, 'key',
        types=str,
        allowed=lok,
    )

    key_mesh2d = coll.dobj[wm3d][key]['mesh2d']
    key_mesh1d_angle = coll.dobj[wm3d][key]['mesh1d_angle']

    # ----------
    # color
    # ----------

    # color
    if color is None:
        color = 'k'
    if not mcolors.is_color_like(color):
        msg = (
            "Arg color must be a valid matplotlib color identifier!\n"
            f"Provided: {color}"
        )
        raise Exception(msg)

    # ----------
    # dleg
    # ----------

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
        key, key_mesh2d, key_mesh1d_angle,
        color, dleg,
    )


# ###########################################
# ###########################################
#                   prepare
# ###########################################


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
