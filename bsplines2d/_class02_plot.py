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

    # ---------------
    # prepare figure
    # --------------

    if dax is None:
        dax = _prepare_figure(
            # keys
            key=key,
            key_mesh2d=key_mesh2d,
            key_mesh1d_angle=key_mesh1d_angle,
            # options
            dmargin=dmargin,
            fs=fs,
        )

    dax = _generic_check._check_dax(dax=dax, main='cross')

    # ---------------
    # plot 2d mesh
    # --------------

    kax = 'cross'
    if dax.get(kax) is not None:
        coll.plot_mesh(
            key_mesh2d,
            crop=crop,
            color=color,
            dax=dax,
            dleg=dleg,
        )

    # ---------------
    # plot angle mesh
    # --------------

    kax = 'hor'
    if dax.get(kax) is not None:
        _plot_mesh1d_angle(
            coll=coll,
            key=key,
            key_mesh2d=key_mesh2d,
            key_mesh1d_angle=key_mesh1d_angle,
            # options
            color=color,
            dax=dax,
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
    if dleg is True:
        defdleg = {
            'bbox_to_anchor': (1.1, 1.),
            'loc': 'upper left',
            'frameon': True,
        }
        dleg = None
    else:
        defdleg = False

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


def _plot_mesh1d_angle(
    coll=None,
    key=None,
    key_mesh2d=None,
    key_mesh1d_angle=None,
    # plotting
    color=None,
    dax=None,
):

    # --------
    # prepare
    # --------

    wm = coll._which_mesh
    phi = np.pi * np.linspace(-1, 1, 101)
    cos = np.cos(phi)
    sin = np.sin(phi)

    # -----------
    # mesh 2d
    # -----------

    knots0 = coll.dobj[wm][key_mesh2d]['knots'][0]
    knots0 = coll.ddata[knots0]['data']
    knots0f = np.repeat(knots0[None, :], phi.size, axis=0)
    knots0f_cos = knots0f * cos[:, None]
    knots0f_sin = knots0f * sin[:, None]

    # edges
    xx_edges = np.r_[knots0f_cos[:, 0], np.nan, knots0f_cos[:, -1]]
    yy_edges = np.r_[knots0f_sin[:, 0], np.nan, knots0f_sin[:, -1]]

    # inside
    nan = np.full(knots0f[0:1, 1:-1].shape, np.nan)
    xx_in = np.concatenate((knots0f_cos[:, 1:-1], nan), axis=0).T.ravel()
    yy_in = np.concatenate((knots0f_sin[:, 1:-1], nan), axis=0).T.ravel()

    # -----------
    # mesh 1d angles
    # -----------

    # angles
    knots = coll.dobj[wm][key_mesh1d_angle]['knots'][0]
    knots = coll.ddata[knots]['data']

    R = np.r_[knots0[0], knots0[-1], np.nan]
    xx_angles = (R[None, :] * np.cos(knots)[:, None]).ravel()
    yy_angles = (R[None, :] * np.sin(knots)[:, None]).ravel()

    # --------
    # plot
    # --------

    ax = dax['hor']['handle']

    # inner and outer edges
    ax.plot(
        xx_edges,
        yy_edges,
        ls='-',
        lw=1.5,
        c=color,
    )

    # inside
    ax.plot(
        xx_in,
        yy_in,
        ls='-',
        lw=1.,
        c=mcolors.to_rgb(color) + (0.5,),
    )

    # angles
    ax.plot(
        xx_angles,
        yy_angles,
        ls='-',
        lw=1.5,
        c=color,
    )

    return


# ################################################
# ################################################
#               Prepare figure
# ################################################


def _prepare_figure(
    # keys
    key=None,
    key_mesh2d=None,
    key_mesh1d_angle=None,
    # options
    dmargin=None,
    fs=None,
):

    # ----------
    # figure
    # ----------

    if fs is None:
        fs = (12, 8)

    if dmargin is None:
        dmargin = {
            'left': 0.08, 'right': 0.95,
            'bottom': 0.08, 'top': 0.9,
            'hspace': 0.1, 'wspace': 0.15,
        }

    fig = plt.figure(figsize=fs)
    gs = gridspec.GridSpec(ncols=2, nrows=1, **dmargin)
    dax = {}

    # ----------
    # axes
    # ----------

    # cross
    ax = fig.add_subplot(gs[:, 0], aspect='equal')
    ax.set_xlabel('R (m)', size=12, fontweight='bold')
    ax.set_ylabel('Z (m)', size=12, fontweight='bold')
    tit = f"{key} - {key_mesh2d}"
    ax.set_title(tit, size=12, fontweight='bold')

    dax['cross'] = {'handle': ax}

    # hor
    ax = fig.add_subplot(gs[:, 1], aspect='equal')
    ax.set_xlabel('X (m)', size=12, fontweight='bold')
    ax.set_ylabel('Y (m)', size=12, fontweight='bold')
    tit = f"{key} - {key_mesh1d_angle}"
    ax.set_title(tit, size=12, fontweight='bold')

    dax['hor'] = {'handle': ax}

    return dax
