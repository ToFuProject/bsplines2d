

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


from .._class03_Bins import Bins as Collection


# #######################################################
# #######################################################
#              Main
# #######################################################


def main():

    # ---------------
    # Instanciate
    # ---------------

    coll = Collection()

    # -----------------------
    # Add meshes and bsplines
    # -----------------------

    # rectangular
    _add_mesh_rect(coll)

    # triangular
    _add_mesh_tri(coll)

    # -----------------------
    # plot
    # -----------------------

    dax = _plot(coll)

    return coll, dax


# #######################################################
# #######################################################
#              Add mesh / bsplines
# #######################################################


# ##############
#   Rect
# ##############

def _add_mesh_rect(coll):

    # ----------
    # mesh
    # ----------

    X = np.linspace(0, 5, 6)
    Y = np.linspace(0, 5, 6)

    km = 'rect'
    coll.add_mesh_2d_rect(
        key=km,
        knots0=X,
        knots1=Y,
    )

    # ----------
    # bsplines
    # ----------

    for deg in [0, 1, 2]:
        coll.add_bsplines(key=km, deg=deg)

    return


# ##############
#   Triangular
# ##############


def _add_mesh_tri(coll):

    # ----------
    # mesh
    # ----------

    X = np.linspace(0.5, 5.5, 9)
    Y = np.linspace(0, 5, 6)

    Xf = np.r_[X[1::2], X[0::2], X[1::2], X[0::2], X[1::2], X[0::2]]
    Yf = np.repeat(Y, np.r_[4, 5, 4, 5, 4, 5])
    knots = np.array([Xf, Yf]).T

    ind = np.array([
        [4, 0, 5],
        [0, 1, 5],
        [5, 1, 6],
        [1, 2, 6],
        [6, 2, 7],
        [2, 3, 7],
        [7, 3, 8],
    ])

    iused = np.unique(ind)

    km = 'mtri'
    coll.add_mesh_2d_tri(
        key=km,
        knots=knots[iused, :],
        indices=ind,
    )

    # ----------
    # bsplines
    # ----------

    for deg in [0, 1]:
        coll.add_bsplines(key=km, deg=deg)

    return


# #######################################################
# #######################################################
#              Plot
# #######################################################


def _plot(coll):

    # -------------
    # prepare keys
    # -------------

    wm = coll._which_mesh
    wbs = coll._which_bsplines
    lkm = sorted(coll.dobj[wm].keys())
    lkbs = sorted(coll.dobj[wbs].keys())
    ldeg = sorted(set([coll.dobj[wbs][kbs]['deg'] for kbs in lkbs]))

    X = coll.ddata[coll.dobj[wm]['rect']['knots'][0]]['data']

    # -------------
    # prepare data
    # -------------

    nn = X.size * 21
    xx = np.linspace(X[0], X[-1], nn)
    xxf = np.repeat(xx[:, None], nn, axis=1)
    yyf = np.repeat(xx[None, :], nn, axis=0)

    # add data
    coll.add_data('xxf', data=xxf, ref=('nx', 'ny'), units='m')
    coll.add_data('yyf', data=yyf, ref=('nx', 'ny'), units='m')

    for kbs in lkbs:

        _ = coll.interpolate(
            keys=None,
            ref_key=kbs,
            x0='xxf',
            x1='yyf',
            grid=False,
            store=True,
            store_keys=f"{kbs}_data",
        )

    # -------------
    # figure
    # -------------

    dmargin = {
        'bottom': 0.08, 'top': 0.9,
        'left': 0.05, 'right': 0.95,
        'wspace': 0.15, 'hspace': 0.15,
    }

    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(len(lkm), len(ldeg), figure=fig, **dmargin)

    # -------------
    # axes
    # -------------

    dax = {}
    ax0 = None
    for im, mtype in enumerate(lkm):
        for ideg, deg in enumerate(ldeg):

            # axes creation
            ax = fig.add_subplot(
                gs[im, ideg],
                aspect='equal',
                sharex=ax0,
                sharey=ax0,
            )

            # labels
            ax.set_xlabel('X (m)', size=12, fontweight='bold')
            ax.set_ylabel('Y (m)', size=12, fontweight='bold')

            # title
            if im == 0:
                ax.set_title(f'deg = {deg}', size=14, fontweight='bold')

            # set ax0
            if ax0 is None:
                ax0 = ax

            # store
            kax = f'{mtype}_{deg}'
            dax[kax] = ax

    # -------------
    # plot meshes
    # -------------

    for im, mtype in enumerate(lkm):
        for ideg, deg in enumerate(ldeg):
            kax = f'{mtype}_{deg}'
            ax = dax[kax]

            coll.plot_mesh(
                key=mtype,
                dax={'cross': {'handle': ax}},
                dleg=False,
            )

    ax0.set_xlim(X[0], X[-1])

    return dax
