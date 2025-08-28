

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import datastock as ds


from .._class03_Bins import Bins as Collection


# #######################################################
# #######################################################
#              Main
# #######################################################


def main(details=True):

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

    dax = _plot(coll, details=details)

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

    X = np.linspace(0, 5, 11)
    Y = np.linspace(0, 5, 6)

    Xf = np.r_[X[::2], X[1:-1:2], X[::2], X[1:-1:2], X[::2], X[1:-1:2]]
    Yf = np.repeat(Y, np.r_[6, 5, 6, 5, 6, 5])
    knots = np.array([Xf, Yf]).T

    ind = np.array([
        [0, 6, 1],
        [1, 7, 2],
        [2, 8, 3],
        [3, 9, 4],
        [4, 10, 5],
        [6, 1, 7],
        [7, 2, 8],
        [8, 3, 9],
        [9, 4, 10],
        [6, 12, 7],
        [7, 13, 8],
        [8, 14, 9],
        [9, 15, 10],
        [11, 6, 12],
        [12, 7, 13],
        [13, 8, 14],
        [14, 9, 15],
        [15, 10, 16],
        [11, 17, 12],
        [12, 18, 13],
        [13, 19, 14],
        [14, 20, 15],
        [15, 21, 16],
        [17, 12, 18],
        [18, 13, 19],
        [19, 14, 20],
        [20, 15, 21],
        [17, 23, 18],
        [18, 24, 19],
        [19, 25, 20],
        [20, 26, 21],
        [22, 17, 23],
        [23, 18, 24],
        [24, 19, 25],
        [25, 20, 26],
        [26, 21, 27],
        [22, 28, 23],
        [23, 29, 24],
        [24, 30, 25],
        [25, 31, 26],
        [26, 32, 27],
        [28, 23, 29],
        [29, 24, 30],
        [30, 25, 31],
        [31, 26, 32],
    ])

    iused = np.unique(ind)

    km = 'tri'
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


def _plot(coll, details=None):

    # -------------
    # check inputs
    # -------------

    details = ds._generic_check._check_var(
        details, 'details',
        types=bool,
        default=True,
    )

    # -------------
    # prepare keys
    # -------------

    wm = coll._which_mesh
    wbs = coll._which_bsplines
    lkm = sorted(coll.dobj[wm].keys())
    lkbs = sorted(coll.dobj[wbs].keys())
    ldeg = sorted(set([coll.dobj[wbs][kbs]['deg'] for kbs in lkbs]))

    X = coll.ddata[coll.dobj[wm]['rect']['knots'][0]]['data']
    Xm = np.mean(X)
    Xm = (Xm, Xm)

    # -------------
    # prepare data
    # -------------

    nn = X.size * 51
    xx = np.linspace(X[0], X[-1], nn)
    xxf = np.repeat(xx[:, None], nn, axis=1)
    yyf = np.repeat(xx[None, :], nn, axis=0)

    dx = 0.5 * (xx[1] - xx[0])
    extent = (xx[0] - dx, xx[-1] + dx, xx[0] - dx, xx[-1] + dx)

    # add data
    coll.add_data('xxf', data=xxf, ref=('nx', 'ny'), units='m')
    coll.add_data('yyf', data=yyf, ref=('nx', 'ny'), units='m')

    # add interpolate
    for kbs in lkbs:

        # shape = coll.dobj[wbs][kbs]['shape']
        # coefs = np.zeros(shape, dtype=float)

        # ind = np.unravel_index(int(np.prod(shape)/2), shape)
        # coefs[ind] = 1.

        # single bs coef
        # kcoef = f"{kbs}_coef"
        # coll.add_data(
            # kcoef,
            # data=data,
            # ref=kbs,
            # units=None,
        # )

        # add interpolate
        _ = coll.interpolate(
            keys=None,
            ref_key=kbs,
            x0='xxf',
            x1='yyf',
            nan0=True,
            grid=False,
            details=True,
            inplace=True,
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

    # -------------
    # plot data
    # -------------

    for im, mtype in enumerate(lkm):
        for ideg, deg in enumerate(ldeg):
            kax = f'{mtype}_{deg}'
            ax = dax[kax]

            # prepare
            kbs = f"{mtype}_bs{deg}"
            # check
            if kbs not in coll.dobj[wbs].keys():
                continue

            kdata = f"{kbs}_data"
            kap0, kap1 = coll.dobj[wbs][kbs]['apex']
            ap0 = coll.ddata[kap0]['data']
            ap1 = coll.ddata[kap1]['data']

            if mtype == 'rect':
                i0 = np.argmin(np.abs(ap0 - Xm[0]))
                i1 = np.argmin(np.abs(ap1 - Xm[1]))
                ind = i0 + ap0.size * i1
            else:
                dist = (ap0 - Xm[0])**2 + (ap1 - Xm[1])**2
                ind = np.argmin(dist)

            # set nan
            if details is True:
                data = coll.ddata[kdata]['data'][:, :, ind]
            else:
                data = np.sum(coll.ddata[kdata]['data'], axis=-1)
            data[data == 0.] = np.nan

            # plot
            ax.imshow(
                data.T,
                origin='lower',
                interpolation='nearest',
                extent=extent,
                vmin=0,
                vmax=1,
            )

    return dax
