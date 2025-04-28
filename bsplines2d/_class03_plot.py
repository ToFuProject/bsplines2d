# -*- coding: utf-8 -*-


# Built-in


# Common
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datastock as ds

# specific
from . import _generic_check
from . import _class01_checks as _checks


# ################################################################
# ################################################################
#                     Main plotting
# ################################################################


def plot_bspline(
    # ressources
    coll=None,
    # inputs
    key=None,
    indbs=None,
    indt=None,
    # parameters
    knots=None,
    cents=None,
    res=None,
    plot_mesh=None,
    val_out=None,
    nan0=None,
    # plot-specific
    cmap=None,
    dax=None,
    dmargin=None,
    fs=None,
    dleg=None,
):

    # --------------
    # check input

    (
        key, keym0, keym, mtype0, mtype,
        indbs, indt,
        knots, cents, knotsi, centsi,
        plot_mesh, cmap, dleg,
    ) = _plot_bspline_check(
        coll=coll,
        key=key,
        indbs=indbs,
        indt=indt,
        knots=knots,
        cents=cents,
        plot_mesh=plot_mesh,
        cmap=cmap,
        dleg=dleg,
    )

    # --------------
    #  Prepare data

    bspline, extent, interp = _plot_bspline_prepare(
        coll=coll,
        key=key,
        keym=keym,
        mtype0=mtype0,
        mtype=mtype,
        indbs=indbs,
        indt=indt,
        knotsi=knotsi,
        centsi=centsi,
        res=res,
        val_out=val_out,
        nan0=nan0,
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

    if plot_mesh is True:
        keym = coll.dobj['bsplines'][key]['mesh']
        if mtype0 == 'polar':
            _ = coll.plot_mesh(key=keym, dleg=False)
        else:
            dax = coll.plot_mesh(key=keym, dax=dax, dleg=False)

    kax = 'cross'
    if dax.get(kax) is not None:
        ax = dax[kax]['handle']

        ax.imshow(
            bspline,
            extent=extent,
            interpolation=interp,
            origin='lower',
            aspect='equal',
            cmap=cmap,
            vmin=0.,
            vmax=1.,
        )

        if mtype0 != 'polar':
            if knots is not False:
                ax.plot(
                    knotsi[0].ravel(),
                    knotsi[1].ravel(),
                    marker='x',
                    ms=6,
                    ls='None',
                    color='k',
                )

            if cents is not False:
                ax.plot(
                    centsi[0].ravel(),
                    centsi[1].ravel(),
                    marker='o',
                    ms=6,
                    ls='None',
                    color='k',
                )

        ax.relim()
        ax.autoscale()

        # --------------
        # dleg

        if dleg is not False:
            ax.legend(**dleg)

    return dax


# ############################################################
# ############################################################
#                           plot bspline
# ############################################################


def _plot_bsplines_get_dRdZ(coll=None, km=None, meshtype=None):
    # Get minimum distances

    if meshtype == 'rect':
        kR, kZ = coll.dobj['mesh'][km]['knots']
        Rk = coll.ddata[kR]['data']
        Zk = coll.ddata[kZ]['data']
        dR = np.min(np.diff(Rk))
        dZ = np.min(np.diff(Zk))

    elif meshtype == 'tri':
        indtri = coll.ddata[coll.dobj['mesh'][km]['ind']]['data']
        kknots = coll.dobj['mesh'][km]['knots']
        Rk = coll.ddata[kknots[0]]['data']
        Zk = coll.ddata[kknots[1]]['data']
        R = Rk[indtri]
        Z = Zk[indtri]
        dist = np.mean(np.array([
            np.sqrt((R[:, 1] - R[:, 0])**2 + (Z[:, 1] - Z[:, 0])**2),
            np.sqrt((R[:, 2] - R[:, 1])**2 + (Z[:, 2] - Z[:, 1])**2),
            np.sqrt((R[:, 2] - R[:, 0])**2 + (Z[:, 2] - Z[:, 0])**2),
        ]))
        dR, dZ = dist, dist

    else:
        km2 = coll.dobj[coll._which_mesh][km]['submesh']
        meshtype = coll.dobj[coll._which_mesh][km2]['type']
        return _plot_bsplines_get_dRdZ(
            coll=coll, km=km2, meshtype=meshtype,
        )

    Rminmax = [Rk.min(), Rk.max()]
    Zminmax = [Zk.min(), Zk.max()]
    return dR, dZ, Rminmax, Zminmax


# ###############################################################
# ###############################################################
#                           checks
# ###############################################################


def _plot_bspline_check(
    coll=None,
    key=None,
    indbs=None,
    indt=None,
    knots=None,
    cents=None,
    plot_mesh=None,
    cmap=None,
    dleg=None,
):

    # key
    (
     which_mesh, which_bsplines, keym, key, cat,
     ) = _checks._get_key_mesh_vs_bplines(
        coll=coll,
        key=key,
        forcecat='bsplines',
    )

    keym0 = coll.dobj[which_bsplines][key][which_mesh]
    mtype0 = coll.dobj[which_mesh][keym0]['type']
    if mtype0 == 'polar':
        keym = coll.dobj[which_mesh][keym0]['submesh']
        mtype = coll.dobj[which_mesh][keym]['type']
    else:
        keym = keym0
        mtype = mtype0

    # knots, cents
    knots = ds._generic_check._check_var(
        knots, 'knots', default=True, types=bool,
    )
    cents = ds._generic_check._check_var(
        cents, 'cents', default=True, types=bool,
    )

    # ind_bspline
    if indbs is not None:
        indbs = coll.select_bsplines(
            key=key,
            ind=indbs,
            returnas='ind',
            return_knots=False,
            return_cents=False,
            crop=False,
        )

    _, knotsi, centsi = coll.select_bsplines(
        key=key,
        ind=indbs,
        returnas='data',
        return_knots=True,
        return_cents=True,
        crop=False,
    )

    # indt
    nt = False
    if mtype0 == 'polar':
        radius2d = coll.dobj[which_mesh][keym0]['radius2d']
        r2d_reft = coll.get_time(key=radius2d)[2]
        if r2d_reft is not None:
            nt = coll.dref[r2d_reft]['size']

    if nt is False:
        indt = None
    else:
        if indt is None:
            indt = 0
        indt = np.atleast_1d(indt).ravel()[0]

    # plot_mesh
    plot_mesh = ds._generic_check._check_var(
        plot_mesh, 'plot_mesh',
        default=True,
        types=bool,
    )

    # cmap
    if cmap is None:
        cmap = 'viridis'

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
        key, keym0, keym, mtype0, mtype,
        indbs, indt,
        knots, cents, knotsi, centsi,
        plot_mesh, cmap, dleg,
    )


# ###############################################################
# ###############################################################
#                           prepare
# ###############################################################


def _plot_bspline_prepare(
    coll=None,
    # keys
    key=None,
    keym0=None,
    keym=None,
    mtype0=None,
    mtype=None,
    # indices
    indbs=None,
    indt=None,
    # options
    res=None,
    knotsi=None,
    centsi=None,
    val_out=None,
    nan0=None,
):

    # check input
    deg = coll.dobj['bsplines'][key]['deg']

    # get dR, dZ
    dR, dZ, _, _ = _plot_bsplines_get_dRdZ(
        coll=coll, km=keym, meshtype=mtype,
    )

    # resolution of sampling
    if res is None:
        if mtype == 'rect':
            res_coef = 0.05
        else:
            res_coef = 0.25
        res = [res_coef*dR, res_coef*dZ]

    # sampling domain
    if mtype0 == 'polar':
        DR = None
        DZ = None
    else:
        knotsiR, knotsiZ = knotsi
        DR = [np.nanmin(knotsiR) + dR*1.e-10, np.nanmax(knotsiR) - dR*1.e-10]
        DZ = [np.nanmin(knotsiZ) + dZ*1.e-10, np.nanmax(knotsiZ) - dZ*1.e-10]

    # sample
    R, Z = coll.get_sample_mesh(
        key=keym,
        res=res,
        DR=DR,
        DZ=DZ,
        mode='abs', grid=True, imshow=True,
    )

    # bspline
    bspline = coll.interpolate_profile2d(
        key=key,
        R=R,
        Z=Z,
        # coefs=coefs,
        indt=indt,
        indbs=indbs,
        details=indbs is not None,
        grid=False,
        nan0=nan0,
        val_out=val_out,
        return_params=False,
    )[0]

    if indbs is None:
        if bspline.ndim == R.ndim + 1:
            assert bspline.shape[1:] == R.shape
            bspline = bspline[0, ...]
    else:
        if bspline.ndim == R.ndim + 1:
            assert bspline.shape[:-1] == R.shape
            bspline = np.nansum(bspline, axis=-1)
        elif bspline.ndim == R.ndim + 2:
            assert bspline.shape[1:-1] == R.shape
            bspline = np.nansum(bspline[0, ...], axis=-1)

    if bspline.shape != R.shape:
        import pdb; pdb.set_trace() # DB
        pass

    # extent
    if mtype0 == 'polar':
        extent = (
            R.min(), R.max(),
            Z.min(), Z.max(),
        )
    else:
        extent = (
            DR[0], DR[1],
            DZ[0], DZ[1],
        )

    # interpolation
    if deg == 0:
        interp = 'nearest'
    elif deg == 1:
        interp = 'bilinear'
    elif deg >= 2:
        interp = 'bicubic'

    return bspline, extent, interp
