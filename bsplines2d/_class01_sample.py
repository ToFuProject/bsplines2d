# -*- coding: utf-8 -*-


# Common
import numpy as np
import datastock as ds
from matplotlib import path


# tofu
from . import _generic_mesh


# ##############################################################
# ##############################################################
#                       Main
# ##############################################################


def sample_mesh(
    coll=None,
    key=None,
    res=None,
    mode=None,
    x0=None,
    x1=None,
    grid=None,
    # options
    Dx0=None,
    Dx1=None,
    submesh=None,
    # output for 2d
    imshow=None,
    return_ind=None,
    in_mesh=None,
    # store
    store=None,
    kx0=None,
    kx1=None,
):

    # -------------
    # check inputs

    (
        key, nd, mtype,
        mode, res, submesh,
        grid, imshow, store,
    ) = _check(
        coll=coll,
        key=key,
        mode=mode,
        res=res,
        submesh=submesh,
        grid=grid,
        imshow=imshow,
        store=store,
        x0=x0,
        x1=x1,
    )

    # ------------
    # sample

    if nd == '1d':

        # check
        knots, Dx = _check_1d(
            coll=coll,
            key=key,
            Dx=Dx0,
            mode=mode,
        )

        # sample
        x0 = _sample_1d(
            res=res,
            mode=mode,
            knots=knots,
            Dx=Dx,
        )

        # prepare storing
        dref, ddata = _store_1d(
            coll=coll,
            key=key,
            x0=x0,
            # store
            store=store,
            kx0=kx0,
        )

    else:

        # check
        knots0, knots1, x0, x1, Dx0, Dx1 = _check_2d(
            coll=coll,
            key=key,
            x0=x0,
            x1=x1,
            Dx0=Dx0,
            Dx1=Dx1,
            mode=mode,
            mtype=mtype,
        )

        # sample
        x0, x1, ind = _sample_2d(
            coll=coll,
            key=key,
            res=res,
            mode=mode,
            knots0=knots0,
            knots1=knots1,
            x0=x0,
            x1=x1,
            Dx0=Dx0,
            Dx1=Dx1,
            # options
            grid=grid,
            imshow=imshow,
            in_mesh=in_mesh,
        )

        # prepare storing
        dref, ddata = _store_2d(
            coll=coll,
            key=key,
            x0=x0,
            x1=x1,
            ind=ind,
            # store
            store=store,
            kx0=kx0,
            kx1=kx1,
        )

    # -------------
    # format output

    # -------------
    # store

    if store is True:

        _store(coll=coll, dref=dref, ddata=ddata)

    # --------
    # return

    return ddata


# ##################################################################
# ##################################################################
#                       checks
# ##################################################################


def _check(
    coll=None,
    key=None,
    mode=None,
    res=None,
    submesh=None,
    grid=None,
    imshow=None,
    store=None,
    x0=None,
    x1=None,
):

    # ---------
    # mesh

    # submesh
    submesh = ds._generic_check._check_var(
        submesh, 'submesh',
        types=bool,
        default=False,
    )

    # key
    wm = coll._which_mesh
    key, _, cat = _generic_mesh._get_key_mesh_vs_bplines(
        coll=coll,
        key=key,
        which=wm,
    )

    nd = coll.dobj[cat][key]['nd']
    mtype = coll.dobj[cat][key]['type']

    if nd == '1d' and submesh is True:
        if coll.dobj[wm][key]['submesh'] is not None:
            key = coll.dobj[wm][key]['submesh']
            nd = coll.dobj[cat][key]['nd']
            mtype = coll.dobj[cat][key]['type']

    # ----------
    # resolution

    # mode
    mode = ds._generic_check._check_var(
        mode, 'mode',
        types=str,
        default='abs',
    )

    # res
    res = _get_res(
        coll=coll,
        key=key,
        res=res,
        nd=nd,
        mtype=mtype,
    )

    # ---------
    # parameters

    if nd == '1d':
        lok = None

    else:

        # grid
        grid = ds._generic_check._check_var(
            grid, 'grid',
            types=bool,
            default=False,
        )

        # imshow
        imshow = ds._generic_check._check_var(
            imshow, 'imshow',
            types=bool,
            default=False,
        )

        lok = [False]
        if grid is False and imshow is False:
            if x0 is None and x1 is None:
                lok.append(True)

    # --------
    # store

    store = ds._generic_check._check_var(
        store, 'store',
        types=bool,
        default=False,
        allowed=lok,
    )

    return (
        key, nd, mtype,
        mode, res, submesh,
        grid, imshow, store,
    )


def _check_1d(
    coll=None,
    key=None,
    Dx=None,
    mode=None,
):

    wm = coll._which_mesh

    # -------------
    # knots

    kknots = coll.dobj[wm][key]['knots'][0]
    knots = coll.ddata[kknots]['data']

    # custom DR or DZ for mode='abs' only
    Dx = _Dx(Dx)

    if Dx is not None and mode != 'abs':
        msg = "Custom Dx can only be provided with mode = 'abs'!"
        raise Exception(msg)

    return knots, Dx


def _check_2d(
    coll=None,
    key=None,
    x0=None,
    x1=None,
    Dx0=None,
    Dx1=None,
    mode=None,
    mtype=None,
):

    wm = coll._which_mesh

    # -----------
    # Parameters

    # R, Z
    if x0 is None and x1 is None:
        pass
    elif x0 is None and np.isscalar(x1):
        pass
    elif x1 is None and np.isscalar(x0):
        pass
    else:
        msg = (
            "For mesh discretisation, (x0, x1) can be either:\n"
            "\t- (None, None): will be created\n"
            "\t- (scalar, None): A vertical line will be created\n"
            "\t- (None, scalar): A horizontal line will be created\n"
        )
        raise Exception(msg)

    # -------------
    # x0, x1

    k0, k1 = coll.dobj[wm][key]['knots']
    if mtype == 'rect':
        knots0 = coll.ddata[k0]['data']
        knots1 = coll.ddata[k1]['data']

        # custom R xor Z for vertical / horizontal lines only
        if x0 is None and x1 is not None:
            x0 = knots0
        if x1 is None and x0 is not None:
            x1 = knots1

    else:
        knots0 = coll.ddata[k0]['data']
        knots1 = coll.ddata[k1]['data']

    # custom DR or DZ for mode='abs' only
    Dx0 = _DRZ(Dx=Dx0, size=None, Dx_name='Dx0')
    Dx1 = _DRZ(
        Dx=Dx1,
        size=len(Dx0) if Dx0 is not None else None,
        Dx_name='Dx1',
    )

    if (Dx0 is not None or Dx1 is not None) and mode != 'abs':
        msg = (
            "Custom Dx0 or Dx1 can only be provided with mode = 'abs'!\n"
            f"\t- mode: {mode}\n"
            f"\t- Dx0: {Dx0}\n"
            f"\t- Dx1: {Dx1}\n"
        )
        raise Exception(msg)

    return knots0, knots1, x0, x1, Dx0, Dx1


# #############
# check Dx
# #############


def _Dx(Dx=None):

    if Dx is None:
        return None

    c0 = (
        hasattr(Dx, '__iter__')
        and len(Dx) == 2
        and all([
            rr is None or (np.isscalar(rr) and not np.isnan(rr))
            for rr in Dx
        ])
    )
    if not c0:
        msg = 'Arg Dx must be an iterable of 2 scalars!'
        raise Exception(msg)

    for ii in range(2):
        if Dx[ii] is None:
            Dx[ii] = -np.inf if ii == 00 else np.inf

    return Dx


def _DRZ(Dx=None, size=None, Dx_name=None):

    if Dx is None:
        return None

    if None not in Dx and np.all(np.isfinite(Dx)) and np.size(Dx) > 2:
        Dx = np.atleast_1d(Dx).ravel()
        if size is not None and Dx.size != size:
            msg = (
                f"Arg {Dx_name} has wrong size!\n"
                f"\t- expected: {size}\n"
                f"\t- Provided: {Dx.size}"
            )
            raise Exception(msg)

    else:
        Dx = _Dx(Dx)

    return Dx


# ################################################################
# ################################################################
#                       utility
# ################################################################


def _get_res(
    coll=None,
    key=None,
    res=None,
    nd=None,
    mtype=None,
):

    wm = coll._which_mesh
    if nd == '1d':

        if res is None:
            kknots = coll.dobj[wm][key]['knots'][0]
            res = np.min(np.diff(coll.ddata[kknots]['data']))

        if not (np.isscalar(res) and res > 0.):
            msg = (
                "Arg res must be a positive float!\n"
                f"Provided: {res}"
            )
            raise Exception(msg)

    else:

        if res is None:
            if mtype == 'rect':
                kR, kZ = coll.dobj[wm][key]['knots']
                res = min(
                    np.min(np.diff(coll.ddata[kR]['data'])),
                    np.min(np.diff(coll.ddata[kZ]['data'])),
                )

            elif mtype == 'tri':
                res = 0.02

        # check len() = 2
        if np.isscalar(res):
            res = [res, res]

        c0 = (
            isinstance(res, list)
            and len(res) == 2
            and all([np.isscalar(rr) and rr > 0 for rr in res])
        )
        if not c0:
            msg = (
                f"Arg res must be a list of 2 positive floats!\n"
                "Provided: {res}"
            )
            raise Exception(msg)

    return res


# ##################################################################
# ##################################################################
#                       sample 1d
# ##################################################################


def _sample_1d(
    res=None,
    mode=None,
    knots=None,
    Dx=None,
):

    kmin, kmax = knots.min(), knots.max()

    if mode == 'abs':
        nx = int(np.ceil((kmax - kmin) / res))
        xx = np.linspace(kmin, kmax, nx)

    else:
        nx = int(np.ceil(1./res))
        kx = np.linspace(0, 1, nx, endpoint=False)[None, :]
        xx = np.concatenate((
            (knots[:-1, None] + kx*np.diff(knots)[:, None]).ravel(),
            knots[-1:],
        ))

    # Dx
    if Dx is not None:
        ind = (xx >= Dx[0]) & (xx < Dx[1])
        xx = xx[ind]

    return xx


# ##################################################################
# ##################################################################
#                       sample 2d
# ##################################################################


def _sample_2d(
    coll=None,
    key=None,
    res=None,
    mode=None,
    knots0=None,
    knots1=None,
    x0=None,
    x1=None,
    Dx0=None,
    Dx1=None,
    # options
    grid=None,
    imshow=None,
    in_mesh=None,
):

    # --------
    # compute

    min0, max0 = knots0.min(), knots0.max()
    min1, max1 = knots1.min(), knots1.max()

    if mode == 'abs':
        if x0 is None:
            n0 = int(np.ceil((max0 - min0) / res[0]))
            x0 = np.linspace(min0, max0, n0)
        if x1 is None:
            n1 = int(np.ceil((max1 - min1) / res[1]))
            x1 = np.linspace(min1, max1, n1)
    else:
        if x0 is None:
            n0 = int(np.ceil(1./res[0]))
            kx0 = np.linspace(0, 1, n0, endpoint=False)[None, :]
            x0 = np.concatenate((
                (knots0[:-1, None] + kx0*np.diff(knots0)[:, None]).ravel(),
                knots0[-1:],
            ))
        if x1 is None:
            n1 = int(np.ceil(1./res[1]))
            kx1 = np.linspace(0, 1, n1, endpoint=False)[None, :]
            x1 = np.concatenate((
                (knots1[:-1, None] + kx1*np.diff(knots1)[:, None]).ravel(),
                knots1[-1:],
            ))

    if np.isscalar(x0):
        x0 = np.full(x1.shape, x0)
    if np.isscalar(x1):
        x1 = np.full(x0.shape, x1)

    # -----------
    # prepare ind

    # x0, x1
    if grid or in_mesh or Dx0 is not None:
        x02 = np.repeat(x0[:, None], x1.size, axis=1)
        x12 = np.repeat(x1[None, :], x0.size, axis=0)

        if in_mesh or Dx0 is not None:
            x0f = x02.ravel()
            x1f = x12.ravel()
            sh = x02.shape

    # ind
    if Dx0 is not None or in_mesh is True:
        ind = np.ones((x0.size, x1.size), dtype=bool)
    else:
        ind = None

    # ---------
    # Dx0, Dx1

    if Dx0 is not None and len(Dx0) == 2:
        ind = ind & (
            (x02 >= Dx0[0]) & (x02 < Dx0[1])
            & (x12 >= Dx1[0]) & (x12 < Dx1[1])
        )

    # -----------
    # indices

    # mesh outline
    if in_mesh is True:
        dout = coll.get_mesh_outline(key=key)
        pa = path.Path(np.array([dout['x0']['data'], dout['x1']['data']]).T)
        ind = ind & pa.contains_points(np.array([x0f, x1f]).T).reshape(sh)

    # Dx0, Dx1
    if Dx0 is not None and len(Dx0) > 2:
        pa = path.Path(np.array([Dx0, Dx1]).T)
        ind = ind & pa.contains_points(np.array([x0f, x1f]).T).reshape(sh)

    # ------------
    # grid

    if grid is True:
        x0, x1 = x02, x12

    # ----------
    # imshow

    if imshow is True:
        if ind is not None:
            ind = ind.T
        if grid is True:
            x0 = x0.T
            x1 = x1.T

    return x0, x1, ind


# ##################################################################
# ##################################################################
#                       store
# ##################################################################


def _store_1d(
    coll=None,
    key=None,
    x0=None,
    # store
    store=None,
    kx0=None,
):

    # -----------
    # check key

    if store is True:
        lout = list(coll.ddata.keys())
    else:
        lout = []

    kx0 = ds._generic_check._check_var(
        kx0, 'kx0',
        types=str,
        default=f'{key}_x0',
        excluded=lout,
    )

    # -----------
    # format dict

    wm = coll._which_mesh
    kknots = coll.dobj[wm][key]['knots'][0]

    # dref
    kr = f'{kx0}_n'
    dref = {kr: {'size': x0.size}}

    # ddata
    ddata = {
        'x0': {
            'key': kx0,
            'data': x0,
            'dim': coll.ddata[kknots]['dim'],
            'quant': coll.ddata[kknots]['quant'],
            'name': coll.ddata[kknots]['name'],
            'units': coll.ddata[kknots]['units'],
        },
    }

    if store:
        ddata['x0']['ref'] = kr

    return dref, ddata


def _store_2d(
    coll=None,
    key=None,
    x0=None,
    x1=None,
    ind=None,
    # store
    store=None,
    kx0=None,
    kx1=None,
):


    # -----------
    # check key

    if store is True:
        lout = list(coll.ddata.keys())
    else:
        lout = []

    kx0 = ds._generic_check._check_var(
        kx0, 'kx0',
        types=str,
        default=f'{key}_x0_temp',
        excluded=lout,
    )

    kx1 = ds._generic_check._check_var(
        kx1, 'kx1',
        types=str,
        default=f'{key}_x1_temp',
        excluded=lout,
    )

    # -----------
    # format dict

    wm = coll._which_mesh
    kk0, kk1 = coll.dobj[wm][key]['knots']

    # dref
    k0r, k1r = f'{kx0}_n', f'{kx1}_n'
    dref = {
        k0r: {'size': x0.size},
        k1r: {'size': x1.size},
    }

    # ddata
    ddata = {
        'x0': {
            'key': kx0,
            'data': x0,
            'dim': coll.ddata[kk0]['dim'],
            'quant': coll.ddata[kk0]['quant'],
            'name': coll.ddata[kk0]['name'],
            'units': coll.ddata[kk0]['units'],
        },
        'x1': {
            'key': kx1,
            'data': x1,
            'dim': coll.ddata[kk1]['dim'],
            'quant': coll.ddata[kk1]['quant'],
            'name': coll.ddata[kk1]['name'],
            'units': coll.ddata[kk1]['units'],
        },
    }

    if ind is not None:
        ddata['ind'] = {
            'key': f'{key}_ind_temp',
            'data': ind,
        }

    if store is True:
        ddata['x0']['ref'] = k0r
        ddata['x1']['ref'] = k1r
        if ind is not None:
            ddata['ind']['ref'] = (k0r, k1r)

    return dref, ddata


def _store(coll=None, dref=None, ddata=None):

    # dref
    for k0, v0 in dref.items():
        coll.add_ref(key=k0, **v0)

    # ddata
    for k0, v0 in ddata.items():
        coll.add_data(**v0)