# -*- coding: utf-8 -*-


# Built-in
import warnings

# Common
import numpy as np
from scipy.spatial import ConvexHull
from matplotlib.path import Path
from contourpy import contour_generator
import datastock as ds


# tofu
from . import _generic_mesh
from . import _utils_bsplines
# from . import _class02_checks as _checks
from . import _class02_bsplines_rect
from . import _class02_bsplines_tri
from . import _class02_bsplines_polar
from . import _class02_bsplines_1d


# #############################################################################
# #############################################################################
#                           main
# #############################################################################


def sample_mesh(
    coll=None,
    key=None,
    res=None,
    mode=None,
    x0=None,
    x1=None,
    Dx0=None,
    Dx1=None,
    grid=None,
    imshow=None,
):

    # -------------
    # check inputs
    
    # key
    key, _, cat = _generic_mesh._get_key_mesh_vs_bplines(
        coll=coll,
        key=key,
        which=coll._which_mesh,
    )
     
    nd = coll.dobj[cat][key]['nd']
    mtype = coll.dobj[cat][key]['type']
    
    # mode
    mode = ds._generic_check._check_var(
        mode, 'mode',
        types=str,
        default='abs',
    )
    
    # res
    if res is None:
        res = _get_sample_mesh_res(coll=coll, keym=key, mtype=mtype)
    
    # ------------
    # sample
    
    if nd == '1d':
        
        # check
        key, mode, Dx, knots = _sample_mesh_check_1d(
            coll=coll,
            key=key,
            res=res,
            mode=mode,
            Dx=Dx0,
        )
        
        # compute
        return sample_mesh_1d(
            res=res,
            mode=mode,
            Dx=Dx,
            knots=knots,
        )
        
    else:

        # check
        (
            key, res, mode, grid, imshow, x0, x1, Dx0, Dx1, x0k, x1k,
        ) = _sample_mesh_check_2d(
            coll=coll,
            key=key,
            res=res,
            mode=mode,
            grid=grid,
            imshow=imshow,
            R=x0,
            Z=x1,
            DR=Dx0,
            DZ=Dx1,
        )
            
        return sample_mesh_2d(
            res=res,
            mode=mode,
            x0=x0,
            x1=x1,
            Dx0=Dx0,
            Dx1=Dx1,
            x0k=x0k,
            x1k=x1k,
            grid=grid,
            imshow=imshow,
        )


# #############################################################################
# #############################################################################
#                           sub-routines
# #############################################################################


def _get_sample_mesh_res(
    coll=None,
    keym=None,
    nd=None,
    mtype=None,
):

    if nd == '1d':
        kknots = coll.dobj[coll._which_msp][keym]['knots'][0]
        res = np.min(np.diff(coll.ddata[kknots]['data']))
        
    elif mtype == 'rect':
        kR, kZ = coll.dobj[coll._which_mesh][keym]['knots']
        res = min(
            np.min(np.diff(coll.ddata[kR]['data'])),
            np.min(np.diff(coll.ddata[kZ]['data'])),
        )
    elif mtype == 'tri':
        res = 0.02

    elif mtype == 'polar':
        keyr2d = coll.dobj[coll._which_mesh][keym]['radius2d']
        keybs0 = coll.ddata[keyr2d]['bsplines']
        keym0 = coll.dobj['bsplines'][keybs0]['mesh']
        mtype0 = coll.dobj[coll._which_mesh][keym0]['type']
        res = _get_sample_mesh_res(coll=coll, keym=keym0, mtype=mtype0)

    else:
        raise Exception("Wrong mtype: {mtype}")
        
    return res


def _sample_mesh_check_1d(
    coll=None,
    key=None,
    res=None,
    mode=None,
    Dx=None,
    mtype=None,
):

    # -----------
    # Parameters

    # for polar mesh => sample underlying mesh
    if len(coll.dobj[coll._which_mesh][key]['shape-c']) > 1:
        msg = "Wrong mesh dimension!"
        raise Exception(msg)

    if not (np.isscalar(res) and res > 0.):
        msg = f"Arg res must be a positive float!\nProvided: {res}"
        raise Exception(msg)

    # -------------
    # knots

    kknots = coll.dobj[coll._which_mesh][key]['knots'][0]
    knots = coll.ddata[kknots]['data']

    # custom DR or DZ for mode='abs' only
    if Dx is not None:
        if mode != 'abs':
            msg = "Custom Dx can only be provided with mode = 'abs'!"
            raise Exception(msg)

            c0 = (
                hasattr(Dx, '__iter__')
                and len(Dx) == 2
                and all([
                    rr is None or (np.isscalar(rr) and np.isfinite(rr))
                    for rr in Dx
                ])
            )
            if not c0:
                msg = 'Arg Dx must be an iterable of 2 scalars!'
                raise Exception(msg)

    if Dx is None:
        Dx = [knots.min(), knots.max()]

    return key, mode, Dx, knots


def _sample_mesh_check_2d(
    coll=None,
    key=None,
    mtype=None,
    res=None,
    mode=None,
    grid=None,
    imshow=None,
    R=None,
    Z=None,
    DR=None,
    DZ=None,
):

    # -----------
    # Parameters

    # for polar mesh => sample underlying mesh
    if mtype == 'polar':
        key = coll.dobj[coll._which_mesh][key]['submesh']
        mtype = coll.dobj[coll._which_mesh][key]['type']

    if np.isscalar(res):
        res = [res, res]
    c0 = (
        isinstance(res, list)
        and len(res) == 2
        and all([np.isscalar(rr) and rr > 0 for rr in res])
    )
    if not c0:
        msg = f"Arg res must be a list of 2 positive floats!\nProvided: {res}"
        raise Exception(msg)

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

    # R, Z
    if R is None and Z is None:
        pass
    elif R is None and np.isscalar(Z):
        pass
    elif Z is None and np.isscalar(R):
        pass
    else:
        msg = (
            "For mesh discretisation, (R, Z) can be either:\n"
            "\t- (None, None): will be created\n"
            "\t- (scalar, None): A vertical line will be created\n"
            "\t- (None, scalar): A horizontal line will be created\n"
        )
        raise Exception(msg)

    # -------------
    # R, Z

    if mtype == 'rect':
        kR, kZ = coll.dobj[coll._which_mesh][key]['knots']
        Rk = coll.ddata[kR]['data']
        Zk = coll.ddata[kZ]['data']

        # custom R xor Z for vertical / horizontal lines only
        if R is None and Z is not None:
            R = Rk
        if Z is None and R is not None:
            Z = Zk
    else:
        kknots = coll.dobj[coll._which_mesh][key]['knots']
        Rk = coll.ddata[kknots[0]]['data']
        Zk = coll.ddata[kknots[1]]['data']

    # custom DR or DZ for mode='abs' only
    if DR is not None or DZ is not None:
        if mode != 'abs':
            msg = "Custom DR or DZ can only be provided with mode = 'abs'!"
            raise Exception(msg)

        for DD, DN in [(DR, 'DR'), (DZ, 'DZ')]:
            if DD is not None:
                c0 = (
                    hasattr(DD, '__iter__')
                    and len(DD) == 2
                    and all([
                        rr is None or (np.isscalar(rr) and np.isfinite(rr))
                        for rr in DD
                    ])
                )
                if not c0:
                    msg = f'Arg {DN} must be an iterable of 2 scalars!'
                    raise Exception(msg)

    if DR is None:
        DR = [Rk.min(), Rk.max()]
    if DZ is None:
        DZ = [Zk.min(), Zk.max()]

    return key, res, mode, grid, imshow, R, Z, DR, DZ, Rk, Zk


def sample_mesh_1d(
    res=None,
    mode=None,
    Dx=None,
    knots=None,
):

    if mode == 'abs':
        nx = int(np.ceil((Dx[1] - Dx[0]) / res))
        xx = np.linspace(Dx[0], Dx[1], nx)

    else:
        nx = int(np.ceil(1./res))
        kx = np.linspace(0, 1, nx, endpoint=False)[None, :]
        xx = np.concatenate((
            (knots[:-1, None] + kx*np.diff(knots)[:, None]).ravel(),
            knots[-1:],
        ))

    return xx
    
    
def sample_mesh_2d(
    res=None,
    mode=None,
    x0=None,
    x1=None,
    Dx0=None,
    Dx1=None,
    x0k=None,
    x1k=None,
    grid=None,
    imshow=None,
):

    # compute
    if mode == 'abs':
        if x0 is None:
            n0 = int(np.ceil((Dx0[1] - Dx0[0]) / res[0]))
            x0 = np.linspace(Dx0[0], Dx0[1], n0)
        if x1 is None:
            n1 = int(np.ceil((Dx1[1] - Dx1[0]) / res[1]))
            x1 = np.linspace(Dx1[0], Dx1[1], n1)
    else:
        if x0 is None:
            n0 = int(np.ceil(1./res[0]))
            kx0 = np.linspace(0, 1, n0, endpoint=False)[None, :]
            x0 = np.concatenate((
                (x0k[:-1, None] + kx0*np.diff(x0k)[:, None]).ravel(),
                x0k[-1:],
            ))
        if x1 is None:
            n1 = int(np.ceil(1./res[1]))
            kx1 = np.linspace(0, 1, n1, endpoint=False)[None, :]
            x1 = np.concatenate((
                (x1k[:-1, None] + kx1*np.diff(x1k)[:, None]).ravel(),
                x1k[-1:],
            ))

    if np.isscalar(x0):
        x0 = np.full(x1.shape, x0)
    if np.isscalar(x1):
        x1 = np.full(x0.shape, x1)

    # ------------
    # grid

    if grid is True:
        n1 = x1.size
        n0 = x0.size
        if imshow is True:
            x0 = np.tile(x0, (n1, 1))
            x1 = np.repeat(x1[:, None], n0, axis=1)
        else:
            x0 = np.repeat(x0[:, None], n1, axis=1)
            x1 = np.tile(x1, (n0, 1))

    return x0, x1