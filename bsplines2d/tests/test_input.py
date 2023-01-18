"""
This module contains tests for tofu.geom in its structured version
"""

# Built-in
import os
import itertools as itt


# Standard
import numpy as np
import matplotlib.pyplot as plt


# tofu-specific



_HERE = os.path.abspath(os.path.dirname(__file__))
_PATH_DATA = os.path.join(_HERE, 'test_data')
_DFNAME = {
    'WEST_tri': 'WEST_trimesh.npz',
    'quad_tri': 'trimesh_quad.npz',
    'WEST_poly': 'WEST_Poly.npz',
}

_DDATA = {}
for k0, fname in _DFNAME.items():
    _PFE = os.path.join(_PATH_DATA, fname)
    dd = dict(np.load(_PFE, allow_pickle=True))
    if k0 == 'WEST':
        _DDATA[k0] = dd
    else:
        _DDATA[k0] = {}
        for k1, v1 in dd.items():
            _DDATA[k0][k1] = v1.tolist()


#######################################################
#
#     Basic instanciation
#
#######################################################


def _add_1d_knots_uniform(bsplines, key=None):
    bsplines.add_mesh_1d(key=key, knots=np.linspace(0, 10, 11))

    
def _add_1d_knots_variable(bsplines, key=None):
    bsplines.add_mesh_1d(key=key, knots=np.r_[0, 1, 4, 7, 10], uniform=False)
    

def _add_rect_uniform(bsplines, key=None):
    # add uniform rect mesh
    bsplines.add_mesh_2d_rect(key=key, domain=[[2, 3], [-1, 1]], res=0.1)


def _add_rect_variable(bsplines, key=None):
    # add variable rect mesh
    bsplines.add_mesh_2d_rect(
        key=key,
        domain=[[2, 2.3, 2.6, 3], [-1, 0., 1]],
        res=[[0.2, 0.1, 0.1, 0.2], [0.2, 0.1, 0.2]],
    )


def _add_rect_variable_crop(bsplines, key=None):
    # add variable rect mesh
    bsplines.add_mesh_2d_rect(
        key=key,
        domain=[[2, 2.3, 2.6, 3], [-1, 0., 1]],
        res=[[0.2, 0.1, 0.1, 0.2], [0.2, 0.1, 0.2]],
        crop_poly=_DDATA['WEST_poly']['Poly'],
    )


def _add_rect_crop_from_knots(bsplines, key=None):
    # add variable rect mesh
    bsplines.add_mesh_2d_rect(
        key=key,
        knots0=np.linspace(2, 3, 11),
        knots1=np.linspace(-1, 1, 11),
        crop_poly=_DDATA['WEST_poly']['Poly'],
    )

    
def _add_rect_variable_crop_from_knots(bsplines, key=None):
    # add variable rect mesh
    bsplines.add_mesh_2d_rect(
        key=key,
        knots0=np.r_[np.linspace(2, 2.4, 5), np.r_[2.5, 2.7, 2.8, 3.]],
        knots1=np.linspace(-1, 1, 11),
        crop_poly=_DDATA['WEST_poly']['Poly'],
    )


def _add_tri_ntri1(bsplines, key=None):
    
    knots = np.array([
        _DDATA['WEST_tri']['pts_x0'],
        _DDATA['WEST_tri']['pts_x1'],
    ]).T
    
    bsplines.add_mesh_2d_tri(
        key=key,
        knots=knots,
        indices=_DDATA['WEST_tri']['indices'],
    )


def _add_tri_ntri2(bsplines, key=None):
    
    knots = np.array([
        _DDATA['quad_tri']['pts_x0'],
        _DDATA['quad_tri']['pts_x1'],
    ]).T
    
    bsplines.add_mesh_2d_tri(
        key=key,
        knots=knots,
        indices=_DDATA['quad_tri']['indices'],
    )


def _add_tri_delaunay(bsplines, key=None):
    
    bsplines.add_mesh_2d_tri(
        key=key,
        pts_x0=_DDATA['quad_tri']['pts_x0'],
        pts_x1=_DDATA['quad_tri']['pts_x1'],
    )


def _add_mesh_1d_subkey_fixed(bsplines, key=None, keybs=None):
    
    ka = bsplines.dobj['bsplines'][keybs]['apex']
    ap = bsplines.ddata[ka]['data']
    rho = 1 - np.exp(-ap**2)
    
    bsplines.add_data(
        key='rho1d',
        data=rho,
        ref=keybs,
    )
    
    bsplines.add_mesh_1d(
        key=key,
        knots=np.linspace(0, 20, 5),
        subkey='rho1d',
    )


def _add_polar1(bsplines, key=None):
    """ Time-independent """

    kR, kZ = bsplines.dobj['bsplines']['m2-bs1']['apex']
    R = bsplines.ddata[kR]['data']
    Z = bsplines.ddata[kZ]['data']
    RR = np.repeat(R[:, None], Z.size, axis=1)
    ZZ = np.repeat(Z[None, :], R.size, axis=0)
    rho = (RR - 2.5)**2/0.08 + (ZZ - 0)**2/0.35

    bsplines.add_data(
        key='rho1',
        data=rho,
        ref='m2-bs1',
        unit='',
        dim='',
        quant='rho',
        name='rho',
    )

    bsplines.add_mesh_2d_polar(
        key=key,
        radius=np.linspace(0, 1.2, 7),
        angle=None,
        radius2d='rho1',
    )


def _add_polar2(bsplines, key=None):
    """ Time-dependent """

    kR, kZ = bsplines.dobj['bsplines']['m2-bs1']['apex']
    R = bsplines.ddata[kR]['data']
    Z = bsplines.ddata[kZ]['data']
    RR = np.repeat(R[:, None], Z.size, axis=1)
    ZZ = np.repeat(Z[None, :], R.size, axis=0)

    rho = (RR - 2.5)**2/0.08 + (ZZ - 0)**2/0.35
    angle = np.arctan2(ZZ/2., (RR - 2.5))

    nt = 11
    t = np.linspace(30, 40, nt)
    rho = rho[None, ...] + 0.1*np.cos(t)[:, None, None]**2
    angle = angle[None, ...] + 0.01*np.sin(t)[:, None, None]**2


    if 'nt' not in bsplines.dref.keys():
        bsplines.add_ref(
            key='nt',
            size=nt,
        )

    if 't' not in bsplines.ddata.keys():
        bsplines.add_data(
            key='t',
            data=t,
            ref=('nt',),
            dim='time',
        )

    if 'rho2' not in bsplines.ddata.keys():
        bsplines.add_data(
            key='rho2',
            data=rho,
            ref=('nt', 'm2-bs1'),
            unit='',
            dim='',
            quant='rho',
            name='rho',
        )

    if 'angle2' not in bsplines.ddata.keys():
        bsplines.add_data(
            key='angle2',
            data=angle,
            ref=('nt', 'm2-bs1'),
            unit='rad',
            dim='',
            quant='angle',
            name='theta',
        )

    # ang
    if key == 'm6':
        ang = np.pi*np.r_[-3./4., -1/4, 0, 1/4, 3/4]
    else:
        ang = None

    # mesh
    bsplines.add_mesh_2d_polar(
        key=key,
        radius=np.linspace(0, 1.2, 7),
        angle=ang,
        radius2d='rho2',
        angle2d='angle2',
    )


def _add_bsplines(bs, key=None, nd=None, kind=None, deg=None, angle=None):

    lkm = _get_mesh(bs, nd=nd, kind=kind)
    if isinstance(deg, int):
        deg = [deg]
    elif deg is None:
        deg = [0, 1, 2, 3]

    ddeg = {
        None: [ii for ii in [0, 1, 2, 3] if ii in deg],
        'rect': [ii for ii in [0, 1, 2, 3] if ii in deg],
        'tri': [ii for ii in [0, 1] if ii in deg],
        'polar': [ii for ii in [0, 1, 2, 3] if ii in deg],
    }

    for km in lkm:
        mtype = bs.dobj[bs._which_mesh][km]['type']
        for dd in ddeg[mtype]:
            bs.add_bsplines(key=k0, deg=dd)

def _get_mesh(bsplines, nd=None, kind=None):
    return [
        k0 for k0, v0 in bsplines.dobj[bsplines._which_mesh].items()
        if nd in [None, v0['nd']]
        and kind in [None, v0['type']]
    ]


def _select_mesh_elements(bsplines, nd=None, kind=None):
    lkm = _get_mesh(bsplines, nd=nd, kind=kind)
    
    if nd == '1d':
        lcrop = [None]
        lind = [None, 0, [0, 3]]
        lneigh = [False]
        
    elif kind == 'rect':
        lcrop = [False, True]
        lind = [None, 0, [0, 3]]
        lneigh = [False, True]
        
    else:
        lcrop = [None]
        lind = [None, 0, [0, 3]]
        lneigh = [False, True]
    
    lel = ['knots', 'cents']
    lret = ['ind', 'data']
    
    for km in lkm:
        
        for comb in itt.product(lcrop, lel, lind, lret, lneigh):
            out = bsplines.select_mesh_elements(
                key=km,
                ind=comb[2],
                elements=comb[1],
                returnas=comb[3],
                return_neighbours=comb[4],
                crop=comb[0],
            )

def _sample_mesh(bsplines, nd=None, kind=None):
    lkm = _get_mesh(bsplines, nd=nd, kind=kind)
    
    lres = [0.1, 0.3]
    lmode = ['abs', 'rel']
    limshow = [False, True]
    
    for km in lkm:
        
        for comb in itt.product(lres, lmode, limshow):
            out = bsplines.get_sample_mesh(
                key=km,
                res=comb[0],
                mode=comb[1],
                grid=None,
                imshow=comb[2]
            )


def _plot_mesh(bsplines, nd=None, kind=None):
    lkm = _get_mesh(bsplines, nd=nd, kind=kind)
    
    lk = [None, 2, [2, 3]]
    lc = [None, 2, [2, 3]]
    if kind == 'rect':
        lcrop = [False, True]
    else:
        lcrop = [False]
    
    for km in lkm:
        
        for comb in itt.product(lk, lc, lcrop):
            dax = bsplines.plot_mesh(
                key=km,
                ind_knot=comb[0],
                ind_cent=comb[1],
                crop=comb[2],
            )
            plt.close('all')


def _add_data_fix(bsplines, key):

    kdata = f'{key}-data-fix'
    shape = bsplines.dobj['bsplines'][key]['shape']
    data = np.random.random(shape)

    bsplines.add_data(
        key=kdata,
        data=data,
        ref=key,
    )
    return kdata


def _add_data_var(bsplines, key):

    kdata = f'{key}-data-var'
    shape = bsplines.dobj['bsplines'][key]['shape']
    t = bsplines.ddata['t']['data']
    tsh = tuple([t.size] + [1 for ii in shape])
    data = np.cos(t.reshape(tsh)) * np.random.random(shape)[None, ...]

    bsplines.add_data(
        key=kdata,
        data=data,
        ref=('nt', key),
    )
    return kdata