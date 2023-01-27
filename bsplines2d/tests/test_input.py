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
    bsplines.add_mesh_1d(key=key, knots=np.linspace(0, 10, 11), units='eV')


def _add_1d_knots_variable(bsplines, key=None):
    bsplines.add_mesh_1d(
        key=key, knots=np.r_[0, 1, 4, 7, 10], uniform=False, units='A',
    )


def _add_rect_uniform(bsplines, key=None):
    # add uniform rect mesh
    bsplines.add_mesh_2d_rect(
        key=key, domain=[[2, 3], [-1, 1]], res=0.1, units='GHz',
    )


def _add_rect_variable(bsplines, key=None):
    # add variable rect mesh
    bsplines.add_mesh_2d_rect(
        key=key,
        domain=[[2, 2.3, 2.6, 3], [-1, 0., 1]],
        res=[[0.2, 0.1, 0.1, 0.2], [0.2, 0.1, 0.2]],
        units='m',
    )


def _add_rect_variable_crop(bsplines, key=None):
    # add variable rect mesh
    bsplines.add_mesh_2d_rect(
        key=key,
        domain=[[2, 2.3, 2.6, 3], [-1, 0., 1]],
        res=[[0.2, 0.1, 0.1, 0.2], [0.2, 0.1, 0.2]],
        crop_poly=_DDATA['WEST_poly']['Poly'],
        units='cm',
    )


def _add_rect_crop_from_knots(bsplines, key=None):
    # add variable rect mesh
    bsplines.add_mesh_2d_rect(
        key=key,
        knots0=np.linspace(2, 3, 11),
        knots1=np.linspace(-1, 1, 11),
        crop_poly=_DDATA['WEST_poly']['Poly'],
        units='mm',
    )


def _add_rect_variable_crop_from_knots(bsplines, key=None):
    # add variable rect mesh
    bsplines.add_mesh_2d_rect(
        key=key,
        knots0=np.r_[np.linspace(2, 2.4, 5), np.r_[2.5, 2.7, 2.8, 3.]],
        knots1=np.linspace(-1, 1, 11),
        crop_poly=_DDATA['WEST_poly']['Poly'],
        units='km',
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
        units='dm',
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
        units='um',
    )


def _add_tri_delaunay(bsplines, key=None):

    bsplines.add_mesh_2d_tri(
        key=key,
        pts_x0=_DDATA['quad_tri']['pts_x0'],
        pts_x1=_DDATA['quad_tri']['pts_x1'],
        units='nm',
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
            bs.add_bsplines(key=km, deg=dd)

def _get_mesh(bsplines, nd=None, kind=None):
    return [
        k0 for k0, v0 in bsplines.dobj[bsplines._which_mesh].items()
        if nd in [None, v0['nd']]
        and kind in [None, v0['type']]
    ]


def _get_bs(bs, nd=None, kind=None):
    return [
        k0 for k0, v0 in bs.dobj[bs._which_bsplines].items()
        if nd in [None, bs.dobj[bs._which_mesh][v0[bs._which_mesh]]['nd']]
        and kind in [None, bs.dobj[bs._which_mesh][v0[bs._which_mesh]]['type']]
    ]


def _get_data(bs, nd=None, kind=None):

    dkd = {}
    dbs, dref = bs.get_dict_bsplines()
    for k0 in sorted(dbs.keys()):
        v0 = dbs[k0]
        for kb in v0.keys():
            km = bs.dobj[bs._which_bsplines][kb][bs._which_mesh]
            c0 = (
                nd in [None, bs.dobj[bs._which_mesh][km]['nd']]
                and kind in [None, bs.dobj[bs._which_mesh][km]['type']]
            )
            if c0 and k0 not in dkd.keys():
                dkd[k0] = {'bs': [kb], 'ref': dref[k0]}
            elif c0:
                dkd[k0]['bs'].append(kb)

    return dkd


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
            _ = bsplines.plot_mesh(
                key=km,
                ind_knot=comb[0],
                ind_cent=comb[1],
                crop=comb[2],
            )
            plt.close('all')


def _add_data_1bs_fix(bs, nd=None, kind=None, remove=None):
    lkb = _get_bs(bs, nd=nd, kind=kind)

    lkd = []
    for kb in lkb:

        lkd.append(f'{kb}_fix')
        shape = bs.dobj[bs._which_bsplines][kb]['shape']
        data = np.random.random(shape)

        bs.add_data(
            key=lkd[-1],
            data=data,
            ref=kb,
        )

        assert bs.ddata[lkd[-1]][bs._which_bsplines] == (kb,)
        bsref = bs.dobj[bs._which_bsplines][kb]['ref']
        assert bs.ddata[lkd[-1]]['ref'] == bsref

    if remove:
        for kd in lkd:
            bs.remove_data(kd)


def _add_data_1bs_arrays(bs, nd=None, kind=None, remove=None):
    lkb = _get_bs(bs, nd=nd, kind=kind)

    nt, nE = 10, 11
    if 'nt' not in bs.dref.keys():
        bs.add_ref(key='nt', size=nt)
        bs.add_ref(key='nE', size=nE)
        bs.add_data(key='t', data=np.linspace(1, 1.1, nt), ref='nt', unit='s')

    lkd = []
    for kb in lkb:

        lkd.append(f'{kb}_fix')
        shape = np.r_[nt, bs.dobj[bs._which_bsplines][kb]['shape'], nE]
        data = np.random.random(shape)

        bs.add_data(
            key=lkd[-1],
            data=data,
            ref=['nt', kb, 'nE'],
        )

        assert bs.ddata[lkd[-1]][bs._which_bsplines] == (kb,)
        bsref = bs.dobj[bs._which_bsplines][kb]['ref']
        ref = tuple(['nt'] + list(bsref) + ['nE'])
        assert bs.ddata[lkd[-1]]['ref'] == ref

    if remove:
        for kd in lkd:
            bs.remove_data(kd)


def _add_data_multibs_arrays(bs, nd=None, kind=None, remove=None):
    lkb = _get_bs(bs, nd=nd, kind=kind)

    nt, nE = 10, 11
    if 'nt' not in bs.dref.keys():
        bs.add_ref(key='nt', size=nt)
        bs.add_ref(key='nE', size=nE)
        bs.add_data(key='t', data=np.linspace(1, 1.1, nt), ref='nt', unit='s')

    lkd = []
    for ii, kb in enumerate(lkb):

        lkd.append(f'{kb}_var')

        kb2 = lkb[(ii + int(len(lkb)/2)) % len(lkb)]
        shape = np.r_[
            nt,
            bs.dobj[bs._which_bsplines][kb]['shape'],
            nE,
            bs.dobj[bs._which_bsplines][kb2]['shape'],
        ]
        data = np.random.random(shape)

        bs.add_data(
            key=lkd[-1],
            data=data,
            ref=['nt', kb, 'nE', kb2],
        )

        if bs.ddata[lkd[-1]][bs._which_bsplines] != (kb, kb2):
            msg = (
                f"Wrong '{bs._which_bsplines}' for ddata['{lkd[-1]}']:\n"
                f"{bs.ddata[lkd[-1]][bs._which_bsplines]} vs {(kb, kb2)}"
            )
            raise Exception(msg)

        bsref = bs.dobj[bs._which_bsplines][kb]['ref']
        bsref2 = bs.dobj[bs._which_bsplines][kb2]['ref']
        ref = tuple(['nt'] + list(bsref) + ['nE'] + list(bsref2))
        if bs.ddata[lkd[-1]]['ref'] != ref:
            msg = (
                f"Wrong '{bs._which_bsplines}' for ddata['{lkd[-1]}']:\n"
                f"{bs.ddata[lkd[-1]]['ref']} vs {ref}"
            )
            raise Exception(msg)

    if remove:
        for kd in lkd:
            bs.remove_data(kd)


def _interpolate(bs, nd=None, kind=None, details=None):
    dkd = _get_data(bs, nd=nd, kind=kind)

    wbs = bs._which_bsplines
    for ii, (kd, vd) in enumerate(dkd.items()):

        ref_key = vd['bs'][0]
        vect = bs.dobj[wbs][ref_key]['apex'][0]
        vect = bs.ddata[vect]['data'][::2]
        vect = np.array([vect, vect])

        # kwdargs
        shapebs = bs.dobj[bs._which_bsplines][ref_key]['shape']
        if nd == '1d':
            kwd = {'x0': vect}
        else:
            kwd = {'x0': vect, 'x1': vect}

        # crop
        if kind == 'rect' and ii % 2 == 0:
            crop = isinstance(bs.dobj[wbs][ref_key]['crop'], str)
        else:
            crop = False

        # indbs_tf
        if kind == 'rect':
            if ii%3 == 0:
                indbs_tf = bs.select_ind(
                    key=ref_key,
                    ind=([1, 2, 5], [4, 5, 3]),
                    crop=crop,
                    returnas='tuple-flat',
                )
                nbs = len(indbs_tf[0])
            else:
                indbs_tf = None
                if crop is True:
                    kcrop = bs.dobj[wbs][ref_key]['crop']
                    nbs = np.sum(bs.ddata[kcrop]['data'])
                else:
                    nbs = np.prod(shapebs)
        else:
            if ii%3 == 0:
                indbs_tf = np.arange(0, np.prod(shapebs), 2)
                nbs = indbs_tf.size
            else:
                indbs_tf = None
                nbs = np.prod(shapebs)

        kwd.update({
            'domain': None,
            'details': details,
            'indbs_tf': indbs_tf,
            'crop': crop,
        })

        # interpolate
        dout, dparam = bs.interpolate(
            keys=kd,
            ref_key=ref_key,
            return_params=True,
            debug=True,
            **kwd,
        )

        # expected shape
        if details is True:
            shape = tuple(np.r_[vect.shape, nbs].astype(int))
        else:
            shape = list(bs.ddata[kd]['shape'])
            ax0 = dparam['axis'][0]
            ax1 = dparam['axis'][-1]
            shape = tuple(np.r_[
                shape[:ax0], vect.shape, shape[ax1+1:]
            ].astype(int))

        # error msg
        if dout[kd]['data'].shape != shape:
            lstr = [
                f"\t- {k0}: {v0}"
                for k0, v0 in dparam.items()
                if k0 not in ['x0', 'x1']
            ]
            msg = (
                f"Wrong shape for '{kd}':\n"
                f"\t- expected: {shape}\n"
                f"\t- obtained: {dout[kd]['data'].shape}\n\n"
                f"\t- nd: {nd}\n"
                f"\t- kind: {kind}\n"
                f"\t- ddata['{kd}']['data'].shape: {bs.ddata[kd]['shape']}\n\n"
                + "\n".join(lstr)
            )
            raise Exception(msg)


def _bin_bs(bs, nd=None, kind=None):
    dkd = _get_data(bs, nd=nd, kind=kind)

    wbs = bs._which_bsplines
    for ii, (kd, vd) in enumerate(dkd.items()):

        if ii % 10 == 0 and len(vd['ref']) > 1:
            ref_key = vd['ref'][0]
            ax = 0
        elif len(vd['ref']) > 1:
            ref_key = vd['bs'][0]
            ax = bs.ddata[kd]['ref'].index(bs.dobj[wbs][ref_key]['ref'][0])
        else:
            continue

        vect = bs.dobj[bs._which_bsplines][vd['bs'][0]]['apex'][0]
        vect = bs.ddata[vect]['data']

        dd = np.abs(np.mean(np.diff(vect)))
        DD = vect[-1] - vect[0]
        nbins = int(DD/dd)
        bins = np.linspace(vect[0] - 0.1*DD, vect[0]+0.5*DD, nbins)

        dout = bs.binning(
            keys=kd,
            ref_key=ref_key,
            bins=bins,
        )

        shape = list(bs.ddata[kd]['shape'])
        shape[ax] = nbins - 1
        shape = tuple(shape)
        if dout[kd]['data'].shape != shape:
            shd = bs.ddata[kd]['data'].shape
            msg = (
                "Binnign of data '{kd}' along ref 'ref_key' has wrong shape:\n"
                f"\t- ddata['{kd}']['data'].shape = {shd}\n"
                f"\t- dout['{kd}']['data'].shape = {dout[kd]['data'].shape}\n"
                f"\t- expected shape = {shape}\n"
            )
            raise Exception(msg)


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
