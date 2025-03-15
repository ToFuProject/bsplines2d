# -*- coding: utf-8 -*-


# Common
import numpy as np
import datastock as ds


from . import _generic_mesh


_ELEMENTS = 'knots'


# ################################################################
# ###############################################################
#                           mesh generic check
# ###############################################################

"""
def _mesh2D_polar_check(
    coll=None,
    radius=None,
    angle=None,
    radius2d=None,
    angle2d=None,
    key=None,
    # parameters
    radius_dim=None,
    radius_quant=None,
    radius_name=None,
    radius_units=None,
    angle_dim=None,
    angle_quant=None,
    angle_name=None,
):

    # key
    key = ds._generic_check._check_var(
        key, 'key',
        types=str,
        excluded=list(coll.dobj.get('mesh', {}).keys())
    )

    # --------------------
    # check / format input

    krk, krc, kkr, kcr = _generic_mesh.names_knots_cents(key=key, knots_name='r')
    kak, kac, kka, kca = _generic_mesh.names_knots_cents(key=key, knots_name='ang')

    # radius data
    radius, _, _ = _mesh1D_check(
        x=radius,
        x_name='radius',
        uniform=False,
    )

    # angle data
    c0 = (
        angle is None
        or (
            hasattr(angle, '__iter__')
            and np.asarray(angle).ndim == 1
            and np.unique(angle).size == np.array(angle).size
            and np.allclose(
                np.unique(np.arctan2(np.sin(angle), np.cos(angle))),
                angle,
            )
        )
    )
    if not c0:
        msg = (
            "Arg angle either\n:"
            "\t- None: radial-only polar mesh"
            "\t- convertible to a 1d increasing array\n"
            "\t\t it must be in radians\n"
            "\t\t it must be in the [-pi; pi] interval\n"
            f"\t- Provided: {angle}"
        )
        raise Exception(msg)

    # extract data
    rknot = np.unique(radius)
    rcent = 0.5*(rknot[1:] + rknot[:-1])

    # radius2d
    dradius = _check_polar_2dquant(
        coll=coll,
        quant2d=radius2d,
        quant2d_name='radius2d',
        dim=radius_dim,
        quant=radius_quant,
        name=radius_name,
        units=radius_units,
    )

    if callable(radius2d):
        keysm = None
    else:
        keysm = coll.dobj['bsplines'][coll.ddata[radius2d]['bsplines']]['mesh']

    if angle is not None:
        aknot = np.unique(np.arctan2(np.sin(angle), np.cos(angle)))
        acent = 0.5*(aknot[1:] + aknot[:-1])
        amid = 0.5*(aknot[-1] + (2.*np.pi + aknot[0]))
        amid = np.arctan2(np.sin(amid), np.cos(amid))
        if amid < acent[0]:
            acent = np.r_[amid, acent]
        else:
            acent = np.r_[acent, amid]

    # -------
    # angle2d

    if angle2d is not None:
        dangle = _check_polar_2dquant(
            coll=coll,
            quant2d=angle2d,
            quant2d_name='angle2d',
            dim=angle_dim,
            quant=angle_quant,
            name=angle_name,
            units='rad',
        )

        # check angle units = rad
        if dangle['units'] != 'rad':
            msg = (
                "Angle units must be rad\n"
                f"\t Provided: {dangle['units']}"
            )
            raise Exception(msg)

        # check angle2d is like radius2d
        c0 = (
            (callable(radius2d) and callable(angle2d))
            or coll._ddata[radius2d]['ref'] == coll._ddata[angle2d]['ref']
        )
        if not c0:
            msg = (
                "radius2d and angle2d must be of the same type, either:\n"
                "\t- both callable\n"
                "\t- both data keys with identical ref!\n"
                f"Provided:\n"
                f"\t- radius2d: {radius2d}\n"
                f"\t- angle2d: {angle2d}\n"
            )
            raise Exception(msg)

    # --------------------
    # prepare dict

    # dref
    dref = {
        krk: {'size': rknot.size},
        krc: {'size': rcent.size},
    }

    if angle is not None:
        dref.update({
            kak: {
                'size': aknot.size,
            },
            kac: {
                'size': acent.size,
            },
        })

    # ddata
    ddata = {
        kkr: {
            'data': rknot,
            'ref': krk,
            **dradius,
        },
        kcr: {
            'data': rcent,
            'ref': krc,
            **dradius,
        },
    }

    if angle is not None:
        ddata.update({
            kka: {
                'data': aknot,
                'ref': kak,
                **dangle,
            },
            kca: {
                'data': acent,
                'ref': kac,
                **dangle,
            },
        })

    # dobj
    if angle is None:
        dmesh = {
            key: {
                'type': 'polar',
                'knots': (kkr,),
                'cents': (kcr,),
                'shape_c': rcent.shape,
                'shape_k': rknot.shape,
                'radius2d': radius2d,
                'angle2d': angle2d,
                'submesh': keysm,
                'crop': False,
            },
        }
    else:
        dmesh = {
            key: {
                'type': 'polar',
                'knots': (kkr, kka),
                'cents': (kcr, kca),
                'shape_c': (rcent.size, acent.size),
                'shape_k': (rknot.size, aknot.size),
                'radius2d': radius2d,
                'angle2d': angle2d,
                'submesh': keysm,
                'crop': False,
            },
        }

    return dref, ddata, dmesh


def _check_polar_2dquant(
    quant2d=None,
    coll=None,
    quant2d_name=None,
    # parameters
    dim=None,
    quant=None,
    name=None,
    units=None,
):

    if coll.dobj.get('bsplines') is not None:
        lok = [
            k0 for k0, v0 in coll.ddata.items()
            if v0.get('bsplines') in coll.dobj['bsplines'].keys()
        ]
    else:
        lok = []

    lc = [
        callable(quant2d),
        isinstance(quant2d, str) and quant2d in lok
    ]
    if not any(lc):
        msg = (
            f"Arg {quant2d_name} must be either:\n"
            f"\t- callable: {quant2d_name} = func(R, Z)\n"
            f"\t- key to existing 2d data in {lok}\n"
            f"Provided: {quant2d}\n"
        )
        raise Exception(msg)

    # quantities
    dquant = {'dim': dim, 'quant': quant, 'name': name, 'units': units}
    if isinstance(quant2d, str):
        for k0 in dquant.keys():
            if dquant[k0] is None:
                dquant[k0] = str(coll.ddata[quant2d][k0])

    return dquant
"""