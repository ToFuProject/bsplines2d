# -*- coding: utf-8 -*-


# Built-in
import warnings


# Common
import numpy as np
from matplotlib.tri import Triangulation as mplTri
import datastock as ds
    

# #############################################################################
# #############################################################################
#                          add data on mesh / bsplines
# #############################################################################


def add_data_meshbsplines_ref(
    coll=None,
    ref=None,
    data=None,
    # ressources
    which_mesh=None,
    which_bsplines=None,
):

    dmesh = coll._dobj.get(which_mesh)
    dbsplines = coll._dobj.get(which_bsplines)
    
    if dmesh is None or dbsplines is None:
        return ref, data

    # ref is str
    if isinstance(ref, str):
        ref = [ref]

    # ref is tuple
    if isinstance(ref, (tuple, list)):

        # ref contains mesh
        rm = [(ii, rr) for ii, rr in enumerate(ref) if rr in dmesh.keys()]
        if len(rm) > 1:
            msg = (
                "ref contains references to several meshes!\n"
                f"\t- ref: {ref}\n"
                f"\t- meshes: {rm}\n"
            )
            raise Exception(msg)

        elif len(rm) == 1:
            ref = list(ref)
            kbs = [
                k0 for k0, v0 in dbsplines.items()
                if v0[which_mesh] == rm[0][1]
            ]
            if len(kbs) == 1:
                ref[rm[0][0]] = kbs[0]
            elif len(kbs) > 1:
                msg = (
                    "ref contains reference to mesh with several bsplines!\n"
                    f"\t- ref: {ref}\n"
                    f"\t- mesh bsplines: {kbs}\n"
                )
                raise Exception(msg)

        # ref contains bsplines
        rbs = [(ii, rr) for ii, rr in enumerate(ref) if rr in dbsplines.keys()]
        if len(rbs) > 1:
            msg = (
                "ref contains references to several bsplines!"
                f"\t- ref: {ref}\n"
                f"\t- splines: {rbs}\n"
            )
            raise Exception(msg)

        elif len(rbs) == 1:
            ref = np.r_[
                ref[:rbs[0][0]],
                dbsplines[rbs[0][1]]['ref'],
                ref[rbs[0][0]+1:],
            ]

            # repeat data if taken from ntri > 1 
            data = _repeat_data_ntri(
                ref=ref,
                rbs1=rbs[0][1],
                refbs=dbsplines[rbs[0][1]]['ref'],
                data=data,
                # mesh
                km=dbsplines[rbs[0][1]][which_mesh],
                dmesh=dmesh,
                dbsplines=dbsplines,
            )

    return tuple(ref), data
    kxk, kxc = f'{key}-{x_name}-nk', f'{key}-{x_name}-nc'
    kkx, kcx = f'{key}-k-{x_name}', f'{key}-c-{x_name}'
    return kxk, kxc, kkx, kcx

    # ------------
    # check inputs

    c0 = hasattr(crop_poly, '__iter__') and len(crop_poly) == 2
    lc = [
        crop_poly is None,
        (
            c0
            and isinstance(crop_poly, tuple)
            and crop_poly[0].__class__.__name__ == 'Config'
            and (isinstance(crop_poly[1], str) or crop_poly[1] is None)
        )
        or crop_poly.__class__.__name__ == 'Config',
        c0
        and all([
            hasattr(cc, '__iter__') and len(cc) == len(crop_poly[0])
            for cc in crop_poly[1:]
        ])
        and np.asarray(crop_poly).ndim == 2
    ]

    if not any(lc):
        msg = (
            "Arg config must be a Config instance!"
        )
        raise Exception(msg)

    # -------------
    # Get polyand domain

    if lc[0]:
        # trivial case
        poly = None

    else:

        # -------------
        # Get poly from input

        if lc[1]:
            # (config, structure name)

            if crop_poly.__class__.__name__ == 'Config':
                config = crop_poly
                key_struct = None
            else:
                config, key_struct = crop_poly

            # key_struct if None
            if key_struct is None:
                lk, ls = zip(*[
                    (ss.Id.Name, ss.dgeom['Surf']) for ss in config.lStructIn
                ])
                key_struct = lk[np.argmin(ls)]

            # poly
            poly = config.dStruct['dObj']['Ves'][key_struct].Poly_closed

        else:

            # make sure poly is np.ndarraya and closed
            poly = np.asarray(crop_poly).astype(float)
            if not np.allclose(poly[:, 0], poly[:, -1]):
                poly = np.concatenate((poly, poly[:, 0:1]))

        # -------------
        # Get domain from poly

        if domain is None:
            domain = [
                [poly[0, :].min(), poly[0, :].max()],
                [poly[1, :].min(), poly[1, :].max()],
            ]

    return domain, poly


# #############################################################################
# #############################################################################
#                           Mesh2DRect - bsplines
# #############################################################################


def _mesh_bsplines(key=None, lkeys=None, deg=None):

    # key
    key = ds._generic_check._check_var(
        key, 'key',
        types=str,
        allowed=lkeys,
    )

    # deg
    deg = ds._generic_check._check_var(
        deg, 'deg',
        types=int,
        default=2,
        allowed=[0, 1, 2, 3],
    )

    # keybs
    keybs = f'{key}_bs{deg}'

    return key, keybs, deg
