# -*- coding: utf-8 -*-


# Common
import numpy as np
import datastock as ds


from . import _generic_mesh
from . import _class01_checks_1d as _checks_1d
    

# #############################################################################
# #############################################################################
#                           mesh generic check
# #############################################################################


def check(
    coll=None,
    key=None,
    # knots
    knots0=None,
    knots1=None,
    knots0_name=None,
    knots1_name=None,
    # automated
    domain=None,
    res=None,
):

    # --------
    # keys
    
    # key
    key = ds._generic_check._obj_key(
        d0=coll._dobj.get(coll._which_mesh, {}),
        short='m',
        key=key,   
    )
    
    # -------------
    # knots vectors
    
    knots0, knots1, res0, res1 = _mesh2DRect_check(
        knots0=knots0,
        knots1=knots1,
        # automated
        domain=domain,
        res=res,
    )
    cents0 = 0.5*(knots0[1:] + knots0[:-1])
    cents1 = 0.5*(knots1[1:] + knots1[:-1])

    # -------------- 
    # to dict

    dref, ddata, dobj = _to_dict(
        coll=coll,
        key=key,
        knots0=knots0,
        knots1=knots1,
        cents0=cents0,
        cents1=cents1,
        knots0_name=knots0_name,
        knots1_name=knots1_name,
    )

    return key, dref, ddata, dobj


# #############################################################################
# #############################################################################
#                           check knots
# #############################################################################



def _mesh2DRect_check(
    knots0=None,
    knots1=None,
    knots0_name=None,
    knots1_name=None,
    domain=None,
    res=None,
):

    # --------------
    # check inputs

    # (domain, res) vs (knots0, knots1)
    lc = [
        domain is not None,
        knots0 is not None and knots1 is not None,
    ]
    if all(lc) or not any(lc):
        msg = (
            "Please provide (domain, res) xor (knots0, knots1), not both:\n"
            "Provided:\n"
            f"\t- domain, res: {domain}, {res}\n"
            f"\t- knots0, knots1: {knots0}, {knots1}\n"
        )
        raise Exception(msg)

    if lc[0]:
        # domain
        c0 = (
            isinstance(domain, list)
            and len(domain) == 2
            and all([
                hasattr(dd, '__iter__') and len(dd) >= 2 for dd in domain
            ])
        )
        if not c0:
            msg = (
                "Arg domain must be a list of 2 iterables of len() >= 2\n"
                f"Provided: {domain}"
            )
            raise Exception(msg)

        # res
        c0 = (
            res is None
            or np.isscalar(res)
            or isinstance(res, list) and len(res) == 2
        )
        if not c0:
            msg = (
                "Arg res must be a int, float or array or a list of 2 such\n"
                f"Provided: {res}"
            )
            raise Exception(msg)

        if np.isscalar(res) or res is None:
            res = [res, res]

        # -------------
        # derive knots

        knots0, res0, _ = _mesh2DRect_X_check(domain[0], res=res[0])
        knots1, res1, _ = _mesh2DRect_X_check(domain[1], res=res[1])

    elif lc[1]:
        
        # knots0
        knots0, res0 = _checks_1d._check_knots(
            knots=knots0,
            knots_name=knots0_name,
            uniform=True,        
        )

        # knots1
        knots1, res1 = _checks_1d._check_knots(
            knots=knots1,
            knots_name=knots1_name,
            uniform=True, 
        )

    return knots0, knots1, res0, res1


def _mesh2DRect_X_check(
    x=None,
    res=None,
):
    """ Returns knots (x) and associated resolution

    res can be:
        - int: numbr of mesh elements desired between knots
        - float: desired average mesh element size
        - array of floats: (one for each x, desired approximate mesh size)

    """

    # ------------
    # Check inputs

    # x
    try:
        x = np.unique(np.ravel(x).astype(float))
    except Exception as err:
        msg = f"x must be convertible to a sorted flat array of floats!\n{err}"
        raise Exception(msg)

    # res
    if res is None:
        res = 10

    lc = [
        isinstance(res, (int, np.int64, np.int32)) and len(x) == 2,
        isinstance(res, (float, np.floating)) and len(x) == 2,
        isinstance(res, (list, tuple, np.ndarray)) and len(x) == len(res),
    ]
    if not any(lc):
        msg = (
            "Arg res must be:\n"
            "\t- int: nb of mesh elements along x\n"
            "\t       requires len(x) = 2\n"
            "\t- float: approximate desired mesh element size along x\n"
            "\t       requires len(x) = 2\n"
            "\t- iterable: approximate desired mesh element sizes along x\n"
            "\t       requires len(x) = len(res)\n"
        )
        raise Exception(msg)

    if lc[0]:
        x_new = np.linspace(x[0], x[1], int(res)+1)
        res_new = res
        indsep = None

    elif lc[1]:
        nb = int(np.ceil((x[1]-x[0]) / res))
        x_new = np.linspace(x[0], x[1], nb+1)
        res_new = np.mean(np.diff(x))
        indsep = None

    else:

        # check conformity
        res = np.ravel(res).astype(float)
        delta = np.diff(x)
        res_sum = res[:-1] + res[1:]
        ind = res_sum > delta + 1.e-14
        if np.any(ind):
            msg = (
                "Desired resolution is not achievable for the following:\n"
                f"res_sum: {res_sum[ind]}\n"
                f"delta  : {delta[ind]}"
            )
            raise Exception(msg)

        # compute nn
        # nn = how many pairs can fit in the interval
        npairs = np.round(delta/res_sum).astype(int)
        res_sum_new = delta / npairs

        fract = res[:-1] / res_sum

        res_new = [None for ii in range(len(x)-1)]
        x_new = [None for ii in range(len(x)-1)]
        for ii in range(len(x)-1):
            res_new[ii] = (
                res_sum_new[ii]
                * np.linspace(fract[ii], 1.-fract[ii], 2*npairs[ii])
            )
            if ii == 0:
                res_add = np.concatenate(([0], np.cumsum(res_new[ii])))
            else:
                res_add = np.cumsum(res_new[ii])
            x_new[ii] = x[ii] + res_add

        indsep = np.cumsum(npairs[:-1]*2)
        res_new = np.concatenate(res_new)
        x_new = np.concatenate(x_new)

    return x_new, res_new, indsep


# #############################################################################
# #############################################################################
#                           to_dict
# #############################################################################


def _to_dict(
    coll=None,
    key=None,
    knots0=None,
    knots1=None,
    cents0=None,
    cents1=None,
    knots0_name=None,
    knots1_name=None,
    res0=None,
    res1=None,
    # attributes
    **kwdargs,
):

    # ---------
    # check
    
    knots0_name = ds._generic_check._check_var(
        knots0_name, 'knots0_name',
        types=str,
        default='x',
    )
    
    knots1_name = ds._generic_check._check_var(
        knots1_name, 'knots1_name',
        types=str,
        default='y',
        excluded=[knots0_name],
    )
    
    # --------------------
    # check / format input

    kkr0, kcr0, kk0, kc0 = _generic_mesh.names_knots_cents(
        key=key,
        knots_name=knots0_name,
    )
    kkr1, kcr1, kk1, kc1 = _generic_mesh.names_knots_cents(
        key=key,
        knots_name=knots1_name,
    )
    
    # attributes
    latt = ['dim', 'quant', 'name', 'units']
    dim, quant, name, units = [kwdargs.get(ss) for ss in latt]

    variable = not (np.isscalar(res0) and np.isscalar(res1))

    # --------------------
    # prepare dict

    # dref
    dref = {
        kkr0: {
            'size': knots0.size,
        },
        kkr1: {
            'size': knots1.size,
        },
        kcr0: {
            'size': cents0.size,
        },
        kcr1: {
            'size': cents1.size,
        },
    }

    # ddata
    ddata = {
        kk0: {
            'data': knots0,
            'units': units,
            # 'source': None,
            'dim': dim,
            'quant': quant,
            'name': name,
            'ref': kkr0,
        },
        kk1: {
            'data': knots1,
            'units': units,
            # 'source': None,
            'dim': dim,
            'quant': quant,
            'name': name,
            'ref': kkr1,
        },
        kc0: {
            'data': cents0,
            'units': units,
            # 'source': None,
            'dim': dim,
            'quant': quant,
            'name': name,
            'ref': kcr0,
        },
        kc1: {
            'data': cents1,
            'units': units,
            # 'source': None,
            'dim': dim,
            'quant': quant,
            'name': name,
            'ref': kcr1,
        },
    }

    # dobj
    dobj = {
        coll._which_mesh: {
            key: {
                'nd': '2d',
                'type': 'rect',
                'knots': (kk0, kk1),
                'cents': (kc0, kc1),
                # 'ref-k': (kRk, kZk),
                # 'ref-c': (kRc, kZc),
                'shape-c': (cents0.size, cents1.size),
                'shape-k': (knots0.size, knots1.size),
                'variable': variable,
                'crop': False,
            },
        }
    }
    
    # additional attributes
    for k0, v0 in kwdargs.items():
        if k0 not in latt:
            dobj[coll._which_mesh][key][k0] = v0
    
    return dref, ddata, dobj


def _from_croppoly(crop_poly=None, domain=None):

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
                poly = np.concatenate((poly, poly[:, 0:1]), axis=1)

        # -------------
        # Get domain from poly

        if domain is None:
            domain = [
                [poly[0, :].min(), poly[0, :].max()],
                [poly[1, :].min(), poly[1, :].max()],
            ]

    return domain, poly