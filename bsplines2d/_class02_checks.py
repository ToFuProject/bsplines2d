

from . import _generic_check


# #############################################
# #############################################
#               main
# #############################################


def main(
    coll=None,
    key=None,
    key_mesh2d=None,
    knots_phi=None,
    # optional
    **kwdargs,
):

    # --------------
    # check inputs
    # --------------

    key, key_mesh2d, knots_phi = _check(
        coll=coll,
        key=key,
        key_mesh2d=key_mesh2d,
        knots_phi=knots_phi,
    )

    # --------------
    # _to_dict
    # --------------

    dref, ddata, dobj = _to_dict(
        coll=coll,
        key=key,
        key_mesh2d=key_mesh2d,
        key_mesh_angle=key_mesh_angle,
        # attributes
        **kwdargs,
    )

    return key, dref, ddata, dobj


# #############################################
# #############################################
#           Check inputs
# #############################################


def _check(
    coll=None,
    key=None,
    key_mesh2d=None,
    knots_phi=None,
):

    # --------------
    # key
    # --------------

    key = _generic_check._obj_key(
        d0=coll.dobj.get(wm3d, {}),
        short='m3d',
        key=key,
        ndigits=2,
    )

    # --------------
    # key_mesh2d
    # --------------

    wm = coll._which_mesh
    lok = list(coll.dobj.get(wm, {}).keys())
    key = _generic_check._check_var(
        key_mesh2d, 'key_mesh2d',
        types=str,
        allowed=lok,
    )

    # --------------
    # knots_phi
    # --------------

    return key, key_mesh2d, knots_phi


# ######################################
# ######################################
#               to_dict
# ######################################


def _to_dict(
    coll=None,
    key=None,
    key_mesh2d=None,
    key_mesh_angle=None,
    # attributes
    **kwdargs,
):

    # ---------
    # prepare

    # attributes
    latt = ['dim', 'quant', 'name', 'units']
    dim, quant, name, units = [kwdargs.get(ss) for ss in latt]

    # -------------
    # prepare dict

    # dref
    dref = {
        kkr: {
            'size': knots.size,
        },
        kcr: {
            'size': cents.size,
        },
    }

    # ddata
    ddata = {
        kk: {
            'data': knots,
            'units': units,
            # 'source': None,
            'dim': dim,
            'quant': quant,
            'name': name,
            'ref': kkr,
        },
        kc: {
            'data': cents,
            'units': units,
            # 'source': None,
            'dim': dim,
            'quant': quant,
            'name': name,
            'ref': kcr,
        },
    }

    # dobj
    wm3d = coll._which_mesh3d
    dobj = {
        wm3d: {
            key: {
                'nd': '3d',
                'type': mtype,
                'mesh2d': key_mesh2d,
                'knots_phi': kphi
            },
        },
    }

    # additional attributes
    for k0, v0 in kwdargs.items():
        if k0 not in latt:
            dobj[wm3d][key][k0] = v0

    return dref, ddata, dobj
