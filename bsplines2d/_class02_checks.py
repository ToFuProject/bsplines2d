

import datastock as ds


# #############################################
# #############################################
#               main
# #############################################


def main(
    coll=None,
    key=None,
    key_mesh2d=None,
    key_mesh1d_angle=None,
    # optional
    **kwdargs,
):

    # --------------
    # check inputs
    # --------------

    key, key_mesh2d, key_mesh1d_angle = _check(
        coll=coll,
        key=key,
        key_mesh2d=key_mesh2d,
        key_mesh1d_angle=key_mesh1d_angle,
    )

    # --------------
    # _to_dict
    # --------------

    dref, ddata, dobj = _to_dict(
        coll=coll,
        key=key,
        key_mesh2d=key_mesh2d,
        key_mesh1d_angle=key_mesh1d_angle,
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
    key_mesh1d_angle=None,
):

    # --------------
    # key
    # --------------

    wm3d = coll._which_mesh3d
    key = ds._generic_check._obj_key(
        d0=coll.dobj.get(wm3d, {}),
        short='m3d',
        key=key,
        ndigits=2,
    )

    # also not already in mesh
    wm = coll._which_mesh
    lout = list(coll.dobj.get(wm, {}).keys())
    if key in lout:
        msg = (
            "Arg key cannot be the same as a pre-existing mesh1d or mesh2d!\n"
            f"\t- Pre-existing meshes: {lout}\n"
            f"\t- Provided: {key}\n"
        )
        raise Exception(msg)

    # --------------
    # key_mesh2d
    # --------------

    lok = [
        k0 for k0, v0 in coll.dobj.get(wm, {}).items()
        if v0['nd'] == '2d'
    ]
    key_mesh2d = ds._generic_check._check_var(
        key_mesh2d, 'key_mesh2d',
        types=str,
        allowed=lok,
        extra_msg="Must refer to a 2d mesh",
    )

    # ----------------
    # key_mesh1d_angle
    # ----------------

    lok = [
        k0 for k0, v0 in coll.dobj.get(wm, {}).items()
        if v0['nd'] == '1d'
        and str(coll.ddata[v0['knots'][0]]['units']) == 'rad'
    ]
    key_mesh1d_angle = ds._generic_check._check_var(
        key_mesh1d_angle, 'key_mesh1d_angle',
        types=str,
        allowed=lok,
        extra_msg="Must refer to a 1d mesh with units='rad'",
    )

    return key, key_mesh2d, key_mesh1d_angle


# ######################################
# ######################################
#               to_dict
# ######################################


def _to_dict(
    coll=None,
    key=None,
    key_mesh2d=None,
    key_mesh1d_angle=None,
    # attributes
    **kwdargs,
):

    # -------------
    # prepare data

    latt = ['dim', 'quant', 'name', 'units']

    wm = coll._which_mesh
    mtype = (
        coll.dobj[wm][key_mesh2d]['type'],
        coll.dobj[wm][key_mesh1d_angle]['type'],
    )

    # -------------
    # prepare dref, ddata

    dref = None
    ddata = None

    # -------------
    # prepare dict

    # dobj
    wm3d = coll._which_mesh3d
    dobj = {
        wm3d: {
            key: {
                'nd': '3d',
                'type': mtype,
                'mesh2d': key_mesh2d,
                'mesh1d_angle': key_mesh1d_angle,
            },
        },
    }

    # additional attributes
    for k0, v0 in kwdargs.items():
        if k0 not in latt:
            dobj[wm3d][key][k0] = v0

    return dref, ddata, dobj
