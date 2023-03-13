# -*- coding: utf-8 -*-


# Built-in
import copy


# Common
import datastock as ds

# local



# ###########################################################
# ###########################################################
#               Main routine
# ###########################################################


def interpolate_all_bsplines(
    coll=None,
    key=None,
    dres=None,
    submesh=None,
):

    # ----------
    # check

    dres, submesh = _check(
        coll=coll,
        key=key,
        dres=dres,
        submesh=submesh,
    )

    # ----------------
    # interpolate loop

    coli = coll
    inplace = False
    wm = coll._which_mesh
    wbs = coll._which_bsplines
    dbs = {}
    for ii, (k0, v0) in enumerate(dres.items()):

        # get 2d mesh
        dout = coli.get_sample_mesh(
            key=v0[wm],
            res=v0['res'],
            mode=v0['mode'],
            grid=False,
            # store
            store=True,
            kx0=None,
            kx1=None,
        )

        # submesh => ref_com
        km = coll.dobj[wbs][k0][wm]
        subkey = coll.dobj[wm][km]['subkey']
        if submesh is True and subkey is not None:
            refsub = coll.ddata[subkey[0]]['ref']
            ref = coll.ddata[key]['ref']
            if refsub[0] in ref:
                ref_com = refsub[0]
            elif refsub[-1] in ref:
                ref_com = refsub[-1]
            else:
                ref_com = None
        else:
            ref_com = None

        # compute
        coll2 = coli.interpolate(
            keys=key,
            ref_key=k0,
            x0=dout['x0']['key'],
            x1=dout.get('x1', {}).get('key'),
            submesh=submesh,
            ref_com=ref_com,
            grid=True,
            details=False,
            # return vs store
            returnas=object,
            return_params=False,
            store=True,
            inplace=inplace,
        )

        # remove / replace
        if ii == 0:
            for k0, v0 in dout.items():
                coli.remove_ref(v0['ref'])


        # remove original data
        if key in coll2.ddata.keys():
            coll2.remove_data(key)

        # rename new data
        keynew = f'{key}_interp'
        coll2.add_data(
            key=key,
            **{k0: v0 for k0, v0 in coll2.ddata[keynew].items()},
        )

        # remove keynew
        coll2.remove_data(keynew)

        # udpate params
        if ii == 0:
            coli = coll2
            inplace = True

        # fill dbs
        dbs[k0] = {
            'ref': tuple([v1['ref'] for v1 in dout.values()]),
        }

    # -----------
    # clean-up

    if wm in coll2.dobj.keys():
        del coll2._dobj[wm]
    if wbs in coll2.dobj.keys():
        del coll2._dobj[wbs]

    for k0, v0 in coll2.ddata.items():
        if wbs in v0.keys():
            del coll2._ddata[k0][wbs]

    return coll2, dbs


# ###########################################################
# ###########################################################
#               check
# ###########################################################


def _check(
    coll=None,
    key=None,
    dres=None,
    submesh=None,
):

    # ----
    # key

    wm = coll._which_mesh
    wbs = coll._which_bsplines

    lok = [
        k0 for k0, v0 in coll.ddata.items()
        if isinstance(v0[wbs], tuple)
    ]

    key = ds._generic_check._check_var(
        key, 'key',
        types=str,
        allowed=lok,
    )

    lbs = coll.ddata[key][wbs]

    # --------------
    # dres

    if dres is None:
        dres = {bs: {'res': None, 'mode': 'abs'} for bs in lbs}

    if isinstance(dres, (int, float)):
        dres = {bs: {'res': dres, 'mode': 'abs'} for bs in lbs}

    # safety check
    c0 = (
        isinstance(dres, dict)
        and all([kk in lbs for kk in dres.keys()])
    )
    if not c0:
        msg = (
            "Arg dres must be a dict with, for each bsplines\n"
            "\t- {'res': float, 'mode': str}\n"
            f"\nFor the following keys ({wbs}): {lbs}\n"
            f"Provided:\n{dres}"
        )
        raise Exception(msg)

    # loop
    for bs in lbs:
        if bs not in dres.keys():
            dres[bs] = {'res': None, 'mode': 'abs'}
        elif isinstance(dres[bs], dict):
            dres[bs] = {
                'res': dres[bs].get('res'),
                'mode': dres[bs].get('mode', 'abs'),
            }
        elif isinstance(dres[bs], (float, int)):
            dres[bs] = {
                'res': dres[bs],
                'mode': 'abs',
            }
        else:
            msg = (
                "Arg dres must be a dict with, for each bsplines\n"
                "\t- {'res': float, 'mode': str}\n"
                f"Provided:\n{dres}"
            )
            raise Exception(msg)

    # --------------
    # rbs vs submesh

    submesh = ds._generic_check._check_var(
        submesh, 'submesh',
        types=bool,
        default=True,
    )

    for ii, bs in enumerate(lbs):
        km = coll.dobj[wbs][bs][wm]
        if submesh and coll.dobj[wm][km]['submesh'] is not None:
            dres[bs][wm] = coll.dobj[wm][km]['submesh']
        else:
            dres[bs][wm] = km

    return dres, submesh
