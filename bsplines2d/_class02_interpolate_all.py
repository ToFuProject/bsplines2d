# -*- coding: utf-8 -*-


# Built-in
import copy
import itertools as itt


# Common
import datastock as ds

# local


# ###########################################################
# ###########################################################
#               Main routine
# ###########################################################


def interpolate_all_bsplines(
    coll=None,
    keys=None,
    # sampling
    dres=None,
    submesh=None,
):

    # ----------
    # check

    keys, coll2, dres, submesh = _check(
        coll=coll,
        keys=keys,
        # sampling
        dres=dres,
        submesh=submesh,
    )

    # ----------------
    # interpolate loop

    wm = coll._which_mesh
    wbs = coll._which_bsplines
    dbs = {}
    for ii, (k0, v0) in enumerate(dres.items()):

        # get 2d mesh
        dout = None
        if v0.get('x0') is None:
            dout = coll2.get_sample_mesh(
                key=v0[wm],
                res=v0['res'],
                mode=v0['mode'],
                grid=False,
                # store
                store=True,
                kx0=None,
                kx1=None,
            )
        else:
            x0str = f'{v0[wm]}_x0_temp'
            dout = {'x0': {'key': x0str, 'ref': f'{x0str}_n'}}
            if v0.get('x1') is not None:
                x1str = f'{v0[wm]}_x1_temp'
                dout = {'xi1': {'key': x1str, 'ref': f'{x1str}_n'}}
                for k1, v1 in dout.items():
                    coll.add_ref(key=v1['ref'], size=v0[k1].size)
                    coll.add_data(
                        key=v1['key'],
                        data=v1[k1],
                        ref=v1['ref'],
                    )

        # compute
        for key in keys:

            if k0 in coll2.ddata[key][wbs]:

                # submesh => ref_com
                km = coll2.dobj[wbs][k0][wm]
                subkey = coll2.dobj[wm][km]['subkey']
                if submesh is True and subkey is not None:
                    refsub = coll2.ddata[subkey[0]]['ref']
                    ref = coll2.ddata[key]['ref']
                    if refsub[0] in ref:
                        ref_com = refsub[0]
                    elif refsub[-1] in ref:
                        ref_com = refsub[-1]
                    else:
                        ref_com = None
                else:
                    ref_com = None

                # interpolate
                coll2 = coll2.interpolate(
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
                    inplace=True,
                )

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

        # fill dbs
        dbs[k0] = {
            'ref': tuple([v1['ref'] for v1 in dout.values()]),
        }

    # -----------
    # clean-up

    if wbs in coll2.dobj.keys():
        coll2.remove_bsplines(key=list(dres.keys()), propagate=True)

    if wm in coll2.dobj.keys():
        coll2.remove_mesh(key=[v0[wm] for v0 in dres.values()], propagate=True)

    return coll2, dbs


# ###########################################################
# ###########################################################
#               check
# ###########################################################


def _check(
    coll=None,
    keys=None,
    # sampling
    knots0=None,
    knots1=None,
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

    if isinstance(keys, str):
        keys = [keys]
    keys = ds._generic_check._check_var_iter(
        keys, 'keys',
        types=(list, tuple),
        types_iter=str,
        allowed=lok,
    )

    lbs = list(itt.chain.from_iterable([coll.ddata[k0][wbs] for k0 in keys]))

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
        if dres.get(bs) is None:
            dres[bs] = {'res': None, 'mode': 'abs'}
        elif isinstance(dres[bs], dict):
            dres[bs] = {
                'res': dres[bs].get('res'),
                'mode': dres[bs].get('mode', 'abs'),
                'x0': dres[bs].get('x0'),
                'x1': dres[bs].get('x1'),
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


    # ----------
    # coll2

    coll2 = coll.extract(keys=keys, vectors=True)

    return keys, coll2, dres, submesh
