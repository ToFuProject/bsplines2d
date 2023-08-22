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
    dunique_mesh_2d=None,
    # parameters
    val_out=None,
    nan0=None,
):
    """ Interpolate along all bsplines for multiple keys

    """

    # ----------
    # check

    keys, coll2, dres, submesh = _check(
        coll=coll,
        keys=keys,
        # sampling
        dres=dres,
        submesh=submesh,
    )

    # ---------
    # unique

    if dunique_mesh_2d is not None:
        for k0, v0 in dunique_mesh_2d.items():
            kref = f"{v0['key']}_n"
            coll2.add_ref(key=kref, size=v0['data'].size)
            coll2.add_data(ref=kref, **v0)
            dunique_mesh_2d[k0]['ref'] = kref
        dout = dunique_mesh_2d

    # ----------------
    # interpolate loop

    dbs = {}
    wm = coll._which_mesh
    wbs = coll._which_bsplines
    for ii, (k0, v0) in enumerate(dres.items()):

        # get 2d mesh
        if dunique_mesh_2d is None:

            if v0.get('x0') is None:
                dout = coll2.get_sample_mesh(
                    key=k0,
                    res=v0['res'],
                    mode=v0['mode'],
                    grid=False,
                    # store
                    store=True,
                    kx0=None,
                    kx1=None,
                )

            else:
                x0str = f'{k0}_x0_temp'
                dout = {'x0': {'key': x0str, 'ref': f'{x0str}_n'}}
                if v0.get('x1') is not None:
                    x1str = f'{k0}_x1_temp'
                    dout['x1'] = {'key': x1str, 'ref': f'{x1str}_n'}

                for k1, v1 in dout.items():
                    coll2.add_ref(key=v1['ref'], size=v0[k1].size)
                    coll2.add_data(
                        key=v1['key'],
                        data=v0[k1],
                        ref=v1['ref'],
                    )

        # compute
        for bs in v0[wbs]:

            if bs not in dbs.keys():
                dbs[bs] = {'ref': tuple([v1['ref'] for v1 in dout.values()])}

            for key in keys:

                if bs not in coll2.ddata[key][wbs]:
                    continue

                # submesh => ref_com
                kms = coll2.dobj[wbs][bs][wm]
                subkey = coll2.dobj[wm][kms]['subkey']
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
                    ref_key=bs,
                    x0=dout['x0']['key'],
                    x1=dout.get('x1', {}).get('key'),
                    submesh=submesh,
                    ref_com=ref_com,
                    grid=True,
                    details=False,
                    # parameters
                    val_out=val_out,
                    nan0=nan0,
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
                    **{k1: v1 for k1, v1 in coll2.ddata[keynew].items()},
                )

                # remove keynew
                coll2.remove_data(keynew)

    # -----------
    # clean-up

    lbs = list(set(itt.chain.from_iterable([v0[wbs] for v0 in dres.values()])))
    if wbs in coll2.dobj.keys():
        coll2.remove_bsplines(key=lbs, propagate=True)

    if wm in coll2.dobj.keys():
        coll2.remove_mesh(key=list(dres.keys()), propagate=True)

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

    lbs = list(set(itt.chain.from_iterable([coll.ddata[k0][wbs] for k0 in keys])))

    # --------------
    # rbs vs submesh

    submesh = ds._generic_check._check_var(
        submesh, 'submesh',
        types=bool,
        default=True,
    )

    dbs = {}
    for ii, bs in enumerate(lbs):
        km = coll.dobj[wbs][bs][wm]
        if submesh is True and coll.dobj[wm][km]['submesh'] is not None:
            km = coll.dobj[wm][km]['submesh']
        dbs[bs] = km

    # --------------
    # derive mesh

    lm = list(set([dbs[k0] for k0 in lbs]))

    # --------------
    # dres

    if dres is None:
        dres = {k0: {'res': None, 'mode': 'abs'} for k0 in lm}

    if isinstance(dres, (int, float)):
        dres = {k0: {'res': dres, 'mode': 'abs'} for k0 in lm}

    # safety check
    c0 = (
        isinstance(dres, dict)
        and all([kk in lm for kk in dres.keys()])
    )
    if not c0:
        msg = (
            "Arg dres must be a dict with, for each mesh\n"
            "\t- {'res': float, 'mode': str}\n"
            f"\nFor the following keys ({wm}): {lm}\n"
            f"Provided:\n{dres}\n"
        )
        raise Exception(msg)

    # loop
    for k0 in lm:
        if dres.get(k0) is None:
            dres[k0] = {'res': None, 'mode': 'abs', 'x0': None, 'x1': None}
        elif isinstance(dres[k0], dict):
            dres[k0] = {
                'res': dres[k0].get('res'),
                'mode': dres[k0].get('mode', 'abs'),
                'x0': dres[k0].get('x0'),
                'x1': dres[k0].get('x1'),
            }
        elif isinstance(dres[k0], (float, int)):
            dres[k0] = {
                'res': dres[k0],
                'mode': 'abs',
                'x0': None,
                'x1': None,
            }
        else:
            msg = (
                "Arg dres must be a dict with, for each mesh\n"
                "\t- {'res': float, 'mode': str}\n"
                f"Provided:\n{dres}"
            )
            raise Exception(msg)

        # add bsplines
        dres[k0][wbs] = [k1 for k1 in lbs if dbs[k1] == k0]

        # particular case
        if dres[k0]['x0'] is None and dres[k0]['res'] is None:
            ldeg1 = [bs for bs in dres[k0][wbs] if coll.dobj[wbs][bs]['deg'] == 1]

            if len(ldeg1) == 1:
                km = coll.dobj[wbs][ldeg1[0]][wm]
                subbs = coll.dobj[wm][km].get('subbs')
                if (
                    submesh is True
                    and subbs is not None
                    and coll.dobj[wbs][subbs]['deg'] == 1
                ):
                    ldeg1[0] = subbs
                elif submesh is False or subbs is None:
                    pass
                else:
                    ldeg1 = None

                # get default x0, x1
                if ldeg1 is not None:
                    kx = coll.dobj[wbs][ldeg1[0]]['apex']
                    dres[k0]['x0'] = coll.ddata[kx[0]]['data']
                    if len(kx) == 2:
                        dres[k0]['x1'] = coll.ddata[kx[1]]['data']

    # -------- DEBUG ------
    # lstr = [
    #     f"\t- {k0}: {v0['x0'].size if v0['x0'] is not None else None} and "
    #     f"{v0['x1'].size if v0['x1'] is not None else None}"
    #     for k0, v0 in dres.items()
    # ]
    # msg = (
    #     "The following resolutions x0 are identified:\n"
    #     + "\n".join(lstr)
    # )
    # print(msg)    # DB
    # ----- DEBUG END ------

    # ----------
    # coll2

    coll2 = coll.extract(keys=keys, vectors=True)

    return keys, coll2, dres, submesh
