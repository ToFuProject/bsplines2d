# -*- coding: utf-8 -*-


# Built-in
import warnings


# Common
import numpy as np
import datastock as ds


# ##############################################################
# ##############################################################
#                          add data on mesh / bsplines
# ###############################################################


def add_data_meshbsplines_ref(
    coll=None,
    ref=None,
    data=None,
):

    dmesh = coll._dobj.get(coll._which_mesh)
    dbsplines = coll._dobj.get(coll._which_bsplines)

    if dmesh is None or dbsplines is None:
        return ref, data

    # ref is str
    if isinstance(ref, str):
        ref = [ref]

    # ref is tuple
    if isinstance(ref, (tuple, list)):

        # ref contains mesh
        rm = [(ii, rr) for ii, rr in enumerate(ref) if rr in dmesh.keys()]
        if len(rm) > 0:

            ref = list(ref)
            for (ki, km) in rm:
                kbs = [
                    k0 for k0, v0 in dbsplines.items()
                    if v0[coll._which_mesh] == km
                ]
                if len(kbs) == 1:
                    ref[ki] = kbs[0]
                elif len(kbs) > 1:
                    msg = (
                        "ref contains reference to mesh with several bsplines!\n"
                        f"\t- ref: {ref}\n"
                        f"\t- mesh bsplines: {kbs}\n"
                    )
                    raise Exception(msg)

        # ref contains bsplines
        rbs = [(ii, rr) for ii, rr in enumerate(ref) if rr in dbsplines.keys()]
        while len(rbs) > 0:

            ii, kb = rbs[0]

            ref = np.r_[
                ref[:ii],
                dbsplines[kb]['ref'],
                ref[ii+1:],
            ]

            rbs = [(ii, rr) for ii, rr in enumerate(ref) if rr in dbsplines.keys()]

            # repeat data if taken from ntri > 1
            data = _repeat_data_ntri(
                ref=ref,
                rbs1=kb,
                refbs=dbsplines[kb]['ref'],
                data=data,
                # mesh
                km=dbsplines[kb][coll._which_mesh],
                dmesh=dmesh,
                dbsplines=dbsplines,
            )

        ref = tuple(ref)

    return ref, data


def _repeat_data_ntri(
    ref=None,
    rbs1=None,
    refbs=None,
    data=None,
    # mesh
    km=None,
    dmesh=None,
    dbsplines=None,
):
    """ If triangular mesh with ntri > 1 => repeat data """

    c0 = (
        dmesh[km]['type'] == 'tri'
        and dmesh[km]['ntri'] > 1
    )
    if c0:
        ntri = dmesh[km]['ntri']
        indr = ref.tolist().index(refbs[0])
        nbs = dbsplines[rbs1]['shape'][0]
        ndata = data.shape[indr]
        if ndata == nbs:
            pass
        elif ndata == nbs / ntri:
            data = np.repeat(data, ntri, axis=indr)
        else:
            msg = (
                "Mismatching data shape vs multi-triangular mesh:\n"
                f"\t- data.shape[tribs] = {ndata}\n"
                f"\t- expected {nbs} / {ntri} = {nbs / ntri}\n"
            )
            raise Exception(msg)

    return data


def _set_data_bsplines(coll=None):

    if coll.dobj.get(coll._which_bsplines) is not None:

        wbs = coll._which_bsplines
        for k0, v0 in coll._ddata.items():

            lbs = [
                (v0['ref'].index(v1['ref'][0]), k1)
                for k1, v1 in coll.dobj[wbs].items()
                if v1['ref'] == tuple([
                    rr for rr in v0['ref']
                    if rr in v1['ref']
                ])
            ]

            if len(lbs) == 0:
                pass
            else:
                # re-order
                ind = np.argsort([bb[0] for bb in lbs])
                lbs = [lbs[ii][1] for ii in ind]

                # store
                coll._ddata[k0]['bsplines'] = tuple(lbs)


# ###############################################################
# ###############################################################
#                           Mesh2DRect - bsplines
# ###############################################################


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


# ###############################################################
# ###############################################################
#                       Remove bsplines
# ###############################################################


def remove_bsplines(coll=None, key=None, propagate=None):

    # ----------
    # check

    # key
    wbs = coll._which_bsplines
    if wbs not in coll.dobj.keys():
        return

    if isinstance(key, str):
        key = [key]
    key = ds._generic_check._check_var_iter(
        key, 'key',
        types=(list, tuple),
        types_iter=str,
        allowed=coll.dobj.get(wbs, {}).keys(),
    )

    # propagate
    propagate = ds._generic_check._check_var(
        propagate, 'propagate',
        types=bool,
        default=True,
    )

    # ---------
    # remove

    for k0 in key:

        # specific data
        ldata = list(set((
            list(coll.dobj[wbs][k0]['apex'])
            + [
                k1 for k1, v1 in coll.ddata.items()
                if k0 in v1[wbs]
            ]
        )))
        for dd in ldata:
            coll.remove_data(dd, propagate=propagate)

        # specific ref
        lref = (
            coll.dobj[wbs][k0]['ref']
            + coll.dobj[wbs][k0]['ref-bs']
        )
        for rr in lref:
            if rr in coll.dref.keys():
                coll.remove_ref(rr, propagate=propagate)

        # obj
        coll.remove_obj(which=wbs, key=k0, propagate=propagate)
