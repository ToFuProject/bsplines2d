# -*- coding: utf-8 -*-


# Built-in
import itertools as itt


# Common
import numpy as np
import astropy.units as asunits
import datastock as ds


# specific


# #############################################################################
# #############################################################################
#                           get operators
# #############################################################################


def get_bsplines_operator(
    coll=None,
    key=None,
    operator=None,
    geometry=None,
    crop=None,
    # store vs return
    store=None,
    returnas=None,
    return_param=None,
    # specific to deg = 0
    centered=None,
    # to return gradR, gradZ, for D1N2 deg 0, for tomotok
    returnas_element=None,
):

    # ---------
    # check

    (
        key, store, returnas,
        crop, cropbs, cropbs_flat, keycropped,
    ) = _check(
        coll=coll,
        key=key,
        operator=operator,
        geometry=geometry,
        crop=crop,
        store=store,
        returnas=returnas,
    )

    # -------
    # compute

    dout = _get_bsplines_operator(
        coll=coll,
        key=key,
        operator=operator,
        geometry=geometry,
        crop=crop,
        cropbs=cropbs,
        cropbs_flat=cropbs_flat,
        keycropped=keycropped,
        store=store,
        returnas=returnas,
        # specific to deg = 0
        centered=centered,
        # to return gradR, gradZ, for D1N2 deg 0, for tomotok
        returnas_element=returnas_element,
    )

    # ------
    # store

    if store is True:

        for k0, v0 in dout.items():

            if operator == 'D1' and None in v0['ref']:
                continue

            coll.add_data(**v0)

    # return
    if returnas is True:

        if return_param:
            dpar = {
                'key': key,
                'keys': list(dout.keys()),
                'operator': operator,
                'geometry': geometry,
                'crop': crop,
            }
            return dout, dpar

        return dout


# ################################################################
# ################################################################
#               check
# ################################################################


def _check(
    coll=None,
    key=None,
    operator=None,
    geometry=None,
    crop=None,
    store=None,
    returnas=None,
):

    # --------
    # key

    wm = coll._which_mesh
    wbs = coll._which_bsplines
    lk = list(coll.dobj.get(wbs, {}).keys())
    key = ds._generic_check._check_var(
        key, 'key',
        types=str,
        allowed=lk,
    )
    keym = coll.dobj[wbs][key][wm]
    nd = coll.dobj[wm][keym]['nd']
    mtype = coll.dobj[wm][keym]['type']

    # --------
    # store

    store = ds._generic_check._check_var(
        store, 'store',
        default=True,
        types=bool,
    )

    # --------
    # returnas

    returnas = ds._generic_check._check_var(
        returnas, 'returnas',
        default=store is False,
        types=bool,
    )

    # --------
    # crop

    crop = ds._generic_check._check_var(
        crop, 'crop',
        default=True,
        types=bool,
    )
    if not (nd == '2d' and mtype == 'rect'):
        crop = False

    # cropbs
    cropbs = coll.dobj['bsplines'][key]['crop']
    keycropped = coll.dobj['bsplines'][key]['ref_bs'][0]
    if cropbs not in [None, False] and crop is True:
        cropbs_flat = coll.ddata[cropbs]['data'].ravel(order='F')
        if coll.dobj['bsplines'][key]['deg'] == 0:
            cropbs = coll.ddata[cropbs]['data']
        keycropped = f"{keycropped}-crop"
    else:
        cropbs = False
        cropbs_flat = False

    return key, store, returnas, crop, cropbs, cropbs_flat, keycropped


# ################################################################
# ################################################################
#                       compute
# ################################################################


def _get_bsplines_operator(
    coll=None,
    key=None,
    operator=None,
    geometry=None,
    crop=None,
    cropbs=None,
    cropbs_flat=None,
    keycropped=None,
    store=None,
    returnas=None,
    # specific to deg = 0
    centered=None,
    # to return gradR, gradZ, for D1N2 deg 0, for tomotok
    returnas_element=None,
):

    # -------------------
    # compute and return

    (
        opmat, operator, geometry,
    ) = coll.dobj['bsplines'][key]['class'].get_operator(
        operator=operator,
        geometry=geometry,
        cropbs_flat=cropbs_flat,
        # specific to deg=0
        cropbs=cropbs,
        centered=centered,
        # to return gradR, gradZ, for D1N2 deg 0, for tomotok
        returnas_element=returnas_element,
    )

    # -----------
    # format dout

    wm = coll._which_mesh
    wbs = coll._which_bsplines
    keym = coll.dobj[wbs][key][wm]
    nd = coll.dobj[wm][keym]['nd']
    mtype = coll.dobj[wm][keym]['type']
    deg = coll.dobj[wbs][key]['deg']

    dout = _dout(
        coll=coll,
        key=key,
        opmat=opmat,
        operator=operator,
        geometry=geometry,
        keycropped=keycropped,
        crop=crop,
        nd=nd,
        mtype=mtype,
        deg=deg,
    )

    return dout


# ###################################################
# ###################################################
#                   ref
# ###################################################


def _dout(
    coll=None,
    key=None,
    opmat=None,
    operator=None,
    geometry=None,
    keycropped=None,
    crop=None,
    nd=None,
    mtype=None,
    deg=None,
):

    # --------
    # get refs

    ref, units = _ref_units(
        coll=coll,
        key=key,
        opmat=opmat,
        operator=operator,
        geometry=geometry,
        keycropped=keycropped,
        nd=nd,
        mtype=mtype,
        deg=deg,
    )

    geom = geometry[:3]

    # ----------
    # matrix types

    nnd = int(nd[0])
    if operator == 'D1':
        kmat = [f'M{ii}' for ii in range(nnd)]
    if operator in ['D0N1']:
        kmat = 'M'
    elif operator in ['D0N2']:
        kmat = 'tMM'
    elif operator in ['D1N2']:
        kmat = [f'tMM{ii}' for ii in range(nnd)]
    elif operator in ['D2N2']:
        lcomb = [] if nd == '1d' else [(0,1)]
        kmat = (
            [f'tMM{ii}{ii}' for ii in range(nnd)]
            + [f'tMM{ii}{jj}' for ii, jj in lcomb]
        )

    # ----------
    # build dout

    dout = {}
    if operator in ['D0N1', 'D0N2']:

        k0 = f'{key}_{operator}_{geom}'
        if crop is True:
            k0 = f'{k0}_crop'

        dout[kmat] = {
            'key': k0,
            'data': opmat,
            'ref': ref,
            'units': units[0],
            'dim': operator,
        }

    elif operator in ['D1', 'D1N2']:

        for ii in range(nnd):

            k0 = f'{key}_{operator}_d{ii}'
            if 'N' in operator:
                k0 = f'{k0}_{geom}'
            if crop is True:
                k0 = f'{k0}_crop'

            dout[kmat[ii]] = {
                'key': k0,
                'data': opmat[ii],
                'ref': ref,
                'units': units[ii],
                'dim': operator,
            }

    elif operator in ['D2N2']:

        for ii, kk in enumerate(kmat):

            k0 = f'{key}_{operator}_d{kk[-2:]}_{geom}'
            if crop is True:
                k0 = f'{k0}_crop'

            dout[kk] = {
                'key': k0,
                'data': opmat[ii],
                'ref': ref,
                'units': units[ii],
                'dim': operator,
            }

    return dout


def _ref_units(
    coll=None,
    key=None,
    opmat=None,
    operator=None,
    geometry=None,
    keycropped=None,
    nd=None,
    mtype=None,
    deg=None,
):

    # --------
    # prepare

    wm = coll._which_mesh
    wbs = coll._which_bsplines
    keym = coll.dobj[wbs][key][wm]

    ref = keycropped
    ref0 = coll.dobj[wbs][key]['ref_bs']

    if deg > 0:
        kbsm1 = f'{keym}_bs{deg-1}'
        if kbsm1 in coll.dobj[wbs].keys():
            rm1 = coll.dobj[wbs][kbsm1]['ref_bs']
            if ref0 != ref:
                rm1 = f'{rm1}_crop'
        else:
            rm1 = None
    else:
        rm1 = ref

    # ----------
    # ref

    if operator == 'D1':
        ref = (rm1, ref)

    elif operator == 'D0N1':
        ref = (ref,)

    elif 'N2' in operator:
        ref = (ref, ref)

    # --------
    # units

    apex = coll.dobj[wbs][key]['apex']
    u0 = coll.ddata[apex[0]]['units']
    if nd == '1d':
        units = [_units(u0, operator, geometry)]

    else:
        u1 = coll.ddata[apex[1]]['units']
        units0 = _units(u0, operator, geometry)
        units1 = _units(u1, operator, 'linear')

        if operator in ['D0N1', 'D0N2']:
            units = [units0 * units1]

        elif operator in ['D1', 'D1N2']:
            units = [units0, units1]

        elif operator in ['D2N2']:
            units = [units0, units0*u0*units1*u1, units1]

    return ref, units


def _units(u0=None, operator=None, geometry=None):

    if str(u0) == '-':
        u0 = asunits.Unit('')
    elif isinstance(u0, str):
        u0 = asunits.Unit(u0)

    if operator == 'D1':
        units = asunits.Unit(1/u0)

    elif operator == 'D0N1':
        if geometry == 'linear':
            units = u0
        else:
            units = u0**2

    elif operator == 'D0N2':
        if geometry == 'linear':
            units = u0
        else:
            units = u0**2

    elif operator == 'D1N2':
        if geometry == 'linear':
            units = asunits.Unit(1/u0)
        else:
            units = asunits.Unit('')

    elif operator == 'D2N2':
        if geometry == 'linear':
            units = asunits.Unit(1/u0)**3
        else:
            units = asunits.Unit(1/u0)**2

    return units


# ###################################################
# ###################################################
#                   Apply
# ###################################################


def apply_operator(
    coll=None,
    # parameters
    key=None,
    keybs=None,
    # operator
    operator=None,
    geometry=None,
    # store
    store=None,
    key_store=None,
    # return
    returnas=None,
):
    """ Apply an operator on desired quantity

    key can be a list of data with the same bsplines

    Notes:
        * Only works for D0N1 so far
        * Works with and w/o cropbs
        * Does not work with subkey so far

    Computes data, ref and units
    Returns as dict
    Optionally store result

    """

    # -------------
    # check inputs
    # -------------

    (
        key, keybs,
        operator,
        store, key_store,
        returnas,
    ) = _apply_operator_check(**locals())

    # -------------
    # get operator
    # -------------

    integ_op = coll.add_bsplines_operator(
        key=keybs,
        operator=operator,
        geometry=geometry,
        store=False,
        returnas=True,
    )

    # -----------------
    # prepare output
    # -----------------

    ddata, daxis, cropbs = _apply_operator_prepare(
        coll=coll,
        key=key,
        keybs=keybs,
        # for units
        operator=operator,
    )

    # -----------------
    # apply
    # -----------------

    if operator == 'D0N1':
        ind = [-1]
        if cropbs is None:
            for k0 in key:
                ddata[k0]['data'][...] = np.tensordot(
                    integ_op['M']['data'],
                    coll.ddata[k0]['data'],
                    (ind, daxis[k0]['axis']),
                )
        else:
            for k0 in key:
                ddata[k0]['data'][...] = np.tensordot(
                    integ_op['M']['data'],
                    coll.ddata[k0]['data'][daxis[k0]['slice']],
                    (ind, daxis[k0]['axis']),
                )
    else:
        raise NotImplementedError()

    # -------------
    # store
    # -------------

    if store is True:
        for k0, v0 in ddata.items():
            coll.add_data(**v0)

    # -------------
    # return
    # -------------

    if returnas is True:
        return ddata


def _apply_operator_check(
    coll=None,
    # parameters
    key=None,
    keybs=None,
    # operator
    operator=None,
    # store
    store=None,
    key_store=None,
    # return
    returnas=None,
    # unused
    **kwdargs,
):

    # --------------
    # keys
    # --------------

    # ------
    # key

    wbs = coll._which_bsplines
    lok = [k0 for k0, v0 in coll.ddata.items() if v0.get(wbs) is not None]
    if isinstance(key, str):
        key = [key]
    if key is None and len(lok) == 1:
        key = lok
    key = ds._generic_check._check_var_iter(
        key, 'key',
        types=list,
        types_iter=str,
        allowed=lok,
    )

    # -----------
    # prepare check on keybs

    lbs = set(itt.chain.from_iterable([coll.ddata[k0][wbs] for k0 in key]))

    # check all have at least one in common
    lbs_all = [
        k0 for k0 in lbs
        if all([k0 in coll.ddata[k1][wbs] for k1 in key])
    ]
    if len(lbs_all) == 0:
        lstr = [
            f"\t- coll.ddata['{k0}']['{wbs}'] = {coll.ddata[k0][wbs]}"
            for k0 in key
        ]
        msg = (
            "All provided keys must have at least one '{wbs}' in common!\n"
            + "\n".join(lstr)
        )
        raise Exception(msg)

    # ------------
    # keybs

    lok = lbs_all
    keybs = ds._generic_check._check_var(
        keybs, 'keybs',
        types=str,
        allowed=lok,
    )

    # -------------------
    # operator
    # --------------------

    lok = ['D0N1']
    operator = ds._generic_check._check_var(
        operator, 'operator',
        types=str,
        allowed=lok,
    )

    # -------------------
    # store vs return
    # --------------------

    # store
    store = ds._generic_check._check_var(
        store, 'store',
        types=bool,
        default=False,
    )

    # returnas
    returnas = ds._generic_check._check_var(
        returnas, 'returnas',
        types=bool,
        default=not store,
    )

    # -------------------
    # key_store
    # --------------------

    # key_store
    if store is True:

        if key_store is None:
            key_store = [f"{k0}_{operator}" for k0 in key]

        if isinstance(key_store, str):
            key_store = [key_store]

        lout = list(coll.ddata.keys())
        key_store = ds._generic_check._check_var_iter(
            key_store, 'key_store',
            types=list,
            types_iter=str,
            excluded=lout,
        )

    else:
        key_store = None

    return (
        key, keybs,
        operator,
        store, key_store,
        returnas,
    )


def _apply_operator_prepare(
    coll=None,
    key=None,
    keybs=None,
    operator=None,
):

    # ----------
    # cropbs

    wbs = coll._which_bsplines
    cropbs = coll.dobj[wbs][keybs].get('crop')
    if isinstance(cropbs, str) and (cropbs in coll.ddata.keys()):
        cropbs = coll.ddata[cropbs]['data']
    else:
        cropbs = None

    # ----------
    # refbs

    refbs = coll.dobj[wbs][keybs]['ref']

    # ----------
    # unitsbs

    unitsbs = [
        coll.ddata[k0]['units']
        for k0 in coll.dobj[wbs][keybs]['apex']
    ]

    # --------------
    # fill dict

    ddata, daxis = {}, {}
    for k0 in key:

        # ref
        ref0 = coll.ddata[k0]['ref']
        ref = tuple([rr for rr in ref0 if rr not in refbs])

        # axis
        axisf = tuple([ii for ii, rr in enumerate(ref0) if rr in refbs])

        # safety check
        assert np.allclose(axisf, axisf[0] + np.arange(len(refbs)))

        # cropbs => reduce to a single dimension
        if len(axisf) > 1 and cropbs is not None:
            axis = (axisf[0],)
        else:
            axis = axisf

        # shape
        shape0 = coll.ddata[k0]['data'].shape
        shape = [ss for ii, ss in enumerate(shape0) if ref0[ii] not in refbs]

        # units
        units0 = asunits.Unit(coll.ddata[k0]['units'])
        if operator == 'D0N1':
            units = units0
            for uu in unitsbs:
                units = units * uu
        else:
            raise NotImplementedError()

         # populate
        ddata[k0] = {
            'data': np.full(shape, np.nan),
            'ref': ref,
            'units': units,
        }

        # slicing
        if cropbs is None:
            sli = None
        else:
            sli = tuple([
                cropbs if ii == axis[0]
                else slice(None)
                for ii in range(len(shape0))
                if ii not in axisf[1:]
            ])

        daxis[k0] = {
            'slice': sli,
            'axis': axis,
        }

    return ddata, daxis, cropbs