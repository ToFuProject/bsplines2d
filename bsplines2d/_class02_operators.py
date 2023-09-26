# -*- coding: utf-8 -*-


# Built-in


# Common
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
    keycropped = coll.dobj['bsplines'][key]['ref-bs'][0]
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
    ref0 = coll.dobj[wbs][key]['ref-bs']

    if deg > 0:
        kbsm1 = f'{keym}_bs{deg-1}'
        if kbsm1 in coll.dobj[wbs].keys():
            rm1 = coll.dobj[wbs][kbsm1]['ref-bs']
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
