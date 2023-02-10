# -*- coding: utf-8 -*-


# Built-in


# Common
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
    store=None,
    returnas=None,
    # specific to deg = 0
    centered=None,
    # to return gradR, gradZ, for D1N2 deg 0, for tomotok
    returnas_element=None,
):

    # -------
    # compute

    (
        opmat, operator, geometry, dim, ref, crop,
        store, returnas, key,
    ) = _get_bsplines_operator(
        coll=coll,
        key=key,
        operator=operator,
        geometry=geometry,
        crop=crop,
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
        if operator == 'D1':
            name = f'{key}-{operator}-dx0'
            if crop is True:
                name = f'{name}-cropped'
            coll.add_data(
                key=name,
                data=opmat[0],
                ref=ref,
                units='',
                name=operator,
                dim=dim,
            )
            name = f'{key}-{operator}-dx1'
            if crop is True:
                name = f'{name}-cropped'
            coll.add_data(
                key=name,
                data=opmat[1],
                ref=ref,
                units='',
                name=operator,
                dim=dim,
            )

        elif operator in ['D0N1', 'D0N2']:
            name = f'{key}-{operator}-{geometry}'
            if crop is True:
                name = f'{name}-cropped'
            coll.add_data(
                key=name,
                data=opmat,
                ref=ref,
                units='',
                name=operator,
                dim=dim,
            )
        elif operator == 'D1N2':
            name = f'{key}-{operator}-dx0-{geometry}'
            if crop is True:
                name = f'{name}-cropped'
            coll.add_data(
                key=name,
                data=opmat[0],
                ref=ref,
                units='',
                name=operator,
                dim=dim,
            )
            name = f'{key}-{operator}-dx1-{geometry}'
            if crop is True:
                name = f'{name}-cropped'
            coll.add_data(
                key=name,
                data=opmat[1],
                ref=ref,
                units='',
                name=operator,
                dim=dim,
            )
        elif operator == 'D2N2':
            name = f'{key}-{operator}-d2x0-{geometry}'
            if crop is True:
                name = f'{name}-cropped'
            coll.add_data(
                key=name,
                data=opmat[0],
                ref=ref,
                units='',
                name=operator,
                dim=dim,
            )
            name = f'{key}-{operator}-d2x1-{geometry}'
            if crop is True:
                name = f'{name}-cropped'
            coll.add_data(
                key=name,
                data=opmat[1],
                ref=ref,
                units='',
                name=operator,
                dim=dim,
            )
            name = f'{key}-{operator}-dx0x1-{geometry}'
            if crop is True:
                name = f'{name}-cropped'
            coll.add_data(
                key=name,
                data=opmat[2],
                ref=ref,
                units='',
                name=operator,
                dim=dim,
            )
        else:
            msg = "Unknown opmat type!"
            raise Exception(msg)

    # return
    if returnas is True:
        return opmat, operator, geometry, dim, ref, crop


# #############################################################################
# #############################################################################
#                           subroutine
# #############################################################################


def _get_bsplines_operator(
    coll=None,
    key=None,
    operator=None,
    geometry=None,
    crop=None,
    store=None,
    returnas=None,
    # specific to deg = 0
    centered=None,
    # to return gradR, gradZ, for D1N2 deg 0, for tomotok
    returnas_element=None,
):

    # check inputs
    lk = list(coll.dobj.get('bsplines', {}).keys())
    key = ds._generic_check._check_var(
        key, 'key',
        types=str,
        allowed=lk,
    )

    store = ds._generic_check._check_var(
        store, 'store',
        default=True,
        types=bool,
    )

    returnas = ds._generic_check._check_var(
        returnas, 'returnas',
        default=store is False,
        types=bool,
    )

    crop = ds._generic_check._check_var(
        crop, 'crop',
        default=True,
        types=bool,
    )

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

    # compute and return
    (
        opmat, operator, geometry, dim,
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

    # cropping
    if operator == 'D1':
        ref = (keycropped, keycropped)
    elif operator == 'D0N1':
        ref = (keycropped,)
    elif 'N2' in operator:
        ref = (keycropped, keycropped)

    return opmat, operator, geometry, dim, ref, crop, store, returnas, key
