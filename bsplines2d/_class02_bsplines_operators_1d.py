# -*- coding: utf-8 -*-


# Built-in


# Common
import numpy as np
import scipy.sparse as scpsp
import datastock as ds


from . import _utils_bsplines_operators as _operators


_LOPERATORS_INT = [
    'D1',
    'D2',
    'D3',
    'D0N1',
    'D0N2',
    'D1N2',
    'D2N2',
    'D3N2',
]


# ###############################################################
# ###############################################################
#                   Mesh2DRect - bsplines - operators
# ###############################################################


def get_operators(
    operator=None,
    geometry=None,
    deg=None,
    knots_mult=None,
    knots_per_bs=None,
    overlap=None,
    sparse_fmt=None,
    # specific to deg = 0
    centered=None,
    # for D1N2 only (to feed to tomotok / Mfr)
    returnas_element=None,
):

    # ------------
    # check inputs

    (
        operator, geometry, sparse_fmt, returnas_element,
    ) = _check(
        deg=deg,
        operator=operator,
        geometry=geometry,
        sparse_fmt=sparse_fmt,
        returnas_element=returnas_element,
    )

    # ------------
    # prepare

    nbs = knots_per_bs.shape[1]
    kn = knots_per_bs
    shape = (nbs, nbs)
    indbs = np.arange(0, nbs)

    if 'N2' in operator and deg >= 1:
        # get intersection indices array
        nbtot = np.sum(overlap >= 0)

        # prepare data and indices arrays
        if operator == 'D0N2':
            data = np.full((nbtot,), np.nan)
            row = np.zeros((nbtot,), dtype=int)
            column = np.zeros((nbtot,), dtype=int)
        elif operator == 'D1N2':
            datadx = np.full((nbtot,), np.nan)
            row = np.zeros((nbtot,), dtype=int)
            column = np.zeros((nbtot,), dtype=int)
        elif operator == 'D2N2':
            datad2x = np.full((nbtot,), np.nan)
            row = np.zeros((nbtot,), dtype=int)
            column = np.zeros((nbtot,), dtype=int)

    # ------------
    # D0 - integral

    geom = geometry
    if operator == 'D0N1':

        if deg == 0:

            opmat = _operators._D0N1_Deg0(kn, geom)

        elif deg == 1:

            opmat = _operators._D0N1_Deg1(kn, geom)

        elif deg == 2:

            in1, in21, in22, in3 = _operators._D0N1_Deg2(kn, geom)
            opmat = (in1 + in21 + in22 + in3)

        elif deg == 3:

            msg = "Integral D0N1 not implemented for deg=3 yet!"
            raise NotImplementedError(msg)

    # ------------
    # D1 - gradient

    elif operator == 'D1':

        # Treat separately discrete case
        if deg == 0:
            grad = _operators._D1_Deg0(
                knots_mult=knots_mult,
                nbs=nbs,
                centered=centered,
            )

        elif deg == 1:
            grad = _operators._D1_Deg1(
                knots_mult=knots_mult,
                nbs=nbs,
            )

        elif deg == 2:
            grad = _operators._D1_Deg2(
                knots_mult=knots_mult,
                nbs=nbs,
            )

        elif deg >= 3:
            raise NotImplementedError()

        opmat = (scpsp.csc_matrix(grad),)

    # ------------
    # D0N2

    elif operator == 'D0N2' and deg == 0:

        ix = _operators._D0N2_Deg0(kn, geom)
        opmat = scpsp.diags(
            [ix],
            [0],
            shape=None,
            format=sparse_fmt,
            dtype=float,
        )

    elif operator == 'D0N2':

        # pre-compute integrals
        if deg == 1:
            ix = _operators._D0N2_Deg1(knots_mult, geom)

        elif deg == 2:
            ix = _operators._D0N2_Deg2(knots_mult, geom)

        elif deg == 3:
            msg = "Integral D0N2 not implemented for deg=3!"
            raise NotImplementedError(msg)

        # set non-diagonal elements
        i0 = 0
        for ii in range(nbs):

            # general case
            overlapi = overlap[:, ii][overlap[:, ii] > ii]

            # diagonal element
            data[i0] = ix[0, ii]
            row[i0] = indbs[ii]
            column[i0] = indbs[ii]
            i0 += 1

            # non-diagonal elements (symmetric)
            for jj in overlapi:

                # store (i, j) and (j, i) (symmetric matrix)
                iix = ix[jj - ii, ii]

                data[i0:i0+2] = iix
                row[i0:i0+2] = (indbs[ii], indbs[jj])
                column[i0:i0+2] = (indbs[jj], indbs[ii])
                i0 += 2

        assert i0 == nbtot
        opmat = scpsp.csc_matrix(
            (data, (row, column)),
            shape=shape,
        )

    # ------------
    # D1N2

    elif operator == 'D1N2':

        # Treat separately discrete case
        if deg == 0:
            grad = _operators._D1_Deg0(
                knots_mult=knots_mult,
                nbs=nbs,
            )

            # surface elements
            if geom == 'linear':
                dx = (knots_mult[1:] - knots_mult[:-1])
            else:
                dx = 0.5*(knots_mult[1:]**2 - knots_mult[:-1]**2)

            # Does not seem to give positive definite matrix, TBC?
            if returnas_element is True:
                opmat = (scpsp.csc_matrix(
                    grad*np.sqrt(dx[:, None])
                ),)

            else:
                opmat = (scpsp.csc_matrix(
                    (grad.T.dot(grad))*(dx[:, None])
                ),)

        else:

            # pre-compute integrals for exact operator deg >= 1
            if deg == 1:
                idx = _operators._D1N2_Deg1(knots_mult, geom)

            elif deg == 2:
                idx = _operators._D1N2_Deg2(knots_mult, geom)

            elif deg == 3:
                msg = "Integral D1N2 not implemented for deg=3!"
                raise NotImplementedError(msg)

            # set non-diagonal elements
            i0 = 0
            for ii in range(nbs):

                # general case
                overlapi = overlap[:, ii][overlap[:, ii] > ii]

                # diagonal element
                datadx[i0] = idx[0, ii]
                row[i0] = indbs[ii]
                column[i0] = indbs[ii]
                i0 += 1

                # non-diagonal elements (symmetric)
                for jj in overlapi:

                    # store (i, j) and (j, i) (symmetric matrix)
                    iidx = idx[jj - ii, ii]

                    datadx[i0:i0+2] = iidx
                    row[i0:i0+2] = (indbs[ii], indbs[jj])
                    column[i0:i0+2] = (indbs[jj], indbs[ii])
                    i0 += 2

            assert i0 == nbtot
            opmat = (scpsp.csc_matrix(
                (datadx, (row, column)),
                shape=shape,
            ),)

    # ------------
    # D2N2

    elif operator == 'D2N2':

        # pre-compute integrals
        if deg in [0, 1]:
            msg = f"degree {deg} too low for operator {operator}"
            raise Exception(msg)
        if deg == 2:
            id2x = _operators._D2N2_Deg2(knots_mult, geom)
            idx = _operators._D1N2_Deg2(knots_mult, geom)
            ix = _operators._D0N2_Deg2(knots_mult, geom)
        elif deg == 3:
            msg = "Integral D2N2 not implemented for deg=3!"
            raise NotImplementedError(msg)

        # set non-diagonal elements
        i0 = 0
        for ii in range(nbs):

                # general case
                overlapi = overlap[:, ii][overlap[:, ii] > ii]

                # diagonal element
                datad2x[i0] = id2x[0, ii]
                row[i0] = indbs[ii]
                column[i0] = indbs[ii]
                i0 += 1

                # non-diagonal elements (symmetric)
                for jj in overlapi:

                    # store (i, j) and (j, i) (symmetric matrix)
                    iid2x = id2x[jj - ii, ii]

                    datad2x[i0:i0+2] = iid2x
                    row[i0:i0+2] = (indbs[ii], indbs[jj])
                    column[i0:i0+2] = (indbs[jj], indbs[ii])
                    i0 += 2

        assert i0 == nbtot
        opmat = (scpsp.csc_matrix(
            (datad2x, (row, column)),
            shape=shape,
        ),)

    # ------------
    # D3N2

    elif operator == 'D3N2' and deg == 3:

        raise NotImplementedError("Integral D3N2 not implemented for deg=3!")

    return opmat, operator, geometry


# ###############################################################
# ###############################################################
#                   checks
# ###############################################################


def _check(
    deg=None,
    operator=None,
    geometry=None,
    sparse_fmt=None,
    returnas_element=None,
):

    # deg
    deg = ds._generic_check._check_var(
        deg, 'deg',
        types=int,
        allowed=[0, 1, 2, 3],
    )

    # operator
    operator = ds._generic_check._check_var(
        operator, 'operator',
        default='D0N1',
        types=str,
        allowed=_LOPERATORS_INT,
    )

    # geometry
    geometry = ds._generic_check._check_var(
        geometry, 'geometry',
        default='toroidal',
        types=str,
        allowed=['toroidal', 'linear'],
    )

    # sparse_fmt
    sparse_fmt = ds._generic_check._check_var(
        sparse_fmt, 'sparse_fmt',
        default='csc',
        types=str,
        allowed=['dia', 'csr', 'csc', 'lil'],
    )

    # returnas_element
    lok = [False]
    if operator == 'D1N2' and deg == 0:
        lok.append(True)
    returnas_element = ds._generic_check._check_var(
        returnas_element, 'returnas_element',
        default=False,
        types=bool,
        allowed=lok,
    )

    return operator, geometry, sparse_fmt, returnas_element
