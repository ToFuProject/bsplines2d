# -*- coding: utf-8 -*-


# Built-in


# Common
import numpy as np
import scipy.sparse as scpsp
import datastock as ds


from . import _utils_bsplines_operators as _operators
from . import _class02_bsplines_operators_1d as _operators_1d


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


def get_mesh2dRect_operators(
    operator=None,
    geometry=None,
    deg=None,
    knotsx_mult=None,
    knotsy_mult=None,
    knotsx_per_bs=None,
    knotsy_per_bs=None,
    overlap=None,
    sparse_fmt=None,
    cropbs_flat=None,
    # specific to deg = 0
    cropbs=None,
    centered=None,
    # for D1N2 only (to feed to tomotok / Mfr)
    returnas_element=None,
):

    # ------------
    # check inputs

    (
        operator, geometry, sparse_fmt, returnas_element,
    ) = _operators_1d._check(
        deg=deg,
        operator=operator,
        geometry=geometry,
        sparse_fmt=sparse_fmt,
        returnas_element=returnas_element,
    )

    # ------------
    # prepare

    nx, ny = knotsx_per_bs.shape[1], knotsy_per_bs.shape[1]
    kR = np.tile(knotsx_per_bs, ny)
    kZ = np.repeat(knotsy_per_bs, nx, axis=1)
    nbs = nx*ny

    if cropbs_flat is None:
        cropbs_flat = False
    if cropbs_flat is not False:
        c0 = (
            isinstance(cropbs_flat, np.ndarray)
            and cropbs_flat.shape == (nbs,)
            and cropbs_flat.dtype == np.bool_
        )
        if not c0:
            msg = (
                f"Arg cropbs_flat must be a bool array of shape {(nbs,)}\n"
                f"Provided: {cropbs_flat.shape}"
            )
            raise Exception(msg)
        nbscrop = cropbs_flat.sum()
        shape = (nbscrop, nbscrop)
        indbs = -np.ones((nbs,), dtype=int)
        indbs[cropbs_flat] = np.arange(0, nbscrop)
    else:
        shape = (nbs, nbs)
        indbs = np.arange(0, nbs)

    if 'N2' in operator and deg >= 1:
        # get intersection indices array
        if cropbs_flat is False:
            nbtot = np.sum(overlap >= 0)
        else:
            ind = cropbs_flat[None, :] & cropbs_flat[overlap]
            nbtot = np.sum(ind)

        # prepare data and indices arrays
        if operator == 'D0N2':
            data = np.full((nbtot,), np.nan)
            row = np.zeros((nbtot,), dtype=int)
            column = np.zeros((nbtot,), dtype=int)
        elif operator == 'D1N2':
            datadR = np.full((nbtot,), np.nan)
            datadZ = np.full((nbtot,), np.nan)
            row = np.zeros((nbtot,), dtype=int)
            column = np.zeros((nbtot,), dtype=int)
        elif operator == 'D2N2':
            datad2R = np.full((nbtot,), np.nan)
            datad2Z = np.full((nbtot,), np.nan)
            datadRZ = np.full((nbtot,), np.nan)
            row = np.zeros((nbtot,), dtype=int)
            column = np.zeros((nbtot,), dtype=int)

    # ------------
    # D0 - integral

    geom = geometry
    if operator == 'D0N1':

        if deg == 0:

            opZ = _operators._D0N1_Deg0(kZ)
            opR = _operators._D0N1_Deg0(kR, geom)
            opmat = opR * opZ

        elif deg == 1:

            opZ = _operators._D0N1_Deg1(kZ)
            opR = _operators._D0N1_Deg1(kR, geom)
            opmat = opR * opZ

        elif deg == 2:

            iZ1, iZ21, iZ22, iZ3 = _operators._D0N1_Deg2(kZ)
            iR1, iR21, iR22, iR3 = _operators._D0N1_Deg2(kR, geom)
            opmat = (
                (iR1 + iR21 + iR22 + iR3)
                * (iZ1 + iZ21 + iZ22 + iZ3)
            )

        elif deg == 3:

            msg = "Integral D0N1 not implemented for deg=3 yet!"
            raise NotImplementedError(msg)

        # crop
        if cropbs_flat is not False:
            opmat = opmat[cropbs_flat]

    # ------------
    # D1 - gradient

    elif operator == 'D1':

        # Treat separately discrete case
        if deg == 0:
            gradR, gradZ = _operators._D1_Deg0_2d(
                knotsx_mult=knotsx_mult,
                knotsy_mult=knotsy_mult,
                cropbs=cropbs,
                cropbs_flat=cropbs_flat,
                nx=nx,
                ny=ny,
                nbs=nbs,
                centered=centered,
            )
            opmat = (
                scpsp.csc_matrix(gradR),
                scpsp.csc_matrix(gradZ),
            )

        elif deg >= 1:
            raise NotImplementedError()

    # ------------
    # D0N2

    elif operator == 'D0N2' and deg == 0:

        iZ = _operators._D0N2_Deg0(kZ, geom)
        iR = _operators._D0N2_Deg0(kR, geom)

        if cropbs_flat is not False:
            iR = iR[cropbs_flat]
            iZ = iZ[cropbs_flat]

        opmat = scpsp.diags(
            [iR*iZ],
            [0],
            shape=None,
            format=sparse_fmt,
            dtype=float,
        )

    elif operator == 'D0N2':

        # pre-compute integrals
        if deg == 1:
            iR = _operators._D0N2_Deg1(knotsx_mult, geom)
            iZ = _operators._D0N2_Deg1(knotsy_mult)

        elif deg == 2:
            iR = _operators._D0N2_Deg2(knotsx_mult, geom)
            iZ = _operators._D0N2_Deg2(knotsy_mult)

        elif deg == 3:
            msg = "Integral D0N2 not implemented for deg=3!"
            raise NotImplementedError(msg)

        # set non-diagonal elements
        i0 = 0
        for ir in range(nx):
            for iz in range(ny):

                iflat = ir + iz*nx
                if cropbs_flat is not False and not cropbs_flat[iflat]:
                    continue

                # general case
                overlapi = overlap[:, iflat][overlap[:, iflat] > iflat]

                # diagonal element
                data[i0] = iR[0, ir] * iZ[0, iz]
                row[i0] = indbs[iflat]
                column[i0] = indbs[iflat]
                i0 += 1

                # non-diagonal elements (symmetric)
                for jflat in overlapi:

                    if cropbs_flat is not False and not cropbs_flat[jflat]:
                        continue

                    jr = jflat % nx
                    jz = jflat // nx

                    # store (i, j) and (j, i) (symmetric matrix)
                    if jr >= ir:
                        iiR = iR[jr - ir, ir]
                    else:
                        iiR = iR[abs(jr - ir), jr]
                    if jz >= iz:
                        iiZ = iZ[jz - iz, iz]
                    else:
                        iiZ = iZ[abs(jz - iz), jz]

                    data[i0:i0+2] = iiR * iiZ
                    row[i0:i0+2] = (indbs[iflat], indbs[jflat])
                    column[i0:i0+2] = (indbs[jflat], indbs[iflat])
                    i0 += 2

        assert i0 == nbtot
        opmat = scpsp.csc_matrix((data, (row, column)), shape=shape)

    # ------------
    # D1N2

    elif operator == 'D1N2':

        # Treat separately discrete case
        if deg == 0:
            gradR, gradZ = _operators._D1_Deg0_2d(
                knotsx_mult=knotsx_mult,
                knotsy_mult=knotsy_mult,
                cropbs=cropbs,
                cropbs_flat=cropbs_flat,
                nx=nx,
                ny=ny,
                nbs=nbs,
            )

            # surface elements
            dZ = np.repeat(knotsy_mult[1:] - knotsy_mult[:-1], nx)
            if geometry == 'linear':
                dR = np.tile(
                    knotsx_mult[1:] - knotsx_mult[:-1],
                    ny,
                )
            else:
                dR = np.tile(
                    0.5*(knotsx_mult[1:]**2 - knotsx_mult[:-1]**2),
                    ny,
                )

            dS = dR*dZ
            if cropbs_flat is not False:
                dS = dS[cropbs_flat]

            # Does not seem to give positive definite matrix, to be checked ?
            if returnas_element is True:
                opmat = (
                    scpsp.csc_matrix(gradR*np.sqrt(dS[:, None])),
                    scpsp.csc_matrix(gradZ*np.sqrt(dS[:, None])),
                )

            else:
                opmat = (
                    scpsp.csc_matrix((gradR.T.dot(gradR))*(dS[:, None])),
                    scpsp.csc_matrix((gradZ.T.dot(gradZ))*(dS[:, None])),
                )

        else:

            # pre-compute integrals for exact operator deg >= 1
            if deg == 1:
                idR = _operators._D1N2_Deg1(knotsx_mult, geom)
                idZ = _operators._D1N2_Deg1(knotsy_mult)
                iR = _operators._D0N2_Deg1(knotsx_mult, geom)
                iZ = _operators._D0N2_Deg1(knotsy_mult)
            elif deg == 2:
                idR = _operators._D1N2_Deg2(knotsx_mult, geom)
                idZ = _operators._D1N2_Deg2(knotsy_mult)
                iR = _operators._D0N2_Deg2(knotsx_mult, geom)
                iZ = _operators._D0N2_Deg2(knotsy_mult)
            elif deg == 3:
                msg = "Integral D1N2 not implemented for deg=3!"
                raise NotImplementedError(msg)

            # set non-diagonal elements
            i0 = 0
            for ir in range(nx):
                for iz in range(ny):

                    iflat = ir + iz*nx
                    if cropbs_flat is not False and not cropbs_flat[iflat]:
                        continue

                    # general case
                    overlapi = overlap[:, iflat][overlap[:, iflat] > iflat]

                    # diagonal element
                    datadR[i0] = idR[0, ir] * iZ[0, iz]
                    datadZ[i0] = iR[0, ir] * idZ[0, iz]
                    row[i0] = indbs[iflat]
                    column[i0] = indbs[iflat]
                    i0 += 1

                    # non-diagonal elements (symmetric)
                    for jflat in overlapi:

                        if cropbs_flat is not False and not cropbs_flat[jflat]:
                            continue

                        jr = jflat % nx
                        jz = jflat // nx

                        # store (i, j) and (j, i) (symmetric matrix)
                        if jr >= ir:
                            iidR = idR[jr - ir, ir]
                            iiR = iR[jr - ir, ir]
                        else:
                            iidR = idR[abs(jr - ir), jr]
                            iiR = iR[abs(jr - ir), jr]
                        if jz >= iz:
                            iidZ = idZ[jz - iz, iz]
                            iiZ = iZ[jz - iz, iz]
                        else:
                            iidZ = idZ[abs(jz - iz), jz]
                            iiZ = iZ[abs(jz - iz), jz]
                        datadR[i0:i0+2] = iidR * iiZ
                        datadZ[i0:i0+2] = iiR * iidZ
                        row[i0:i0+2] = (indbs[iflat], indbs[jflat])
                        column[i0:i0+2] = (indbs[jflat], indbs[iflat])
                        i0 += 2

            assert i0 == nbtot
            opmat = (
                scpsp.csc_matrix((datadR, (row, column)), shape=shape),
                scpsp.csc_matrix((datadZ, (row, column)), shape=shape),
            )

    # ------------
    # D2N2

    elif operator == 'D2N2':

        # pre-compute integrals
        if deg in [0, 1]:
            msg = f"degree {deg} too low for operator {operator}"
            raise Exception(msg)
        if deg == 2:
            id2R = _operators._D2N2_Deg2(knotsx_mult, geom)
            id2Z = _operators._D2N2_Deg2(knotsy_mult)
            idR = _operators._D1N2_Deg2(knotsx_mult, geom)
            idZ = _operators._D1N2_Deg2(knotsy_mult)
            iR = _operators._D0N2_Deg2(knotsx_mult, geom)
            iZ = _operators._D0N2_Deg2(knotsy_mult)
        elif deg == 3:
            msg = "Integral D2N2 not implemented for deg=3!"
            raise NotImplementedError(msg)

        # set non-diagonal elements
        i0 = 0
        for ir in range(nx):
            for iz in range(ny):

                iflat = ir + iz*nx
                if cropbs_flat is not False and not cropbs_flat[iflat]:
                    continue

                # general case
                overlapi = overlap[:, iflat][overlap[:, iflat] > iflat]

                # diagonal element
                datad2R[i0] = id2R[0, ir] * iZ[0, iz]
                datad2Z[i0] = iR[0, ir] * id2Z[0, iz]
                datadRZ[i0] = idR[0, ir] * idZ[0, iz]
                row[i0] = indbs[iflat]
                column[i0] = indbs[iflat]
                i0 += 1

                # non-diagonal elements (symmetric)
                for jflat in overlapi:

                    if cropbs_flat is not False and not cropbs_flat[jflat]:
                        continue

                    jr = jflat % nx
                    jz = jflat // nx

                    # store (i, j) and (j, i) (symmetric matrix)
                    if jr >= ir:
                        iid2R = id2R[jr - ir, ir]
                        iidR = idR[jr - ir, ir]
                        iiR = iR[jr - ir, ir]
                    else:
                        iid2R = id2R[abs(jr - ir), jr]
                        iidR = idR[abs(jr - ir), jr]
                        iiR = iR[abs(jr - ir), jr]
                    if jz >= iz:
                        iid2Z = id2Z[jz - iz, iz]
                        iidZ = idZ[jz - iz, iz]
                        iiZ = iZ[jz - iz, iz]
                    else:
                        iid2Z = id2Z[abs(jz - iz), jz]
                        iidZ = idZ[abs(jz - iz), jz]
                        iiZ = iZ[abs(jz - iz), jz]
                    datad2R[i0:i0+2] = iid2R * iiZ
                    datad2Z[i0:i0+2] = iiR * iid2Z
                    datadRZ[i0:i0+2] = iidR * iidZ
                    row[i0:i0+2] = (indbs[iflat], indbs[jflat])
                    column[i0:i0+2] = (indbs[jflat], indbs[iflat])
                    i0 += 2

        assert i0 == nbtot
        opmat = (
            scpsp.csc_matrix((datad2R, (row, column)), shape=shape),
            scpsp.csc_matrix((datad2Z, (row, column)), shape=shape),
            scpsp.csc_matrix((datadRZ, (row, column)), shape=shape),
        )

    # ------------
    # D3N2

    elif operator == 'D3N2' and deg == 3:

        raise NotImplementedError("Integral D3N2 not implemented for deg=3!")

    return opmat, operator, geometry