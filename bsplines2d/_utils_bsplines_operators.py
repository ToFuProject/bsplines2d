# -*- coding: utf-8 -*-


import numpy as np
import datastock as ds


# ###################################################
# ###################################################
#               D0N1
# ###################################################


def _D0N1_Deg0(kn, geom='linear'):
    if geom  == 'linear':
        return kn[1, :] - kn[0, :]
    else:
        return 0.5 * (kn[1, :]**2 - kn[0, :]**2)


def _D0N1_Deg1(kn, geom='linear'):
    if geom  == 'linear':
        return 0.5 * (kn[2, :] - kn[0, :])
    else:
        return (
            kn[2, :]**2 - kn[0, :]**2
            + kn[1, :]*(kn[2, :]-kn[0, :])
        ) / 6.


def _D0N1_Deg2(kn, geom='linear'):

    # --------
    # prepare

    # initialize
    in1 = np.zeros(kn.shape[1], dtype=float)
    in21 = np.zeros(kn.shape[1], dtype=float)
    in22 = np.zeros(kn.shape[1], dtype=float)
    in3 = np.zeros(kn.shape[1], dtype=float)

    # prepare indices
    i0 = kn[2, :] != kn[0, :]
    i1 = kn[3, :] != kn[1, :]

    # ------------------
    # compute components

    if geom == 'linear':
        in1[i0] = (
            (kn[1, i0] - kn[0, i0])**2
            / (3.*(kn[2, i0] - kn[0, i0]))
        )

        in21[i0] = (
            (
                kn[2, i0]**2
                - 2. * kn[1, i0]**2
                + kn[1, i0] * kn[2, i0]
                + 3. * kn[0, i0] * (kn[1, i0] - kn[2, i0])
            )
            / (6. * (kn[2, i0] - kn[0, i0]))
        )

        in22[i1] = (
            (
                -2. * kn[2, i1]**2
                + kn[1, i1]**2
                + kn[1, i1] * kn[2, i1]
                + 3. * kn[3, i1] * (kn[2, i1] - kn[1, i1])
            )
            / (6.*(kn[3, i1] - kn[1, i1]))
        )

        in3[i1] = (
            (kn[3, i1] - kn[2, i1])**2
            / (3.*(kn[3, i1] - kn[1, i1]))
        )

    else:
        in1[i0] = (
            (
                3.*kn[1, i0]**3
                + kn[0, i0]**3
                - 5.*kn[0, i0] * kn[1, i0]**2
                + kn[0, i0]**2 * kn[1, i0]
            )
            / (12. * (kn[2, i0] - kn[0, i0]))
        )
        in21[i0] = (
            (
                kn[2, i0]**3
                - 3.*kn[1, i0]**3
                + kn[1, i0]**2 * kn[2, i0]
                + kn[1, i0] * kn[2, i0]**2
                - 2.*kn[0, i0] * kn[2, i0]**2
                - 2.*kn[0, i0] * kn[1, i0] * kn[2, i0]
                + 4.*kn[0, i0] * kn[1, i0]**2
            )
            / (12. * (kn[2, i0] - kn[0, i0]))
        )
        in22[i1] = (
            (
                -3.*kn[2, i1]**3
                + kn[1, i1]**3
                + kn[1, i1] * kn[2, i1]**2
                + kn[1, i1]**2 * kn[2, i1]
                + 4.*kn[2, i1]**2 * kn[3, i1]
                - 2.*kn[1, i1]*kn[2, i1]*kn[3, i1]
                - 2.*kn[1, i1]**2 * kn[3, i1]
            )
            / (12. * (kn[3, i1] - kn[1, i1]))
        )
        in3[i1] = (
            (
                kn[3, i1]**3
                + 3.*kn[2, i1]**3
                - 5.*kn[2, i1]**2 * kn[3, i1]
                + kn[2, i1]*kn[3, i1]**2
            ) / (12. * (kn[3, i1] - kn[1, i1]))
        )

    return in1, in21, in22, in3


# #####################################################
# #####################################################
#               D0N2 - Deg 0
# #####################################################


def _D0N2_Deg0(knots, geom='linear'):
    if geom == 'linear':
        return knots[1, :] - knots[0, :]
    else:
        return 0.5 * (knots[1, :]**2 - knots[0, :]**2)


# #####################################################
# #####################################################
#               D0N2 - Deg 1
# #####################################################


def _D0N2_Deg1(knots, geom='linear'):

    return np.array([
        _D0N2_Deg1_full(
            knots[:-2],
            knots[1:-1],
            knots[2:],
            geom,
        ),
        _D0N2_Deg1_2(
            knots[1:-1],
            knots[2:],
            geom,
        ),
    ])


def _D0N2_Deg1_full(k0, k1, k2, geom='linear'):
    """ from 1d knots, return int_0^2 x b**2(x) dx """
    if geom == 'linear':
        return (k2 - k0) / 3.
    else:
        intt = np.zeros((k0.size,))
        intt[1:] += (
            (
                3. * k1**3
                - 5.*k0*k1**2
                + k1*k0**2
                + k0**3
            )[1:]
            / (12. * (k1 - k0))[1:]
        )
        intt[:-1] = (
            + (
                3.*k1**3
                - 5.*k2*k1**2
                + k1*k2**2
                + k2**3
            )[:-1]
            / (12. * (k2 - k1))[:-1]
        )
        return intt


def _D0N2_Deg1_2(k1, k2, geom='linear'):
    """ from 1d knots, return int_0^2 x b**2(x) dx """
    if geom == 'linear':
        return (k2 - k1) / 6.
    else:
        return (k2**2 - k1**2) / 12.


# #####################################################
# #####################################################
#               D0N2 - Deg 2
# #####################################################


def _D0N2_Deg2(knots, geom='linear'):

    if geom == 'linear':
        ffull = _D0N2_Deg2_full_linear
        f3 = _D0N2_Deg2_3_linear
        f2 = _D0N2_Deg2_2_linear
    else:
        ffull = _D0N2_Deg2_full_toroidal
        f3 = _D0N2_Deg2_3_toroidal
        f2 = _D0N2_Deg2_2_toroidal

    integ = np.array([
        ffull(
            knots[:-3],
            knots[1:-2],
            knots[2:-1],
            knots[3:],
        ),
        f3(
            knots[:-3],
            knots[1:-2],
            knots[2:-1],
            knots[3:],
            np.r_[knots[4:], np.nan],
        ),
        f2(
            knots[1:-2],
            knots[2:-1],
            knots[3:],
            np.r_[knots[4:], np.nan],
        ),
    ])
    return integ


def _D0N2_Deg2_full_linear(k0, k1, k2, k3):
    """ from 1d knots, return int_0^3 b**2(x) dx """
    intt = np.zeros((k0.size,))
    intt[1:] += (
        (k1 - k0)[1:]**3 / (5.*(k2 - k0)[1:]**2)
        + (k2 - k1)[1:]
        * (
            10.*k0**2 + 6.*k1**2 + 3.*k1*k2 + k2**2 - 5.*k0*(3.*k1 + k2)
        )[1:] / (30.*(k2 - k0)[1:]**2)
    )
    intt[1:-1] += (
        (k2 - k1)[1:-1]
        * (
            -3.*k1**2 - 4.*k1*k2 - 3.*k2**2 + 5.*k0*(k1 + k2 - 2.*k3)
            + 5.*k3*(k1 + k2)
        )[1:-1] / (60.*(k2 - k0)*(k3 - k1))[1:-1]
    )
    intt[:-1] += (
        (k2 - k1)[:-1]
        * (
            10.*k3**2 + 6.*k2**2 + 3.*k1*k2 + k1**2 - 5.*k3*(3.*k2 + k1)
        )[:-1] / (30.*(k3 - k1)[:-1]**2)
        + (k3 - k2)[:-1]**3 / (5.*(k3 - k1)[:-1]**2)
    )
    return intt


def _D0N2_Deg2_full_toroidal(k0, k1, k2, k3):
    """ from 1d knots, return int_0^3 b**2(x) dx """
    intt = np.zeros((k0.size,))
    intt[1:] += (
        (5.*k1 + k0)[1:]*(k1 - k0)[1:]**3 / (30.*(k2 - k0)[1:]**2)
        + (k2 - k1)[1:]
        * (
            10*k1**3 + 6.*k1**2*k2 + 3.*k1*k2**2
            + k2**3 + 5.*k0**2*(3.*k1 + k2)
            - 4.*k0*(6.*k1**2 + 3.*k1*k2 + k2**2)
        )[1:] / (60.*(k2 - k0)**2)[1:]
    )
    intt[:-1] += (
        (5.*k2 + k3)[:-1]*(k3 - k2)[:-1]**3 / (30.*(k3 - k1)[:-1]**2)
        + (k2 - k1)[1:]
        * (
            10*k2**3 + 6.*k2**2*k1 + 3.*k2*k1**2
            + k1**3 + 5.*k3**2*(3.*k2 + k1)
            - 4.*k3*(6.*k2**2 + 3.*k2*k1 + k1**2)
        )[:-1] / (60.*(k3 - k1)**2)[:-1]
    )
    intt[1:-1] += (
        (k2 - k1)[1:-1]
        * (
            - 2.*k1**3 - 2.*k2**3
            - 3.*k1*k2*(k1 + k2)
            - 5.*k0*k3*(k1 + k2)
            + (k0 + k3)*(3.*k2**2 + 4.*k1*k2 + 3.*k1**2)
        )[1:-1] / (30.*(k3 - k1)*(k2 - k0))[1:-1]
    )

    return intt


def _D0N2_Deg2_3_linear(k0, k1, k2, k3, k4):
    """ from 1d knots, return int_0^3 b**2(x) dx """
    intt = np.zeros((k0.size,))
    intt[1:-1] += (
        (3.*k2 + 2.*k1 - 5.*k0)[1:-1]*(k2 - k1)[1:-1]**2
        / (60.*(k3 - k1)*(k2 - k0))[1:-1]
    )
    intt[:-2] += (
        + (5.*k4 - 2.*k3 - 3.*k2)[:-2]*(k3 - k2)[:-2]**2
        / (60.*(k4 - k2)*(k3 - k1))[:-2]
    )
    intt[:-1] += (
        + (5.*k3 - 4.*k2 - k1)[:-1]*(k2 - k1)[:-1]**2
        / (20.*(k3 - k1)**2)[:-1]
        + (4.*k2 + k3 - 5.*k1)[:-1]*(k3 - k2)[:-1]**2
        / (20.*(k3 - k1)**2)[:-1]
    )
    return intt


# TBC, seems too low compared to full toroidal
def _D0N2_Deg2_3_toroidal(k0, k1, k2, k3, k4):
    """ from 1d knots, return int_0^3 b**2(x) dx """
    intt = np.zeros((k0.size,))
    intt[:-1] = (
        (k2 - k1)[:-1]**2
        * (-10*k2**2 - 4*k1*k2 - k1**2 + 3*k3*(4*k2 + k1))[:-1]
        / (60.*(k3 - k1)**2)[:-1]
        + (k3 - k2)[:-1]**2
        * (k3**2 + 4*k3*k2 + 10*k2**2 - 3*k1*(k3 + 4*k2))[:-1]
        / (60*(k3 - k1)**2)[:-1]
    )
    intt[1:-1] = (
        (k2 - k1)[1:-1]**2
        * (2*k2**2 + 2*k1*k2 + k1**2 - k0*(3.*k2 + 2.*k1))[1:-1]
        / (60.*(k3 - k1)*(k2 - k0))[1:-1]
    )
    intt[:-2] = (
        + (k3 - k2)[:-2]**2
        * (-k3**2 - 2*k3*k2 - 2*k2**2 + k4*(2*k3 + 3*k2))[:-2]
        / (60*(k4 - k2)*(k3 - k1))[:-2]
    )
    return intt


def _D0N2_Deg2_2_linear(k1, k2, k3, k4):
    """ from 1d knots, return int_0^3 b**2(x) dx """
    intt = np.zeros((k1.size,))
    intt[:-2] = (
        (k3 - k2)[:-2]**3
        / (30.*(k4 - k2)*(k3 - k1))[:-2]
    )
    return intt


def _D0N2_Deg2_2_toroidal(k1, k2, k3, k4):
    """ from 1d knots, return int_0^3 b**2(x) dx """
    intt = np.zeros((k1.size,))
    intt[:-2] = (
        (k3 + k2)[:-2]*(k3 - k2)[:-2]**3
        / (60.*(k4 - k2)[:-2]*(k3 - k1)[:-2])
    )
    return intt


# ###############################################################
# ###############################################################
#               D1 - deg = 0 - discrete
# ###############################################################


def _D1_Deg0(
    knots_mult=None,
    nbs=None,
    centered=None,
):
    """ Discrete apprmixation of the gradient for pixels

    Centered when possible
    Non-centered otherwise

    """
    # ------------------
    # check / initialize

    centered = ds._generic_check._check_var(
        centered, 'centered',
        types=bool,
        default=False,
    )

    # initialize output
    datadx = np.zeros((nbs, nbs), dtype=float)

    # positions of centers
    cents = 0.5*(knots_mult[1:] + knots_mult[:-1])

    # --------
    # compute

    if centered is True:
        # Determine points that have 2 neighbours in R (centered)
        n2x = np.ones((nbs,), dtype=bool)
        n2x[0] = False
        n2x[-1] = False

        # points with a neighbours at higher R
        npx = np.ones((nbs,), dtype=bool) & (~n2x)
        npx[-1] = False

        # points with a neighbours at lower R
        nmx = np.ones((nbs,), dtype=bool) & (~n2x) & (~npx)
        nmx[0] = False

        # iterate on each type of point in R
        for ii in n2x.nonzero()[0]:
            dxi = 1./(cents[ii + 1] - cents[ii - 1])
            datadx[ii, ii - 1] = -dxi
            datadx[ii, ii + 1] = dxi
    else:
        # points with a neighbours at higher R
        npx = np.ones((nbs,), dtype=bool)
        npx[-1] = False
        # points with a neighbours at lower R
        nmx = np.ones((nbs,), dtype=bool) & (~npx)
        nmx[0] = False

    # -----------------------------
    # iterate on each type of point

    for ii in npx.nonzero()[0]:
        dxi = 1./(cents[ii + 1] - cents[ii])
        datadx[ii, ii] = -dxi
        datadx[ii, ii + 1] = dxi

    for ii in nmx.nonzero()[0]:
        dxi = 1./(cents[ii] - cents[ii - 1])
        datadx[ii, ii - 1] = -dxi
        datadx[ii, ii] = dxi

    return datadx


def _D1_Deg0_2d(
    knotsx_mult=None,
    knotsy_mult=None,
    cropbs=None,
    cropbs_flat=None,
    nx=None,
    ny=None,
    nbs=None,
    centered=None,
):
    """ Discrete apprmixation of the gradient for pixels

    Centered when possible
    Non-centered otherwise

    """

    # check input
    centered = ds._generic_check._check_var(
        centered, 'centered',
        types=bool,
        default=False,
    )

    # initialize output
    datadR = np.zeros((nbs, nbs), dtype=float)
    datadZ = np.zeros((nbs, nbs), dtype=float)
    if cropbs is False:
        cropbs = np.ones((nx, ny), dtype=bool)

    # positions of centers
    centsR = 0.5*(knotsx_mult[1:] + knotsx_mult[:-1])
    centsZ = 0.5*(knotsy_mult[1:] + knotsy_mult[:-1])

    if centered is True:
        # Determine points that have 2 neighbours in R (centered)
        n2R = np.zeros(cropbs.shape, dtype=bool)
        n2R[1:-1, :] = cropbs[1:-1, :] & cropbs[2:, :] & cropbs[:-2, :]
        # points with a neighbours at higher R
        npR = cropbs & (~n2R)
        npR[-1, :] = False
        npR[:-1, :] &= cropbs[1:, :]
        # points with a neighbours at lower R
        nmR = cropbs & (~n2R) & (~npR)
        nmR[0, :] = False
        nmR[1:, :] &= cropbs[:-1, :]

        # Determine points that have 2 neighbours in Z (centered)
        n2Z = np.zeros(cropbs.shape, dtype=bool)
        n2Z[:, 1:-1] = cropbs[:, 1:-1] & cropbs[:, 2:] & cropbs[:, :-2]
        n2Z[:, 1:-1] = n2Z[:, 1:-1] & n2Z[:, 2:] & n2Z[:, :-2]
        # points with a neighbours at higher Z
        npZ = cropbs & (~n2Z)
        npZ[:, -1] = False
        npZ[:, :-1] &= cropbs[:, 1:]
        # points with a neighbours at lower Z
        nmZ = cropbs & (~n2Z) & (~npZ)
        nmZ[:, 0] = False
        nmZ[:, 1:] &= cropbs[:, :-1]

        # iterate on each type of point in R
        for ir, iz in zip(*n2R.nonzero()):
            iflat = ir + iz*nx
            dRi = 1./(centsR[ir + 1] - centsR[ir - 1])
            datadR[iflat, iflat - 1] = -dRi
            datadR[iflat, iflat + 1] = dRi
        for ir, iz in zip(*n2Z.nonzero()):
            iflat = ir + iz*nx
            dZi = 1./(centsZ[iz + 1] - centsZ[iz - 1])
            datadZ[iflat, iflat - nx] = -dZi
            datadZ[iflat, iflat + nx] = dZi
    else:
        # points with a neighbours at higher R
        npR = np.copy(cropbs)
        npR[-1, :] = False
        npR[:-1, :] &= cropbs[1:, :]
        # points with a neighbours at lower R
        nmR = cropbs & (~npR)
        nmR[0, :] = False
        nmR[1:, :] &= cropbs[:-1, :]

        # points with a neighbours at higher Z
        npZ = np.copy(cropbs)
        npZ[:, -1] = False
        npZ[:, :-1] &= cropbs[:, 1:]
        # points with a neighbours at lower Z
        nmZ = cropbs & (~npZ)
        nmZ[:, 0] = False
        nmZ[:, 1:] &= cropbs[:, :-1]

    # iterate on each type of point in R
    for ir, iz in zip(*npR.nonzero()):
        iflat = ir + iz*nx
        dRi = 1./(centsR[ir + 1] - centsR[ir])
        datadR[iflat, iflat] = -dRi
        datadR[iflat, iflat + 1] = dRi
    for ir, iz in zip(*nmR.nonzero()):
        iflat = ir + iz*nx
        dRi = 1./(centsR[ir] - centsR[ir-1])
        datadR[iflat, iflat - 1] = -dRi
        datadR[iflat, iflat] = dRi

    # iterate on each type of point in Z
    for ir, iz in zip(*npZ.nonzero()):
        iflat = ir + iz*nx
        dZi = 1./(centsZ[iz + 1] - centsZ[iz])
        datadZ[iflat, iflat] = -dZi
        datadZ[iflat, iflat + nx] = dZi
    for ir, iz in zip(*nmZ.nonzero()):
        iflat = ir + iz*nx
        dZi = 1./(centsZ[iz] - centsZ[iz - 1])
        datadZ[iflat, iflat - nx] = -dZi
        datadZ[iflat, iflat] = dZi

    # crop and return
    if cropbs_flat is False:
        return datadR, datadZ
    else:
        return (
            datadR[cropbs_flat, :][:, cropbs_flat],
            datadZ[cropbs_flat, :][:, cropbs_flat],
        )


def _D1_Deg1(
    knots_mult=None,
    nbs=None,
):

    grad = np.zeros((nbs-1, nbs), dtype=float)

    dx = np.diff(knots_mult[1:-1])
    for ii in range(nbs - 1):
        grad[ii, ii+1] = 1./dx[ii]
        grad[ii, ii] = -1./dx[ii]

    return grad


def _D1_Deg2(
    knots_mult=None,
    nbs=None,
):
    grad = np.zeros((nbs-1, nbs), dtype=float)

    raise NotImplementedError()


# ###############################################################
# ###############################################################
#               D1N2 - exact
# ###############################################################


def _D1N2_Deg1(knots, geometry='linear'):

    if geometry == 'linear':
        ffull = _D1N2_Deg1_full_linear
        f2 = _D1N2_Deg1_2_linear
    else:
        ffull = _D1N2_Deg1_full_toroidal
        f2 = _D1N2_Deg1_2_toroidal

    integ = np.array([
        ffull(
            knots[:-2],
            knots[1:-1],
            knots[2:],
        ),
        f2(
            knots[1:-1],
            knots[2:]
        ),
    ])
    return integ


def _D1N2_Deg1_full_linear(k0, k1, k2):
    intt = np.zeros((k0.size,))
    intt[1:] += 1. / (k1 - k0)[1:]
    intt[:-1] += 1. / (k2 - k1)[:-1]
    return intt


def _D1N2_Deg1_full_toroidal(k0, k1, k2):
    intt = np.zeros((k0.size,))
    intt[1:] += (k1 + k0)[1:] / (2.*(k1 - k0))[1:]
    intt[:-1] += (k2 + k1)[:-1] / (2.*(k2 - k1))[:-1]
    return intt


def _D1N2_Deg1_2_linear(k1, k2):
    intt = np.zeros((k1.size,))
    intt[:-1] = -1. / (k2 - k1)[:-1]
    return intt


def _D1N2_Deg1_2_toroidal(k1, k2):
    intt = np.zeros((k1.size,))
    intt[:-1] = - (k2 + k1)[:-1] / (2.*(k2 - k1))[:-1]
    return intt


def _D1N2_Deg2(knots, geometry='linear'):

    if geometry == 'linear':
        ffull = _D1N2_Deg2_full_linear
        f3 = _D1N2_Deg2_3_linear
        f2 = _D1N2_Deg2_2_linear
    else:
        ffull = _D1N2_Deg2_full_toroidal
        f3 = _D1N2_Deg2_3_toroidal
        f2 = _D1N2_Deg2_2_toroidal

    integ = np.array([
        ffull(
            knots[:-3],
            knots[1:-2],
            knots[2:-1],
            knots[3:],
        ),
        f3(
            knots[:-3],
            knots[1:-2],
            knots[2:-1],
            knots[3:],
            np.r_[knots[4:], np.nan],
        ),
        f2(
            knots[1:-2],
            knots[2:-1],
            knots[3:],
            np.r_[knots[4:], np.nan],
        ),
    ])
    return integ


def _D1N2_Deg2_full_linear(k0, k1, k2, k3):
    intt = np.zeros((k0.size,))
    intt[1:] += 4.*(k1 - k0)[1:] / (3.*(k2 - k0)[1:]**2)
    intt[:-1] += 4.*(k3 - k2)[:-1] / (3.*(k3 - k1)[:-1]**2)
    intt[1:-1] += (
        4.*(k2 - k1)[1:-1]
        * (
            k2**2 + k2*k1 + k1**2 + k3**2 + k0*k3 + k0**2
            - k3*(k2 + 2.*k1) - k0*(2.*k2 + k1)
        )[1:-1]
        / (3.*(k3 - k1)[1:-1]**2*(k2 - k0)[1:-1]**2)
    )
    return intt


def _D1N2_Deg2_full_toroidal(k0, k1, k2, k3):
    intt = np.zeros((k0.size,))
    intt[1:] += (3.*k1 + k0)[1:]*(k1 - k0)[1:] / (3.*(k2 - k0)[1:]**2)
    intt[:-1] += (3.*k2 + k3)[:-1]*(k3 - k2)[:-1] / (3.*(k3 - k1)[:-1]**2)
    intt[1:-1] += (
        (k2 - k1)[1:-1]
        * (
            3.*(k2 + k1)*(k2**2 + k1**2)
            + k3**2*(k2 + 3.*k1)
            + k0**2*(3.*k2 + k1)
            - 2.*k3*(k2**2 + 2.*k2*k1 + 3.*k1**2)
            - 2.*k0*(3.*k2**2 + 2.*k2*k1 + k1**2)
            + 2.*k3*k0*(k2 + k1)
        )[1:-1]
        / (3.*(k3 - k1)[1:-1]**2*(k2 - k0)[1:-1]**2)
    )
    return intt


def _D1N2_Deg2_3_linear(k0, k1, k2, k3, k4):
    intt = np.zeros((k0.size,))
    intt[1:-1] += (
        2.*(k2 - k1)[1:-1]
        * (k3 - 2.*k2 - k1 + 2.*k0)[1:-1]
        / (3.*(k3 - k1)**2*(k2 - k0))[1:-1]
    )
    intt[:-2] += (
        2.*(k3 - k2)[:-2]
        * (-2.*k4 + k3 + 2.*k2 - k1)[:-2]
        / (3.*(k4 - k2)*(k3 - k1)**2)[:-2]
    )
    return intt


def _D1N2_Deg2_3_toroidal(k0, k1, k2, k3, k4):
    intt = np.zeros((k0.size,))
    intt[1:-1] += (
        (k2 - k1)[1:-1]
        * (
            - (3.*k2**2 + 2.*k2*k1 + k1**2)
            + k3*(k2 + k1)
            + k0*(3.*k2 + k1)
        )[1:-1]
        / (3.*(k3 - k1)**2*(k2 - k0))[1:-1]
    )
    intt[:-2] += (
        (k3 - k2)[:-2]
        * (
            k3**2 + 2.*k2*k3 + 3.*k2**2
            - k4*(k3 + 3.*k2)
            - k1*(k3 + k2)
        )[:-2]
        / (3.*(k4 - k2)*(k3 - k1)**2)[:-2]
    )
    return intt


def _D1N2_Deg2_2_linear(k1, k2, k3, k4):
    intt = np.zeros((k1.size,))
    intt[:-2] += -(
        2.*(k3 - k2)[:-2]
        / (3.*(k4 - k2)*(k3 - k1))[:-2]
    )
    return intt


def _D1N2_Deg2_2_toroidal(k1, k2, k3, k4):
    intt = np.zeros((k1.size,))
    intt[:-2] += -(
        (k3 + k2)[:-2]*(k3 - k2)[:-2]
        / (3.*(k4 - k2)*(k3 - k1))[:-2]
    )
    return intt


# ###############################################################
# ###############################################################
#               D2N2
# ###############################################################


def _D2N2_Deg2(knots, geometry='linear'):

    if geometry == 'linear':
        ffull = _D2N2_Deg2_full_linear
        f3 = _D2N2_Deg2_3_linear
        f2 = _D2N2_Deg2_2_linear
    else:
        ffull = _D2N2_Deg2_full_toroidal
        f3 = _D2N2_Deg2_3_toroidal
        f2 = _D2N2_Deg2_2_toroidal

    integ = np.array([
        ffull(
            knots[:-3],
            knots[1:-2],
            knots[2:-1],
            knots[3:],
        ),
        f3(
            knots[:-3],
            knots[1:-2],
            knots[2:-1],
            knots[3:],
            np.r_[knots[4:], np.nan],
        ),
        f2(
            knots[1:-2],
            knots[2:-1],
            knots[3:],
            np.r_[knots[4:], np.nan],
        ),
    ])
    return integ


def _D2N2_Deg2_full_linear(k0, k1, k2, k3):
    intt = np.zeros((k0.size,))
    intt[2:] += 4. / ((k2 - k0)**2*(k1 - k0))[2:]
    intt[:-2] += 4. / ((k3 - k2)*(k3 - k1)**2)[:-2]
    intt[1:-1] += (
        4.*(k3 + k2 - k1 - k0)[1:-1]**2
        / ((k3 - k1)**2*(k2 - k1)*(k2 - k0)**2)[1:-1]
    )
    return intt


def _D2N2_Deg2_full_toroidal(k0, k1, k2, k3):
    intt = np.zeros((k0.size,))
    intt[2:] += 2.*(k0 + k0)[2:] / ((k2 - k0)**2*(k1 - k0))[2:]
    intt[:-2] += 2.*(k3 + k2)[2:] / ((k3 - k2)*(k3 - k1)**2)[:-2]
    intt[1:-1] += (
        (2.*(k3 + k2 - k1 - k0)**2*(k2 + k1))[1:-1]
        / ((k3 - k1)**2*(k2 - k1)*(k2 - k0)**2)[1:-1]
    )
    return intt


def _D2N2_Deg2_3_linear(k0, k1, k2, k3, k4):
    intt = np.zeros((k0.size,))
    intt[1:-1] += (
        - 4.*(k3 + k2 - k1 - k0)[1:-1]
        / ((k3 - k1)**2*(k2 - k1)*(k2 - k0))[1:-1]
    )
    intt[:-2] += (
        - 4.*(k4 + k3 - k2 - k1)[:-2]
        / ((k4 - k2)*(k3 - k2)*(k3 - k1)**2)[:-2]
    )
    return intt


def _D2N2_Deg2_3_toroidal(k0, k1, k2, k3, k4):
    intt = np.zeros((k0.size,))
    intt[1:-1] += (
        - 2.*((k3 + k2 - k1 - k0)*(k2 + k1))[1:-1]
        / ((k3 - k1)**2*(k2 - k1)*(k2 - k0))[1:-1]
    )
    intt[:-2] += (
        - 2.*((k4 + k3 - k2 - k1)*(k3 + k2))[:-2]
        / ((k4 - k2)*(k3 - k2)*(k3 - k1)**2)[:-2]
    )
    return intt


def _D2N2_Deg2_2_linear(k1, k2, k3, k4):
    intt = np.zeros((k1.size,))
    intt[:-2] += 4. / ((k4 - k2)*(k3 - k2)*(k3 - k1))[:-2]
    return intt


def _D2N2_Deg2_2_toroidal(k1, k2, k3, k4):
    intt = np.zeros((k1.size,))
    intt[:-2] += 2.*(k3 + k2)[:-2] / ((k4 - k2)*(k3 - k2)*(k3 - k1))[:-2]
    return intt