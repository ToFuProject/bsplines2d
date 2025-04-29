

import numpy as np
import datastock as ds


# ##############################################################
# ##############################################################
#                   Get knots per bsplines 1d
# ##############################################################


def _get_knots_per_bs(
    knots,
    deg=None,
    returnas=None,
    return_unique=None,
    period=None,
):

    # -----------
    # check input
    # -----------

    returnas = ds._generic_check._check_var(
        returnas, 'returnas',
        types=str,
        default='data',
        allowed=['ind', 'data'],
    )
    return_unique = ds._generic_check._check_var(
        return_unique, 'return_unique',
        types=bool,
        default=False,
    )

    if period is None:
        period = False

    # --------
    # prepare
    # --------

    nkpbs = 2 + deg
    size = knots.size

    if period is False:
        if size < 1 - deg:
            msg = (
                "For the desired degree ({deg}), "
                f"a minimum of {2 - deg} poloidal knots is necessary\n"
                f"Provided: {knots}"
            )
            raise Exception(msg)
        nbs = size - 1 + deg

    else:
        if size < nkpbs - 1:
            msg = (
                f"For the desired degree ({deg}), with a periodic mesh, "
                f"a minimum of {nkpbs - 1} knots is necessary\n"
                f"Provided: {knots}"
            )
            raise Exception(msg)

        nbs = size
        if deg == 0 and size == nkpbs - 1:
            msg = "Using 2 pts for a deg = 0 bsplines leads to bspline!"
            raise Exception(msg)

    # --------------------
    # compute unique knots
    # --------------------

    if return_unique:
        if deg == 0:
            knots_per_bs = np.arange(0, size)
        elif deg == 1:
            knots_per_bs = np.r_[0, np.arange(0, size), size-1]
        elif deg == 2:
            knots_per_bs = np.r_[0, 0, np.arange(0, size), size-1, size-1]
        elif deg == 3:
            knots_per_bs = np.r_[
                0, 0, 0, np.arange(0, size), size-1, size-1, size-1,
            ]

    # --------------------
    # compute multiplicity
    # --------------------

    else:
        knots_per_bs = np.zeros((nkpbs, nbs), dtype=int)

        # -------
        # deg = 0

        if deg == 0:
            if period is False:
                knots_per_bs[:, :] = np.array([
                    np.arange(0, size-1),
                    np.arange(1, size),
                ])
            else:
                if cents[0] > knots[0]:
                    knots_per_bs[:, :] = np.array([
                        np.arange(0, size),
                        np.r_[np.arange(1, size), 0],
                    ])
                else:
                    knots_per_bs[:, :] = np.array([
                        size-1, np.arange(0, size-1),
                        np.arange(0, size),
                    ])

        # -------
        # deg = 1

        elif deg == 1:
            if period is False:
                knots_per_bs[:, 1:-1] = np.array([
                    np.arange(0, size-2),
                    np.arange(1, size-1),
                    np.arange(2, size),
                ])
                knots_per_bs[:, 0] = [0, 0, 1]
                knots_per_bs[:, -1] = [-2, -1, -1]

            else:
                if size == nkpbs - 1:
                    knots_per_bs[:, :] = np.array([
                        np.r_[0, 1],
                        np.r_[1, 0],
                        np.r_[0, 1],
                    ])
                else:
                    knots_per_bs[:, :] = np.array([
                        np.arange(0, size),
                        np.r_[np.arange(1, size), 0],
                        np.r_[np.arange(2, size), 0, 1],
                    ])

        # -------
        # deg = 2

        elif deg == 2:

            if period is False:
                knots_per_bs[:, 2:-2] = np.array([
                    np.arange(0, size-3),
                    np.arange(1, size-2),
                    np.arange(2, size-1),
                    np.arange(3, size),
                ])
                knots_per_bs[:, 0] = [0, 0, 0, 1]
                knots_per_bs[:, 1] = [0, 0, 1, 2]
                knots_per_bs[:, -2] = [-3, -2, -1, -1]
                knots_per_bs[:, -1] = [-2, -1, -1, -1]

            else:
                if size == nkpbs - 1:
                    knots_per_bs[:, :] = np.array([
                        np.arange(0, size),
                        np.r_[np.arange(1, size), 0],
                        np.r_[np.arange(2, size), 0, 1],
                        np.arange(0, size),
                    ])
                else:
                    knots_per_bs[:, :] = np.array([
                        np.arange(0, size),
                        np.r_[np.arange(1, size), 0],
                        np.r_[np.arange(2, size), 0, 1],
                        np.r_[np.arange(3, size), 0, 1, 2],
                    ])

        # -------
        # deg = 3

        elif deg == 3:

            if period is False:
                knots_per_bs[:, 3:-3] = np.array([
                    np.arange(0, size-4),
                    np.arange(1, size-3),
                    np.arange(2, size-2),
                    np.arange(3, size-1),
                    np.arange(4, size),
                ])
                knots_per_bs[:, 0] = [0, 0, 0, 0, 1]
                knots_per_bs[:, 1] = [0, 0, 0, 1, 2]
                knots_per_bs[:, 2] = [0, 0, 1, 2, 3]
                knots_per_bs[:, -3] = [-4, -3, -2, -1, -1]
                knots_per_bs[:, -2] = [-3, -2, -1, -1, -1]
                knots_per_bs[:, -1] = [-2, -1, -1, -1, -1]

            else:
                if size == nkpbs - 1:
                    knots_per_bs[:, :] = np.array([
                        np.arange(0, size),
                        np.r_[np.arange(1, size), 0],
                        np.r_[np.arange(2, size), 0, 1],
                        np.r_[np.arange(3, size), 0, 1, 2],
                        np.arange(0, size),
                    ])
                else:
                    knots_per_bs[:, :] = np.array([
                        np.arange(0, size),
                        np.r_[np.arange(1, size), 0],
                        np.r_[np.arange(2, size), 0, 1],
                        np.r_[np.arange(3, size), 0, 1, 2],
                        np.r_[np.arange(4, size), 0, 1, 2, 3],
                    ])

    # ----------
    # return

    if returnas == 'data':
        knots_per_bs = knots[knots_per_bs]

    if return_unique:
        return knots_per_bs, nbs
    else:
        return knots_per_bs


# ##############################################################
# ##############################################################
#                   Get cents per bsplines 1d
# ##############################################################


def _get_cents_per_bs(
    cents,
    knots=None,
    deg=None,
    returnas=None,
    period=None,
):

    # ------------
    # check inputs
    # ------------

    returnas = ds._generic_check._check_var(
        returnas, 'returnas',
        types=str,
        default='data',
        allowed=['ind', 'data'],
    )

    if period is None:
        period = False

    # -------
    # prepare
    # -------

    nkpbs = 1 + deg
    size = cents.size
    if period is False:
        nbs = size + deg
    else:
        nbs = size
    cents_per_bs = np.zeros((nkpbs, nbs), dtype=int)

    # -------
    # compute
    # -------

    # -------
    # deg = 0

    if deg == 0:
        cents_per_bs[0, :] = np.arange(0, size)

    # -------
    # deg = 1

    elif deg == 1:
        if period is False:
            cents_per_bs[...] = np.array([
                np.r_[0, np.arange(0, size-1), -1],
                np.r_[0, np.arange(1, size), -1],
            ])
        else:
            cents_per_bs[...] = np.array([
                np.arange(0, size),
                np.r_[np.arange(1, size), 0],
            ])
            if cents[0] < knots[0]:
                cents_per_bs += 1
                cents_per_bs = cents_per_bs % size

    # -------
    # deg = 2

    elif deg == 2:
        if period is False:
            cents_per_bs[:, 2:-2] = np.array([
                np.arange(0, size-2),
                np.arange(1, size-1),
                np.arange(2, size),
            ])
            cents_per_bs[:, 0] = [0, 0, 0]
            cents_per_bs[:, 1] = [0, 0, 1]
            cents_per_bs[:, -2] = [-2, -1, -1]
            cents_per_bs[:, -1] = [-1, -1, -1]

        else:
            cents_per_bs[...] = np.array([
                np.arange(0, size),
                np.r_[np.arange(1, size-1), 0],
                np.r_[np.arange(2, size), 0, 1],
            ])
            if cents[0] < knots[0]:
                cents_per_bs += 1
                cents_per_bs = cents_per_bs % size

    # -------
    # deg = 3

    elif deg == 3:
        if period is False:
            cents_per_bs[:, 3:-3] = np.array([
                np.arange(0, size-3),
                np.arange(1, size-2),
                np.arange(2, size-1),
                np.arange(3, size),
            ])
            cents_per_bs[:, 0] = [0, 0, 0, 0]
            cents_per_bs[:, 1] = [0, 0, 0, 1]
            cents_per_bs[:, 2] = [0, 0, 1, 2]
            cents_per_bs[:, -3] = [-3, -2, -1, -1]
            cents_per_bs[:, -2] = [-2, -1, -1, -1]
            cents_per_bs[:, -1] = [-1, -1, -1, -1]
        else:
            cents_per_bs[...] = np.array([
                np.arange(0, size),
                np.r_[np.arange(1, size-1), 0],
                np.r_[np.arange(2, size), 0, 1],
                np.r_[np.arange(3, size), 0, 1, 2],
            ])
            if cents[0] < knots[0]:
                cents_per_bs += 1
                cents_per_bs = cents_per_bs % size

    # ------
    # return

    if returnas == 'data':
        cents_per_bs = cents[cents_per_bs]

    return cents_per_bs


# ###############################################################
# ###############################################################
#                   Get apex positions per bsplines 1d
# ###############################################################


def _get_apex_per_bs(
    knots=None,
    knots_per_bs=None,
    cents_per_bs=None,
    deg=None,
    period=None,
):

    # -------
    # prepare

    if period is None:
        period = False

    nkpbs, nbs = knots_per_bs.shape

    # -------------
    # compute basis

    if nkpbs % 2 == 0:
        ii = int(nkpbs/2) - 1
        apex = cents_per_bs[ii, :]

    else:
        ii = int((nkpbs-1) / 2)
        apex = knots_per_bs[ii, :]

    # ------
    # adjust

    if period is False:
        # manage edges
        if deg == 1:
            apex[:deg] = knots[0]
            apex[-deg:] = knots[-1]
        elif deg == 2:
            apex[:deg] = [knots[0], 0.5*(knots[0] + knots[1])]
            apex[-deg:] = [0.5*(knots[-2] + knots[-1]), knots[-1]]
        elif deg == 3:
            apex[:deg] = [knots[0], 0.5*(knots[0]+knots[1]), knots[1]]
            apex[-deg:] = [knots[-2], 0.5*(knots[-2]+knots[-1]), knots[-1]]

    return apex


# ###################################################
# ###################################################
#       Get apex positions per bsplines 1d
# ####################################################


def _get_shapes_ind(
    axis=None,
    shape_c=None,
    shape_x=None,
):

    # shape of values
    shape_v = tuple(np.r_[
        shape_c[:axis[0]], shape_x, shape_c[axis[-1]+1:]
    ].astype(int))

    # axis_v
    axis_v = axis[0] + np.arange(len(shape_x))

    # ind for coefs
    ind_c = np.arange(0, len(shape_c))
    ind_c[axis[-1]:] -= len(axis)

    # ind for values
    ind_v = np.arange(0, len(shape_v))
    ind_v[axis_v[-1]:] -= len(axis_v)

    # shape of other dimensions (common to coefs and values)
    shape_o = tuple([
        ss for ii, ss in enumerate(shape_c)
        if ii not in axis
    ])

    return shape_v, axis_v, ind_c, ind_v, shape_o


def _get_slice_cx(
    axis=None,
    shape=None,
    ind_cv=None,
    reverse=None,
):
    """ Return coefs slicing

    For mtype = 'tri', revser=True
        => return a list

    for other cases, return a function that returns a tuple

    """

    if reverse is True:
            sli = [
                ind_cv if ii == axis[0] else slice(None)
                for ii in range(len(shape))
                if ii not in axis[1:]
            ]

    else:
        def sli(ind, axis=axis, shape=shape, ind_cv=ind_cv):
            return tuple([
                slice(None) if ii in axis else ind[ind_cv[ii]]
                for ii in range(len(shape))
            ])

    return sli


def _get_slice_out(axis=None, shape_c=None):

    def func(indout, axis=axis, shape_c=shape_c):
        return tuple([
            indout if ii == axis[0] else slice(None)
            for ii in range(len(shape_c))
            if ii not in axis[1:]
        ])

    return func
