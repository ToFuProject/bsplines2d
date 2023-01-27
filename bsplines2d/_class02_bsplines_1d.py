# -*- coding: utf-8 -*-


# Built-in
import itertools as itt


import numpy as np
import scipy.interpolate as scpinterp


from . import _utils_bsplines


# #############################################################################
# #############################################################################
#                   class
# #############################################################################


class UnivariateSpline():
    """ Subclass for tofu

    Defined from knots (unique) and deg
    coefs set to 1 by default

    Used self.set_coefs() to update
    """

    def __init__(self, knots=None, deg=None):

        assert deg in [0, 1, 2, 3], deg

        # get knots pr bs
        self._get_knots_per_bs_for_basis_elements(
            knots=knots,
            deg=deg,
        )

        # get nbs
        self.shapebs = (self.nbs,)

        # deg
        self.deg = deg

    def _get_knots_per_bs_for_basis_elements(
        self,
        knots=None,
        deg=None,
    ):

        # ------------------------------------
        # get knots per bs in radius direction

        knots_per_bs = _utils_bsplines._get_knots_per_bs(
            knots, deg=deg, returnas='data',
        )
        knots_with_mult, nbs = _utils_bsplines._get_knots_per_bs(
            knots, deg=deg, returnas='data', return_unique=True,
        )

        if nbs != knots_per_bs.shape[1]:
            msg = "Inconsistent nb. of splines"
            raise Exception(msg)

        # ----------------
        # Pre-compute bsplines basis elements

        lbs = [
            scpinterp.BSpline.basis_element(
                knots_per_bs[:, ii],
                extrapolate=False,
            )
            for ii in range(nbs)
        ]

        # ----------------
        # bsplines centers

        cents_per_bs = _utils_bsplines._get_cents_per_bs(
            0.5*(knots[1:] + knots[:-1]),
            deg=deg,
            returnas='data',
        )

        # ----------------
        # bsplines apex

        apex_per_bs = _utils_bsplines._get_apex_per_bs(
            knots=knots,
            knots_per_bs=knots_per_bs,
            deg=deg,
        )

        # ------
        # store

        self.knots = knots
        self.knots_with_mult = knots_with_mult
        self.knots_per_bs = knots_per_bs
        self.cents_per_bs = cents_per_bs
        self.apex_per_bs = apex_per_bs
        self.nbs = nbs
        self.lbs = lbs

    def _check_shape(self, x0=None, coefs=None, axis=None):
        return tuple(np.r_[
            coefs.shape[:axis], x0.shape, coefs.shape[axis+1]
        ])

    def _check_coefs(self, coefs=None, axis=None):
        """ None for ev_details, (nt, shapebs) for sum """
        assert coefs.shape[axis] == self.nbs

    def __call__(
        self,
        # interp points
        x0=None,
        # coefs
        coefs=None,
        axis=None,
        # options
        val_out=None,
        deriv=None,
        # slicing
        sli_c=None,
        sli_v=None,
        sli_o=None,
        shape_v=None,
        shape_o=None,
        # for compatibility (unused)
        **kwdargs,
    ):
        """ Assumes

        coefs.shape = (..., nbs, ...)

        """

        # ------------
        # check inputs

        # axis
        axis = axis[0]

        # coefs
        self._check_coefs(coefs=coefs, axis=axis)

        # x0, shape, val
        val = np.zeros(shape_v, dtype=float)

        # ------------
        # compute

        for ind in itt.product(*[range(aa) for aa in shape_o]):

            # slices
            slic = sli_c(ind)
            sliv = sli_v(ind)

            # call be called on any shape of x0
            val[sliv] = scpinterp.BSpline(
                self.knots_with_mult,
                coefs[slic],
                self.deg,
                axis=0,
                extrapolate=False,
            )(x0, nu=deriv)

        # clean out-of-mesh
        if val_out is not False:
            slio = sli_o((x0 < self.knots[0]) | (x0 > self.knots[-1]))
            val[slio] = val_out

        return val

    # TBF
    def ev_details(
        self,
        # coordinates
        x0=None,
        # options
        indbs_tf=None,
        deriv=None,
        val_out=None,
        # for compatibility (unused)
        **kwdargs,
    ):
        """ Assumes

        indbs_tf = flat array of int indices

        """

        # ------------
        # check inputs

        if indbs_tf is None:
            nbs = self.nbs
        else:
            indbs_tf = np.atleast_1d(indbs_tf).ravel()
            assert 'int' in indbs_tf.dtype.name, indbs_tf
            assert np.unique(indbs_tf).size == indbs_tf.size
            nbs = indbs_tf.size

        # -------
        # prepare

        shape = tuple(np.r_[x0.shape, nbs])
        val = np.full(shape, val_out)

        # ------------
        # compute

        ni = 0
        for ii in range(self.nbs):

            if indbs_tf is not None and ii not in indbs_tf:
                continue

            iok = (
                (x0 >= self.knots_per_bs[0, ii])
                & ((x0 < self.knots_per_bs[-1, ii]))
            )
            if np.any(iok):
                val[iok, ni] = self.lbs[ii](x0[iok], nu=deriv)
            ni += 1

        return val

    def get_constraints_deriv(
        self,
        deriv=None,
        x0=None,
        val=None,
    ):
        """
        To set constraints on a derivative
        Return indices of bsplines + coefs + offset

        Assumes:
            - deriv in ['deriv0', 'deriv1']
            - rad and val are 1d arrays of the same shape

        return as flattened nbsplines indexing

        """
        # ------------
        # check inputs

        ld = ['deriv0', 'deriv1']
        if deriv not in ld:
            msg = f"Arg deriv must be in {ld}!\n Provided: {deriv}"
            raise Exception(msg)

        # --------
        # compute

        # coefs per radius per bs (nrad, nbs)
        ideriv = int(deriv[-1])
        vv = self.ev_details(x0=x0, deriv=ideriv, val_out=0.)

        # check conflicts
        indok = (vv != 0)
        if np.unique(indok, axis=0).shape[0] < indok.shape[0]:
            msg = f"Conflicting constraints on {deriv}:\n{indok}"
            raise Exception(msg)

        coefs = vv
        offset = np.repeat(val[:, None], self.nbs, axis=1)

        return indok, coefs, offset

    # -----------------
    # operator methods
    # -----------------

    def get_overlap(self):
        return _get_overlap(
            deg=self.degrees[0],
            knots=self.knots_per_bs,
            shapebs=self.shapebs,
        )

    def get_operator(
        self,
        operator=None,
        geometry=None,
        cropbs_flat=None,
        # specific to deg = 0
        cropbs=None,
        centered=None,
        # to return gradR, gradZ, for D1N2 deg 0, for tomotok
        returnas_element=None,
    ):
        """ Get desired operator """

        msg = (
            "Operator not implemented yet for 1d bsplines!"
        )
        raise NotImplementedError(msg)


# #############################################################################
# #############################################################################
#                       Mesh2DPolar - bsplines - overlap
# #############################################################################


def _get_overlap(
    deg=None,
    knots=None,
    shapebs=None,
):
    raise NotImplementedError()


# #############################################################################
# #############################################################################
#                   Main
# #############################################################################


def get_bs_class(
    deg=None,
    knots=None,
    coll=None,
):

    # ----------------
    # Define functions

    return UnivariateSpline(
        knots=knots,
        deg=deg,
    )
