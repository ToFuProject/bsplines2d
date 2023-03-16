# -*- coding: utf-8 -*-


# Built-in
import itertools as itt


import numpy as np
import scipy.interpolate as scpinterp


from . import _utils_bsplines
from . import _class02_bsplines_operators_1d as _operators_1d


# ################################################################
# ################################################################
#                   class
# ################################################################


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

    def _check_coefs(self, coefs=None, axis=None):
        """ None for ev_details, (nt, shapebs) for sum """
        if coefs.shape[axis[0]] != self.nbs:
            msg = (
                "Arg coefs has wrong shape!\n"
                f"\t- coefs.shape = {coefs.shape}\n"
                f"\t- coefs.shape[{axis[0]}] != {self.nbs}\n"
            )
            raise Exception(msg)

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
        sli_x=None,
        sli_v=None,
        sli_o=None,
        indokx0=None,
        shape_v=None,
        shape_o=None,
        dref_com=None,
        # for compatibility (unused)
        **kwdargs,
    ):
        """ Assumes

        coefs.shape = (..., nbs, ...)

        """

        # ------------
        # check inputs

        # coefs
        self._check_coefs(coefs=coefs, axis=axis)

        # x0, shape, val
        val = np.zeros(shape_v, dtype=float)

        # ------------
        # compute

        for ind in itt.product(*[range(aa) for aa in shape_o]):

            # slices
            slic = sli_c(
                ind,
                axis=axis,
                ddim=coefs.ndim,
            )

            slix = sli_x(
                ind,
                indokx0=indokx0,
                ix=dref_com['ix'],
                iother=dref_com['iother'],
            )

            sliv = sli_v(
                ind,
                indokx0=indokx0,
                ddim=coefs.ndim,
                axis=axis,
                ix=dref_com['ix'],
                iother=dref_com['iother'],
            )

            # call be called on any shape of x0
            val[sliv] = scpinterp.BSpline(
                self.knots_with_mult,
                coefs[slic],
                self.deg,
                axis=0,
                extrapolate=False,
            )(x0[slix], nu=deriv)

        # clean out-of-mesh
        if dref_com['ix'] is None and val_out is not False:
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

    # -----------------
    # constraints methods
    # -----------------

    def get_constraints_out_rlim(
        self,
        rlim=None,
        rm=None,
        lim=None,
    ):
        """
        Return indices of bslines fully out of rlim (min or max)

        Assumes rlim is a 1d array, one per time step
        """


        # ------------
        # check inputs

        if rm not in ['rmin', 'rmax']:
            msg = f"Invalid rm!\nShould be in ['rmax', 'rmin']\nProvided: {rm}"
            raise Exception(msg)

        if np.isscalar(rlim):
            rlim = np.r_[rlim]

        if not (isinstance(rlim, np.ndarray) and rlim.ndim == 1):
            msg = "rlim must ba a 1d array!"
            raise Exception(msg)

        # ------------
        # get index of bsplines out of limits

        kpbs = self.knots_per_bs
        if rm == 'rmax':
            if lim == 'allout':
                ind = np.all(
                    kpbs[None, :, :] > rlim[:, None, None],
                    axis=1,
                )
            elif lim == 'outer':
                ind = np.sum(
                    kpbs[None, :, :] > rlim[:, None, None],
                    axis=1,
                ) > 1
            elif lim == 'inner':
                ind = np.any(
                    kpbs[None, :, :] > rlim[:, None, None],
                    axis=1,
                )

        else:
            if lim == 'allin':
                ind = np.all(
                    kpbs[None, :, :] < rlim[:, None, None],
                    axis=1,
                )
            elif lim == 'inner':
                ind = np.sum(
                    kpbs[None, :, :] < rlim[:, None, None],
                    axis=1,
                ) > 1
            elif lim == 'outer':
                ind = np.any(
                    kpbs[None, :, :] < rlim[:, None, None],
                    axis=1,
                )

        # if self.knotsa is None:
            # pass
        # else:
            # ind = np.repeat(ind, self.nbs_a_per_r, axis=1)

        # ------------
        # get coefs / offset

        offset = np.zeros((self.nbs,), dtype=float)

        return ind, offset

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
            deg=self.deg,
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

        return _operators_1d.get_operators(
            deg=self.deg,
            operator=operator,
            geometry=geometry,
            knots_mult=self.knots_with_mult,
            knots_per_bs=self.knots_per_bs,
            overlap=self.get_overlap(),
            # specific to deg = 0
            centered=centered,
            # to return gradR, gradZ, for D1N2 deg 0, for tomotok
            returnas_element=returnas_element,
        )


# ################################################################
# ################################################################
#                       overlap
# ################################################################


def _get_overlap(
    deg=None,
    knots=None,
    shapebs=None,
):
    # nb of overlapping, inc. itself in 1d
    nbs, = shapebs
    ind0 = np.arange(0, nbs)

    # complete
    ntot = 2*deg + 1

    add0 = np.arange(-deg, deg+1)

    inter = ind0[None, :] + add0[:, None]

    # purge
    ineg = (inter < 0) | (inter >= nbs)
    inter[ineg] = -1

    return inter


# ################################################################
# ################################################################
#                   Main
# ################################################################


def get_bs_class(
    deg=None,
    knots=None,
):

    # ----------------
    # Define functions

    return UnivariateSpline(
        knots=knots,
        deg=deg,
    )
