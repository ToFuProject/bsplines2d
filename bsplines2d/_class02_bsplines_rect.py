# -*- coding: utf-8 -*-


# Built-in
import itertools as itt
import warnings


# Common
import numpy as np
import scipy.interpolate as scpinterp


# specific
from . import _utils_bsplines
from . import _class02_bsplines_operators_rect


if hasattr(scpinterp._bspl, 'evaluate_spline'):
    evaluate_spline = scpinterp._bspl.evaluate_spline

else:
    msg = (
        "\n\n"
        "bsplines2d using a new version of scipy"
        " with no scpinterp._bspl.evaluate_spline()\n"
        "Instead using scpinterp._bspl.evaluate_ndspline()\n"
        "Prototypal and not thoroughly tested!\n"
    )
    warnings.warn(msg)

    def evaluate_spline(t, c, k, xp, nu, extrapolate, out):
        ndim = 1
        c1 = c.reshape(c.shape[:ndim] + (-1,))
        num_c_tr = c1.shape[-1]
        strides_c1 = [stride // c.dtype.itemsize for stride in c.strides]
        indices_k1d = np.unravel_index(np.arange((k+1)**ndim), (k+1,)*ndim)[0][:, None]
        return scpinterp._bspl.evaluate_ndbspline(
            xp[:, None],
            t[None, :],
            np.array([t.size], dtype=np.int32),
            np.array([k], dtype=np.int32),
            np.array([nu], dtype=np.int32),
            extrapolate,
            c.ravel(),
            num_c_tr,
            np.array(strides_c1, dtype=np.intp),
            indices_k1d,
            out,
        )


# ################################################################
# ################################################################
#                       BivariateSplineRect - scipy subclass
# ################################################################


class BivariateSplineRect(scpinterp.BivariateSpline):
    """ Subclass for tofu

    Defined from knots (unique) and deg
    coefs set to 1 by default

    Used self.set_coefs() to update
    """

    def __init__(self, knots0=None, knots1=None, deg=None, shapebs=None):

        assert np.allclose(np.unique(knots0), knots0)
        assert np.allclose(np.unique(knots1), knots1)
        assert deg in [0, 1, 2, 3]

        # get knots pr bs
        self._get_knots_per_bs_for_basis_elements(
            knots0=knots0,
            knots1=knots1,
            deg=deg,
        )

        # full knots with multiplicity
        knots0, nbs0 = _utils_bsplines._get_knots_per_bs(
            knots0, deg=deg, returnas='data', return_unique=True,
        )
        knots1, nbs1 = _utils_bsplines._get_knots_per_bs(
            knots1, deg=deg, returnas='data', return_unique=True,
        )

        coefs = np.ones((nbs0*nbs1,), dtype=float)

        self.__nbs = (nbs0, nbs1)
        self.tck = [knots0, knots1, coefs]
        self.degrees = [deg, deg]

        # shapebs
        self.shapebs = shapebs

    def _get_knots_per_bs_for_basis_elements(
        self,
        knots0=None,
        knots1=None,
        deg=None,
    ):

        # added for details
        knots_per_bs_x0 = _utils_bsplines._get_knots_per_bs(
            knots0, deg=deg, returnas='data',
        )
        knots_per_bs_x1 = _utils_bsplines._get_knots_per_bs(
            knots1, deg=deg, returnas='data',
        )

        self.knots_per_bs_x0 = knots_per_bs_x0
        self.knots_per_bs_x1 = knots_per_bs_x1

        if deg == 0:
            pass
        else:
            knots_per_bs_x0 = np.concatenate(
                (
                    np.tile(knots_per_bs_x0[0, :] - 1, (deg, 1)),
                    knots_per_bs_x0,
                    np.tile(knots_per_bs_x0[-1, :] + 1, (deg, 1)),
                ),
                axis=0,
            )
            knots_per_bs_x1 = np.concatenate(
                (
                    np.tile(knots_per_bs_x1[0, :] - 1, (deg, 1)),
                    knots_per_bs_x1,
                    np.tile(knots_per_bs_x1[-1, :] + 1, (deg, 1)),
                ),
                axis=0,
            )

        self.knots_per_bs_x0_pad = np.asfortranarray(knots_per_bs_x0)
        self.knots_per_bs_x1_pad = np.asfortranarray(knots_per_bs_x1)

    def set_coefs(
        self,
        coefs=None,
        cropbs_neg_flat=None,
    ):

        if coefs.shape == self.shapebs:
            self.tck[2][...] = coefs.ravel()
        elif coefs.shape == (self.nbs,):
            self.tck[2][...] = coefs
        else:
            msg = f"Wrong coefs shape!\nProvided: {coefs.shape}"
            raise Exception(msg)

        # ------------
        # crop and set

        if cropbs_neg_flat is not None:
            self.tck[2][cropbs_neg_flat] = 0.

    def __call__(
        self,
        # interp points
        x0=None,
        x1=None,
        # derivatives
        deriv=None,
        # coefs
        coefs=None,
        axis=None,
        # crop
        crop=None,
        cropbs=None,
        crop_path=None,
        # options
        val_out=None,
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

        # prepare
        val = np.zeros(shape_v, dtype=float)
        cropbs_neg_flat = (~cropbs).ravel() if crop else None

        # interpolate
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

            self.set_coefs(
                coefs=coefs[slic],
                cropbs_neg_flat=cropbs_neg_flat,
            )

            # can be called on any shape of x0, x1?
            val[sliv] = super().__call__(
                x0[slix],
                x1[slix],
                dx=int(deriv[0]),
                dy=int(deriv[1]),
                grid=False,
            )

        # clean out-of-mesh
        if dref_com['ix'] is None and val_out is not False:

            # out of cropped mesh
            if crop_path is not None:
                if x0.ndim == 1:
                    indout = ~crop_path.contains_points(np.array([x0, x1]).T)
                else:
                    indout = ~crop_path.contains_points(
                        np.array([x0.ravel(), x1.ravel()]).T
                    ).reshape(x0.shape)
            else:
                indout = (
                    (x0 < self.tck[0][0]) | (x0 > self.tck[0][-1])
                    | (x1 < self.tck[1][0]) | (x1 > self.tck[1][-1])
                )

            slio = sli_o(indout)
            val[slio] = val_out

        return val

    def ev_details(
        self,
        x0=None,
        x1=None,
        #derivatives
        deriv=None,
        # others
        indbs_tf=None,
        crop=None,
        cropbs=None,
        # for compatibility (unused)
        **kwdargs,
    ):
        """
        indbs_tf = (ar0, ar1)
            tuple of 2 flat arrays of int (for R and Z)
        """

        # -----------
        # prepare

        deg = self.degrees[0]
        nbs = indbs_tf[0].size
        shape = x0.shape
        x0 = np.ascontiguousarray(x0.ravel(), dtype=float)
        x1 = np.ascontiguousarray(x1.ravel(), dtype=float)
        coef = np.zeros((deg + 4, 1), dtype=float)
        coef[deg] = 1.
        out1 = np.full((x0.size, 1), np.nan)

        # -----------
        # compute

        val = np.zeros(tuple(np.r_[x0.shape, nbs]))
        indtot = np.arange(0, nbs)

        iz_u = np.unique(indbs_tf[1])

        for iz in iz_u:

            evaluate_spline(
                self.knots_per_bs_x1_pad[:, iz],
                coef,
                self.degrees[1],
                x1,
                deriv[0],
                False,
                out1,
            )

            indok1 = ~np.isnan(out1)
            if not np.any(indok1):
                continue
            indok0 = np.copy(indok1)

            ind0 = indbs_tf[1] == iz
            i0 = indbs_tf[0][ind0]
            for ii, ii0 in enumerate(i0):

                msg = (
                    f'\t bspline {indtot[ind0][ii]} / {nbs}   ({x0.size} pts)'
                )
                print(msg, end='\r', flush=True)

                if ii > 0:
                    indok0[...] = indok1

                out0 = np.full((indok1.sum(), 1), np.nan)

                evaluate_spline(
                    self.knots_per_bs_x0_pad[:, ii0],
                    coef,
                    self.degrees[0],
                    x0[indok1[:, 0]],
                    deriv[1],
                    False,
                    out0,
                )

                ixok = ~np.isnan(out0)
                if not np.any(ixok):
                    continue

                indok0[indok1] = ixok[:, 0]
                val[indok0[:, 0], indtot[ind0][ii]] = (out0[ixok]*out1[indok0])

        if shape != x0.shape:
            val = np.reshape(val, tuple(np.r_[shape, -1]))

        return val

    def get_overlap(self):
        return _get_overlap(
            deg=self.degrees[0],
            knots0=self.knots_per_bs_x0,
            knots1=self.knots_per_bs_x1,
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
        return _class02_bsplines_operators_rect.get_mesh2dRect_operators(
            deg=self.degrees[0],
            operator=operator,
            geometry=geometry,
            knotsx_mult=self.tck[0],
            knotsy_mult=self.tck[1],
            knotsx_per_bs=self.knots_per_bs_x0,
            knotsy_per_bs=self.knots_per_bs_x1,
            overlap=self.get_overlap(),
            cropbs_flat=cropbs_flat,
            # specific to deg = 0
            cropbs=cropbs,
            centered=centered,
            # to return gradR, gradZ, for D1N2 deg 0, for tomotok
            returnas_element=returnas_element,
        )


# ################################################################
# ################################################################
#                       Mesh2DRect - bsplines - overlap
# ################################################################


def _get_overlap(
    deg=None,
    knots0=None,
    knots1=None,
    shapebs=None,
):
    # nb of overlapping, inc. itself in 1d
    nbs0, nbs1 = shapebs
    ind00 = np.tile(np.arange(0, nbs0), nbs1)
    ind10 = np.repeat(np.arange(0, nbs1), nbs0)

    # complete
    ntot = 2*deg + 1

    add0= np.tile(np.arange(-deg, deg+1), ntot)
    add1 = np.repeat(np.arange(-deg, deg+1), ntot)

    inter0 = ind00[None, :] + add0[:, None]
    inter1 = ind10[None, :] + add1[:, None]

    # purge
    inter = inter0 + inter1*nbs0
    indneg = (
        (inter0 < 0) | (inter0 >= nbs0) | (inter1 < 0) | (inter1 >= nbs1)
    )
    inter[indneg] = -1

    return inter


# #############################################################################
# #############################################################################
#                           Mesh2DRect - bsplines
# #############################################################################


def get_bs2d_x01(deg=None, knots0=None, knots1=None):

    # ----------------
    # get knots per bspline, nb of bsplines...

    knots_per_bs_x0 = _utils_bsplines._get_knots_per_bs(
        knots0, deg=deg, returnas='data',
    )
    knots_per_bs_x1 = _utils_bsplines._get_knots_per_bs(
        knots1, deg=deg, returnas='data',
    )
    # nbkbs = knots_per_bs_R.shape[0]
    shapebs = (knots_per_bs_x0.shape[1], knots_per_bs_x1.shape[1])

    # ----------------
    # get centers of bsplines

    x0bs_apex = _utils_bsplines._get_apex_per_bs(
        knots=knots0,
        knots_per_bs=knots_per_bs_x0,
        deg=deg
    )
    x1bs_apex = _utils_bsplines._get_apex_per_bs(
        knots=knots1,
        knots_per_bs=knots_per_bs_x1,
        deg=deg
    )
    return shapebs, x0bs_apex, x1bs_apex, knots_per_bs_x0, knots_per_bs_x1


def get_bs_class(
    deg=None,
    knots0=None,
    knots1=None,
    shapebs=None,
    # knots_per_bs_R=None,
    # knots_per_bs_Z=None,
):

    # ----------------
    # Define functions

    return BivariateSplineRect(
        knots0=knots0,
        knots1=knots1,
        deg=deg,
        shapebs=shapebs,
    )