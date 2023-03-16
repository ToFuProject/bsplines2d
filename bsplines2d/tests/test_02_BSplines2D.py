"""
This module contains tests for tofu.geom in its structured version
"""

# Built-in
import os
import itertools as itt
import warnings


# Standard
import numpy as np
import matplotlib.pyplot as plt


# specific
from . import test_input


#######################################################
#
#     Setup and Teardown
#
#######################################################


def setup_module():
    pass

def teardown_module():
    pass


#######################################################
#
#     object mesh2D
#
#######################################################


"""
class Test01_BSplines():

    @classmethod
    def setup_class(cls):
        pass

    def setup_method(self):
        bsplines = tfd.Collection()

        # add rect mesh
        _add_rect_uniform(bsplines)
        _add_rect_variable(bsplines)
        _add_rect_variable_crop(bsplines)

        # add tri mesh
        _add_tri_ntri1(bsplines)
        _add_tri_ntri2(bsplines)

        # add bsplines
        _add_bsplines(bsplines)

        # add polar mesh
        _add_polar1(bsplines)
        _add_polar2(bsplines)

        # add bsplines for polar meshes
        _add_bsplines(bsplines, kind=['polar'])

        # Add polar with variable poloidal discretization
        _add_polar2(bsplines, key='m7')
        _add_bsplines(
            bsplines,
            key=['m7'],
            angle=np.pi*np.r_[-3./4., -1/4, 0, 1/4, 3/4],
        )

        # store
        self.obj = bsplines
        self.lm = list(bsplines.dobj['mesh'].keys())
        self.lbs = list(bsplines.dobj['bsplines'].keys())

    def teardown_method(self):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def test04_select_bsplines(self):

        lind0 = [None, ([0, 2], [0, 4]), [0, 2, 4], ([0, 2, 4], [0, 2, 3])]
        lind1 = [None, [1], 1, [0, 1]]
        for ii, k0 in enumerate(self.lbs):

            km = self.obj.dobj['bsplines'][k0]['mesh']
            if len(self.obj.dobj['bsplines'][k0]['shape']) == 2:
                lind = lind0
            else:
                lind = lind1

            if self.obj.dobj['mesh'][km]['type'] == 'polar':
                return_cents = False
                return_knots = False
            else:
                return_cents = None if ii == 1 else bool(ii%3)
                return_knots = None if ii == 2 else bool(ii%2)

            out = self.obj.select_bsplines(
                key=k0,
                ind=lind[ii%len(lind)],
                returnas='ind' if ii%3 == 0 else 'data',
                return_cents=return_cents,
                return_knots=return_knots,
            )

    def test08_interpolate_profile2d_details_vs_sum(self):

        x = np.linspace(2.2, 2.8, 5)
        y = np.linspace(-0.5, 0.5, 5)
        x = np.tile(x, (y.size, 1))
        y = np.tile(y, (x.shape[1], 1)).T

        for ii, k0 in enumerate(self.lbs):

            keym = self.obj.dobj['bsplines'][k0]['mesh']
            mtype = self.obj.dobj['mesh'][keym]['type']

            val, t, ref = self.obj.interpolate_profile2d(
                key=k0,
                R=x,
                Z=y,
                coefs=None,
                indbs=None,
                indt=None,
                grid=False,
                details=True,
                reshape=True,
                res=None,
                crop=None,
                nan0=ii % 2 == 0,
                imshow=False,
                return_params=False,
            )

            crop = self.obj.dobj['bsplines'][k0].get('crop', False)
            nbs = np.prod(self.obj.dobj['bsplines'][k0]['shape'])
            if isinstance(crop, str):
                nbs = self.obj.ddata[crop]['data'].sum()

            vshap0 = tuple(np.r_[x.shape, nbs])
            if mtype == 'polar':
                # radius2d can be time-dependent => additional dimension
                vshap = val.shape[-len(vshap0):]
            else:
                vshap = val.shape
            assert vshap == vshap0, val.shape

            val_sum, t, ref = self.obj.interpolate_profile2d(
                key=k0,
                R=x,
                Z=y,
                coefs=None,
                indbs=None,
                indt=None,
                grid=False,
                details=False,
                reshape=True,
                res=None,
                crop=None,
                nan0=False,
                val_out=0.,
                imshow=False,
                return_params=False,
            )

            if mtype == 'polar':
                # radius2d can be time-dependent => additional dimension
                vshap_sum = val_sum.shape[-len(x.shape):]
            else:
                vshap_sum = val_sum.shape
            assert vshap_sum == x.shape, val_sum.shape
            assert (val.ndim == x.ndim + 2) == (val_sum.ndim == x.ndim + 1), [val.shape, val_sum.shape]

            indok = np.isfinite(val_sum)
            indok[indok] = val_sum[indok] != 0

            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # Does not work for rect mesh
            # because of knots padding used in func_details
            # Due to scpinterp._bspl.evaluate_spline()...
            if mtype in ['tri', 'polar']:   # To be debugged
                assert np.allclose(
                    val_sum[indok],
                    np.nansum(val, axis=-1)[indok],
                    equal_nan=True,
                )
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    def test10_plot_bsplines(self):

        li0 = [None, ([1, 2], [2, 1]), (1, 1), [1, 2, 4]]
        li1 = [None, [1, 2], (1, 1), [1, 2, 4]]
        for ii, k0 in enumerate(self.lbs):

            km = self.obj.dobj['bsplines'][k0]['mesh']
            if len(self.obj.dobj['mesh'][km]['shape-c']) == 2:
                li = li0
            else:
                li = li1

            if self.obj.dobj['mesh'][km]['type'] == 'polar':
                plot_mesh = False
            else:
                plot_mesh = True

            dax = self.obj.plot_bsplines(
                key=k0,
                indbs=li[ii%len(li)],
                knots=bool(ii%3),
                cents=bool(ii%2),
                res=0.05,
                plot_mesh=plot_mesh,
            )
            plt.close('all')

    def test12_add_bsplines_operator(self):
        lkey = ['m0-bs0', 'm1-bs1', 'm2-bs2']
        lop = ['D0N1', 'D0N2', 'D1N2', 'D2N2']
        lgeom = ['linear', 'toroidal']
        lcrop = [False, True]

        dfail = {}
        for ii, k0 in enumerate(self.lbs):

            km = self.obj.dobj['bsplines'][k0]['mesh']
            if self.obj.dobj['mesh'][km]['type'] == 'tri':
                continue
            elif self.obj.dobj['mesh'][km]['type'] == 'polar':
                continue

            for comb in itt.product(lop, lgeom, lcrop):
                deg = self.obj.dobj['bsplines'][k0]['deg']

                if deg == 3 and comb[0] in ['D0N1', 'D0N2', 'D1N2', 'D2N2']:
                    continue

                # only test exact operators
                if int(comb[0][1]) > deg:
                    # except deg = 0 D1N2
                    if deg == 0 and comb[0] == 'D1N2':
                        pass
                    else:
                        continue
                try:
                    self.obj.add_bsplines_operator(
                        key=k0,
                        operator=comb[0],
                        geometry=comb[1],
                        crop=comb[2],
                    )
                except Exception as err:
                    dfail[k0] = (
                        f"key {k0}, op '{comb[0]}', geom '{comb[1]}': "
                        + str(err)
                    )

        # Raise error if any fail
        if len(dfail) > 0:
            lstr = [f'\t- {k0}: {v0}' for k0, v0 in dfail.items()]
            msg = (
                "The following operators failed:\n"
                + "\n".join(lstr)
            )
            raise Exception(msg)
"""
