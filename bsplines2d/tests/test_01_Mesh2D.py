"""
This module contains tests for tofu.geom in its structured version
"""

# Built-in
import itertools as itt


# Standard
import numpy as np


# specific
from . import test_input
from .. import _class01_checks_2d_rect as _checks
from .._class03_Bins import Bins as BSplines2D


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
#     checking routines
#
#######################################################


class Test00_check_routines():

    def test00_mesh2DRect_X_check(self):

        lx = [[1, 2], [1, 2, 3, 4]]
        lres = [None, 10, 0.1, [0.1, 0.2], [0.1, 0.2, 0.3, 0.1]]

        for comb in itt.product(lx, lres):
            if hasattr(lres, '__iter__') and len(lres) != len(lx):
                continue
            x, res, ind = _checks._mesh2DRect_X_check(
                x=[1, 2, 3, 4],
                res=10,
            )
            if hasattr(lres, '__iter__'):
                assert x.size == np.unique(x).size == res.size + 1


#######################################################
#
#     Fixed meshes (1d, rect and tri, no submesh)
#
#######################################################


class Test01_BSplines2D():

    @classmethod
    def setup_class(cls):
        cls.bs = BSplines2D()

    @classmethod
    def setup_method(self):
        pass

    @classmethod
    def teardown_method(self):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    # #################
    #  mesh 1d (fixed)
    # #################

    def test01_mesh1d_from_knots_uniform(self):
        test_input._add_1d_knots_uniform(self.bs)

    def test02_mesh1d_from_knots_variable(self):
        test_input._add_1d_knots_variable(self.bs)

    # ##############
    #  mesh 2d rect
    # ##############

    def test03_add_mesh_rect_uniform(self):
        test_input._add_rect_uniform(self.bs)

    def test04_add_mesh_rect_variable(self):
        test_input._add_rect_variable(self.bs)

    def test05_add_mesh_rect_variable_crop(self):
        test_input._add_rect_variable_crop(self.bs)

    def test06_add_mesh_rect_variable_crop(self):
        test_input._add_rect_variable_crop(self.bs)

    def test07_add_mesh_rect_variable_crop(self):
        test_input._add_rect_variable_crop(self.bs)

    # ##############
    #  mesh 2d tri
    # ##############

    def test08_add_mesh_tri_ntri1(self):
        test_input._add_tri_ntri1(self.bs)

    def test09_add_mesh_tri_ntri2(self):
        test_input._add_tri_ntri2(self.bs)

    def test10_add_mesh_tri_delaunay(self):
        test_input._add_tri_delaunay(self.bs)

    # ##############
    #  select
    # ##############

    def test11_select_mesh_element_1d(self):
        test_input._select_mesh_elements(self.bs, nd='1d', kind=None)

    def test12_select_mesh_element_rect(self):
        test_input._select_mesh_elements(self.bs, nd='2d', kind='rect')

    def test13_select_mesh_element_tri(self):
        test_input._select_mesh_elements(self.bs, nd='2d', kind='tri')

    # ##############
    #  sample
    # ##############

    def test14_sample_mesh_1d(self):
        test_input._sample_mesh(self.bs, nd='1d', kind=None)

    def test15_sample_mesh_rect(self):
        test_input._sample_mesh(self.bs, nd='2d', kind='rect')

    def test16_sample_mesh_tri(self):
        test_input._sample_mesh(self.bs, nd='2d', kind='tri')

    # ##############
    #  plot
    # ##############

    def test17_plot_mesh_1d(self):
        test_input._plot_mesh(self.bs, nd='1d', kind=None)

    def test18_plot_mesh_rect(self):
        test_input._plot_mesh(self.bs, nd='2d', kind='rect')

    def test19_plot_mesh_tri(self):
        test_input._plot_mesh(self.bs, nd='2d', kind='tri')

    # ##############
    #  add bsplines
    # ##############

    def test20_add_bsplines_1d(self):
        test_input._add_bsplines(self.bs, nd='1d')

    def test21_add_bsplines_2d_rect(self):
        test_input._add_bsplines(self.bs, kind='rect')

    def test22_add_bsplines_2d_tri(self):
        test_input._add_bsplines(self.bs, kind='tri')

    # ##############
    #  select bsplines
    # ##############

    def test23_select_bsplines_1d(self):
        test_input._select_bsplines(self.bs, nd='1d', kind=None)

    def test24_select_bsplines_2d_rect(self):
        test_input._select_bsplines(self.bs, nd='2d', kind='rect')

    def test25_select_bsplines_2d_tri(self):
        test_input._select_bsplines(self.bs, nd='2d', kind='tri')

    # ###############
    #  add data vs bs
    # ###############

    def test26_add_data_1bs_fix_1d(self, remove=True):
        test_input._add_data_1bs_fix(self.bs, nd='1d', kind=None, remove=remove)

    def test27_add_data_1bs_fix_2d_rect(self, remove=True):
        test_input._add_data_1bs_fix(self.bs, nd='2d', kind='rect', remove=remove)

    def test28_add_data_1bs_fix_2d_tri(self, remove=True):
        test_input._add_data_1bs_fix(self.bs, nd='2d', kind='tri', remove=remove)

    def test29_add_data_1bs_arrays_1d(self, remove=False):
        test_input._add_data_1bs_arrays(self.bs, nd='1d', kind=None, remove=remove)

    def test30_add_data_1bs_arrays_2d_rect(self, remove=False):
        test_input._add_data_1bs_arrays(self.bs, nd='2d', kind='rect', remove=remove)

    def test31_add_data_1bs_arrays_2d_tri(self, remove=False):
        test_input._add_data_1bs_arrays(self.bs, nd='2d', kind='tri', remove=remove)

    def test32_add_data_multibs_arrays(self, remove=False):
        test_input._add_data_multibs_arrays(self.bs, remove=remove)

    # ##############
    # interp bs
    # ##############

    def test33_interpolate_bsplines_1d(self):
        test_input._interpolate(self.bs, nd='1d', kind=None, details=False)

    def test34_interpolate_bsplines_1d_details(self):
        test_input._interpolate(self.bs, nd='1d', kind=None, details=True)

    def test35_interpolate_bsplines_2d_rect(self):
        test_input._interpolate(self.bs, nd='2d', kind='rect', details=False)

    def test36_interpolate_bsplines_2d_rect_details(self):
        test_input._interpolate(self.bs, nd='2d', kind='rect', details=True)

    def test37_interpolate_bsplines_2d_tri(self):
        test_input._interpolate(self.bs, nd='2d', kind='tri', details=False)

    def test38_interpolate_bsplines_2d_tri_details(self):
        test_input._interpolate(self.bs, nd='2d', kind='tri', details=True)

    # ##############
    # binning 1d
    # ##############

    def test39_binning_1d(self):
        test_input._bin_bs(self.bs, nd='1d', kind=None)

    # ##############
    # plot bsplines
    # ##############

    def test40_plot_bsplines_1d(self):
        pass

    def test41_plot_bsplines_2d_rect(self):
        pass

    def test42_plot_bsplines_2d_tri(self):
        pass

    # ####################
    # add mesh with subkey
    # ####################

    def test43_add_mesh_1d_subkey_1d(self):
        test_input._add_mesh_1d_subkey(self.bs, nd='1d', kind=None)

    def test44_add_mesh_1d_subkey_rect(self):
        test_input._add_mesh_1d_subkey(self.bs, nd='2d', kind='rect')

    def test45_add_mesh_1d_subkey_tri(self):
        test_input._add_mesh_1d_subkey(self.bs, nd='2d', kind='tri')

    def test46_add_mesh_2d_rect_subkey_rect(self):
        test_input._add_mesh_2d_rect_subkey(self.bs, nd='2d', kind='rect')

    def test47_add_mesh_2d_rect_subkey_tri(self):
        test_input._add_mesh_2d_rect_subkey(self.bs, nd='2d', kind='tri')

    def test48_add_mesh_2d_rect_var_subkey_rect(self):
        test_input._add_mesh_2d_rect_subkey(self.bs, nd='2d', kind='rect')

    # ################################
    # add bsplines on mesh with subkey
    # ################################

    def test49_add_bsplines_subkey(self):
        test_input._add_bsplines(self.bs, subkey=True)

    # ################################
    # add data on bsplines with subkey
    # ################################

    def test50_add_data_subkey(self):
        test_input._add_data_multibs_arrays(
            self.bs,
            nd=None,
            kind=None,
            subbs=True,
            remove=False,
        )

    # ################################
    # interpolate data with subkey
    # ################################

    def test51_interpolate_subkey_1d(self):
        test_input._interpolate(
            self.bs,
            nd='1d',
            kind=None,
            submesh=True,
        )

    # def test48_interpolate_data_subkey_from_subkey(self):
        # test_input._interpolate_from_subkey(
            # self.bs,
        # )

    # ################################
    # interpolate data with subkey
    # ################################

    def test53_plot_as_profile2d(self):
        test_input._plot_as_profile2d(
            self.bs,
            nd='2d',
            kind=None,
        )

    def test54_plot_as_profile2d_compare(self):
        test_input._plot_as_profile2d_compare(
            self.bs,
            nd='2d',
            kind=None,
        )

    # ################################
    # operators
    # ################################

    def test55_operators_1d(self):
        test_input._get_operators(
            self.bs,
            nd='1d',
            kind=None,
        )

    def test56_operators_2d_rect(self):
        test_input._get_operators(
            self.bs,
            nd='2d',
            kind='rect',
        )

    def test57_operators_2d_tri(self):
        pass

    def test58_operators_1d_subkey(self):
        pass
