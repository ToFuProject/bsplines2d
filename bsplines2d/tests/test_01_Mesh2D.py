"""
This module contains tests for tofu.geom in its structured version
"""

# Built-in
import os
import itertools as itt


# Standard
import numpy as np


# specific
from . import test_input
from .. import _class01_checks_2d_rect as _checks
from .._class04_Bins import Bins as Collection
from .. import _saveload


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


class Test01_Collection():

    @classmethod
    def setup_class(cls):
        cls.bs = Collection()

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

    def test03_mesh1d_from_knots_angle(self):
        test_input._add_1d_knots_angle(self.bs)

    # ##############
    #  mesh 2d rect
    # ##############

    def test04_add_mesh_rect_uniform(self):
        test_input._add_rect_uniform(self.bs)

    def test05_add_mesh_rect_variable(self):
        test_input._add_rect_variable(self.bs)

    def test06_add_mesh_rect_variable_crop(self):
        test_input._add_rect_variable_crop(self.bs)

    def test07_add_mesh_rect_crop_from_knots(self):
        test_input._add_rect_crop_from_knots(self.bs)

    def test08_add_mesh_rect_variable_crop_from_knots(self):
        test_input._add_rect_variable_crop_from_knots(self.bs)

    # ##############
    #  mesh 2d tri
    # ##############

    def test09_add_mesh_tri_ntri1(self):
        test_input._add_tri_ntri1(self.bs)

    def test10_add_mesh_tri_ntri2(self):
        test_input._add_tri_ntri2(self.bs)

    def test11_add_mesh_tri_delaunay(self):
        test_input._add_tri_delaunay(self.bs)

    # ##############
    #  mesh 3d cyl
    # ##############

    def test12_add_mesh_3d_cyl(self):
        test_input._add_3d_cyl(self.bs)

    # ##############
    #  select
    # ##############

    def test13_select_mesh_element_1d(self):
        test_input._select_mesh_elements(self.bs, nd='1d', kind=None)

    def test14_select_mesh_element_rect(self):
        test_input._select_mesh_elements(self.bs, nd='2d', kind='rect')

    def test15_select_mesh_element_tri(self):
        test_input._select_mesh_elements(self.bs, nd='2d', kind='tri')

    # ##############
    #  mesh outline
    # ##############

    def test16_get_mesh_outline_rect(self):
        test_input._get_mesh_outline_rect(self.bs)

    def test17_get_mesh_outline_tri(self):
        test_input._get_mesh_outline_tri(self.bs)

    # ##############
    #  sample
    # ##############

    def test18_sample_mesh_1d(self):
        test_input._sample_mesh(self.bs, nd='1d', kind=None)

    def test19_sample_mesh_rect(self):
        test_input._sample_mesh(self.bs, nd='2d', kind='rect')

    def test20_sample_mesh_tri(self):
        test_input._sample_mesh(self.bs, nd='2d', kind='tri')

    # ##############
    #  plot
    # ##############

    def test21_plot_mesh_1d(self):
        test_input._plot_mesh(self.bs, nd='1d', kind=None)

    def test22_plot_mesh_rect(self):
        test_input._plot_mesh(self.bs, nd='2d', kind='rect')

    def test23_plot_mesh_tri(self):
        test_input._plot_mesh(self.bs, nd='2d', kind='tri')

    # ##############
    #  add bsplines
    # ##############

    def test24_add_bsplines_1d(self):
        test_input._add_bsplines(self.bs, nd='1d')

    def test25_add_bsplines_2d_rect(self):
        test_input._add_bsplines(self.bs, kind='rect')

    def test26_add_bsplines_2d_tri(self):
        test_input._add_bsplines(self.bs, kind='tri')

    # ##############
    #  select bsplines
    # ##############

    def test27_select_bsplines_1d(self):
        test_input._select_bsplines(self.bs, nd='1d', kind=None)

    def test28_select_bsplines_2d_rect(self):
        test_input._select_bsplines(self.bs, nd='2d', kind='rect')

    def test29_select_bsplines_2d_tri(self):
        test_input._select_bsplines(self.bs, nd='2d', kind='tri')

    # ###############
    #  add data vs bs
    # ###############

    def test30_add_data_1bs_fix_1d(self, remove=True):
        test_input._add_data_1bs_fix(self.bs, nd='1d', kind=None, remove=remove)

    def test31_add_data_1bs_fix_2d_rect(self, remove=True):
        test_input._add_data_1bs_fix(self.bs, nd='2d', kind='rect', remove=remove)

    def test32_add_data_1bs_fix_2d_tri(self, remove=True):
        test_input._add_data_1bs_fix(self.bs, nd='2d', kind='tri', remove=remove)

    def test33_add_data_1bs_arrays_1d(self, remove=False):
        test_input._add_data_1bs_arrays(self.bs, nd='1d', kind=None, remove=remove)

    def test34_add_data_1bs_arrays_2d_rect(self, remove=False):
        test_input._add_data_1bs_arrays(self.bs, nd='2d', kind='rect', remove=remove)

    def test35_add_data_1bs_arrays_2d_tri(self, remove=False):
        test_input._add_data_1bs_arrays(self.bs, nd='2d', kind='tri', remove=remove)

    def test36_add_data_multibs_arrays(self, remove=False):
        test_input._add_data_multibs_arrays(self.bs, remove=remove)

    # ##############
    # interp bs
    # ##############

    def test37_interpolate_bsplines_1d(self):
        test_input._interpolate(self.bs, nd='1d', kind=None, details=False)

    def test38_interpolate_bsplines_1d_details(self):
        test_input._interpolate(self.bs, nd='1d', kind=None, details=True)

    def test39_interpolate_bsplines_2d_rect(self):
        test_input._interpolate(self.bs, nd='2d', kind='rect', details=False)

    def test40_interpolate_bsplines_2d_rect_details(self):
        test_input._interpolate(self.bs, nd='2d', kind='rect', details=True)

    def test41_interpolate_bsplines_2d_tri(self):
        test_input._interpolate(self.bs, nd='2d', kind='tri', details=False)

    def test42_interpolate_bsplines_2d_tri_details(self):
        test_input._interpolate(self.bs, nd='2d', kind='tri', details=True)

    # ##############
    # binning 1d
    # ##############

    def test43_binning_1d(self):
        test_input._bin_bs(self.bs, nd='1d', kind=None)

    # ##############
    # plot bsplines
    # ##############

    def test44_plot_bsplines_1d(self):
        pass

    def test45_plot_bsplines_2d_rect(self):
        pass

    def test46_plot_bsplines_2d_tri(self):
        pass

    # ####################
    # add mesh with subkey
    # ####################

    def test47_add_mesh_1d_subkey_1d(self):
        test_input._add_mesh_1d_subkey(self.bs, nd='1d', kind=None)

    def test48_add_mesh_1d_subkey_rect(self):
        test_input._add_mesh_1d_subkey(self.bs, nd='2d', kind='rect')

    def test49_add_mesh_1d_subkey_tri(self):
        test_input._add_mesh_1d_subkey(self.bs, nd='2d', kind='tri')

    def test50_add_mesh_2d_rect_subkey_rect(self):
        test_input._add_mesh_2d_rect_subkey(self.bs, nd='2d', kind='rect')

    def test51_add_mesh_2d_rect_subkey_tri(self):
        test_input._add_mesh_2d_rect_subkey(self.bs, nd='2d', kind='tri')

    def test52_add_mesh_2d_rect_var_subkey_rect(self):
        test_input._add_mesh_2d_rect_subkey(self.bs, nd='2d', kind='rect')

    # ################################
    # add bsplines on mesh with subkey
    # ################################

    def test53_add_bsplines_subkey(self):
        test_input._add_bsplines(self.bs, subkey=True)

    # ################################
    # add data on bsplines with subkey
    # ################################

    def test54_add_data_subkey(self):
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

    def test55_interpolate_subkey_1d(self):
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

    def test56_plot_as_profile2d(self):
        test_input._plot_as_profile2d(
            self.bs,
            nd='2d',
            kind=None,
        )

    def test57_plot_as_profile2d_compare(self):
        test_input._plot_as_profile2d_compare(
            self.bs,
            nd='2d',
            kind=None,
        )

    # ################################
    # operators
    # ################################

    def test58_operators_1d(self):
        test_input._get_operators(
            self.bs,
            nd='1d',
            kind=None,
        )

    def test59_operators_2d_rect(self):
        test_input._get_operators(
            self.bs,
            nd='2d',
            kind='rect',
        )

    def test60_operators_2d_tri(self):
        pass

    def test61_operators_1d_subkey(self):
        pass

    # ################################
    # saving / loading
    # ################################

    def test62_saveload_equal(self):

        # save
        pfe = self.bs.save(return_pfe=True)

        # load
        out = _saveload.load(pfe)

        # remove file
        os.remove(pfe)

        # equal
        assert self.bs == out

    def test63_saveload_coll(self):

        # save
        pfe = self.bs.save(return_pfe=True)

        # load
        coll = Collection()
        coll = _saveload.load(pfe, coll=coll)

        # remove file
        os.remove(pfe)

        # equal
        assert self.bs == coll
