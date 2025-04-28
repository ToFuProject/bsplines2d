# -*- coding: utf-8 -*-


# Common
import datastock as ds


# local
from ._class01_Mesh2D import Mesh2D as Previous
from . import _class02_checks as _checks


__all__ = ['Mesh3D']


# #############################################
# #############################################
#               Class
# #############################################


class Mesh3D(Previous):

    _which_mesh3d = 'mesh3d'

    _show_in_summary = 'all'
    _dshow = dict(ds.DataStock._dshow)
    _dshow.update({
        _which_mesh3d: [
            'nd',
            'type',
            'mesh2d',
            'mesh1d_angle',
        ],
    })

    def add_mesh_3d_cyl(
        self,
        key=None,
        key_mesh2d=None,
        key_mesh1d_angle=None,
        # direct addition of bsplines
        deg=None,
        # additional attributes
        **kwdargs,
    ):
        """ Add an 3d mesh by key

        The mesh is defined by:
            - key_mesh2d: a pre-existing 2d mesh (rect or tri)
            - key_mesh1d_angle: a pre-existing 1d angle mesh

        If deg is provided, immediately adds a bsplines

        Example:
        --------
                >>> import bsplines2d as bs2
                >>> coll = bs2.Collection()
                # add 1d angle mesh
                >>> coll.add_mesh_1d(
                    'm1dangle',
                    knots=np.linspace(-np.pi, np.pi, 10),
                    units='rad',
                )
                # add 2d mesh
                >>> coll.add_mesh_2d_rect(
                    'm2d',
                    knots0=np.linspace(1, 2, 10),
                    knots1=np.linspace(-1, 1, 10),
                )
                # add 3d mesh
                >>> coll.add_mesh_3d_cyl(
                    key='m3d',
                    key_mesh2d='m2d',
                    key_mesh1d_angle='m1dangle',
                )
        """

        # check input data and get input dicts
        key, dref, ddata, dobj = _checks.main(
            coll=self,
            # key
            key=key,
            # mesh2d
            key_mesh2d=key_mesh2d,
            # mesh angle
            key_mesh1d_angle=key_mesh1d_angle,
            # additional attributes
            **kwdargs,
        )

        # update dicts
        self.update(dref=dref, ddata=ddata, dobj=dobj)

        # optional bspline
        if deg is not None:
            self.add_bsplines(key=key, deg=deg)

    # -----------
    # plot
    # -----------

    def plot_mesh(
        self,
    ):
        """ Plot the selected 3d mesh
        """

        return _plot.main(
            coll=self,
        )



