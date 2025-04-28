# -*- coding: utf-8 -*-


# Common
import datastock as ds


# local
from ._class01_Mesh2D import Mesh2D as Previous
from . import _generic_mesh
from . import _class01_checks_1d as _checks_1d
from . import _class01_checks_2d_rect as _checks_2d_rect
from . import _class01_checks_2d_tri as _checks_2d_tri
# from . import _class01_checks_2d_polar as _checks_2d_polar


from . import _class01_cropping as _cropping
from . import _class01_select as _select
from . import _class01_sample as _sample
from . import _class01_outline as _outline
from . import _class01_plot as _plot


__all__ = ['Mesh3D']


# #############################################################################
# #############################################################################
#
# #############################################################################


class Mesh3D(Previous):

    _which_mesh3d = 'mesh3d'

    _show_in_summary = 'all'
    _dshow = dict(ds.DataStock._dshow)

    def add_mesh_3d_cyl(
        self,
        key=None,
        key_mesh2d=None,
        knots_phi=None,
        # direct addition of bsplines
        deg=None,
        # additional attributes
        **kwdargs,
    ):

        """ Add an 3d mesh by key

        The mesh is defined by:
            - key_mesh2d: a pre-existing 2d mesh (rect or tri)
            - knots_phi: a vector of knots on angle phi (toroidal)

        If deg is provided, immediately adds a bsplines

        Example:
        --------
                >>> import bsplines2d as bs2
                >>> coll = bs2.Collection()
                >>> coll.add_mesh_2d_rect(
                    'm2d',
                    knots0=np.linspace(1, 2, 10),
                    knots1=np.linspace(-1, 1, 10),
                )
                >>> coll.add_mesh_3d_cyl(
                    key='m3d',
                    key_mesh2d='m2d',
                    knots_phi=np.linspace(-np.pi, np.pi, 6),
                )

        """

        # check input data and get input dicts
        key, dref, ddata, dobj = _checks_1d.check(
            coll=self,
            # key
            key=key,
            # mesh knots
            knots=knots,
            uniform=uniform,
            # defined from pre-existing bsplines
            subkey=subkey,
            # additional attributes
            **kwdargs,
        )

        # update dicts
        self.update(dref=dref, ddata=ddata, dobj=dobj)

        # optional bspline
        if deg is not None:
            self.add_bsplines(key=key, deg=deg)

