# -*- coding: utf-8 -*-


# Common
import datastock as ds


# local
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


__all__ = ['Mesh2D']


# #############################################################################
# #############################################################################
#
# #############################################################################


class Mesh2D(ds.DataStock):

    _which_mesh = 'mesh'
    _which_bsplines = 'bsplines'

    _show_in_summary = 'all'
    _dshow = dict(ds.DataStock._dshow)
    _dshow.update({
        _which_mesh: [
            'nd',
            'type',
            'variable',
            'crop',
            'crop_thresh',
            'ntri',
            'knots',
            'cents',
            'shape_k',
            'shape_c',
            'ind',
            'subkey',
            f'nb {_which_bsplines}',
        ],
    })

    def add_mesh_1d(
        self,
        key=None,
        knots=None,
        uniform=None,
        # defined from pre-existing bsplines
        subkey=None,
        # direct addition of bsplines
        deg=None,
        # additional attributes
        **kwdargs,
    ):

        """ Add an 1d mesh by key

        The mesh is defined by a strictly increasing vector of knots (edges)
        The knots can be defined from a pre-existing 2d bsplines.

        If deg is provided, immediately adds a bsplines

        Example:
        --------
                >>> from bsplines2d import BSplines2D
                >>> bs = BSplines2D()
                >>> bs.add_mesh_1d(key='m0', knots=np.linspace(0, 10, 11))

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

    def add_mesh_2d_rect(
        self,
        # rectangular mesh
        key=None,
        # knots
        knots0=None,
        knots1=None,
        domain=None,
        res=None,
        uniform0=None,
        uniform1=None,
        # cropping
        crop_poly=None,
        thresh_in=None,
        remove_isolated=None,
        # defined from pre-existing bsplines
        subkey0=None,
        subkey1=None,
        # direct addition of bsplines
        deg=None,
        **kwdargs,
    ):
        """ Add a mesh by key and domain / resolution

        Can create a 2d rectangular or triangular mesh:
            - knots0: sequence of unique floats
            - knots1: sequence of unique floats
            # or
            - domain: a list of 2 iterables of len() >= 2
            - res: float


        Can optionally be cropped by a closed polygon crop_poly, that can be:
            - a (2, N) np.narray of (R, Z) coordinates
            - a tuple (Config, key_struct) to designate a struct poly

        Args thresh_in and remove_isolated control the level of cropping:
            - thresh_in:
            - remove_isolated:

        If deg is provided, immediately adds a bsplines

        """

        # get domain, poly from crop_poly
        if crop_poly is not None:
            domain, poly = _checks_2d_rect._from_croppoly(
                crop_poly=crop_poly,
                domain=domain,
            )
        else:
            poly = None

        # check input data and get input dicts
        key, dref, ddata, dobj = _checks_2d_rect.check(
            coll=self,
            # key
            key=key,
            # knots
            knots0=knots0,
            knots1=knots1,
            # rectangular
            domain=domain,
            res=res,
            uniform0=uniform0,
            uniform1=uniform1,
            # defined from pre-existing bsplines
            subkey0=subkey0,
            subkey1=subkey1,
            # attributes
            **kwdargs,
        )

        # update dicts
        self.update(dref=dref, ddata=ddata, dobj=dobj)

        # optional bspline
        if deg is not None:
            self.add_bsplines(key=key, deg=deg)

        # optional cropping
        if poly is not None:
            self.crop(
                key=key,
                crop=poly,
                thresh_in=thresh_in,
                remove_isolated=remove_isolated,
            )

    def add_mesh_2d_tri(
        self,
        # rectangular mesh
        key=None,
        # triangular mesh
        knots=None,
        indices=None,
        # from points and polygon
        pts_x0=None,
        pts_x1=None,
        # cropping
        crop_poly=None,
        thresh_in=None,
        remove_isolated=None,
        # direct addition of bsplines
        deg=None,
        **kwdargs,
    ):
        """ Add a mesh by key and domain / resolution

        Can create a 2d tringular mesh:
            - knots: (nknots, 2) array of (R, Z) coordinates
            - cents: (ncents, 3 or 4) array of int indices

        If deg is provided, immediately adds a bsplines

        """

        # get domain, poly from crop_poly
        if crop_poly is not None:
            poly = _checks_2d_rect._from_croppoly(
                crop_poly=crop_poly,
                domain=None,
            )[1]
        else:
            poly = None

        # check input data and get input dicts
        key, dref, ddata, dobj = _checks_2d_tri.check(
            coll=self,
            # key
            key=key,
            # triangular
            knots=knots,
            indices=indices,
            # from pts and polygon
            pts_x0=pts_x0,
            pts_x1=pts_x1,
            # crop poly
            poly=poly,
            # others
            **kwdargs,
        )

        # update dicts
        self.update(dref=dref, ddata=ddata, dobj=dobj)

        # optional bspline
        if deg is not None:
            self.add_bsplines(key=key, deg=deg)

    # def add_mesh_2d_polar(
        # self,
        # # polar mesh
        # key=None,
        # radius=None,
        # angle=None,
        # # Defined on
        # radius2d=None,
        # angle2d=None,
        # # res for contour discontinuity of angle2d
        # res=None,
        # # parameters
        # radius_dim=None,
        # radius_quant=None,
        # radius_name=None,
        # radius_units=None,
        # angle_dim=None,
        # angle_quant=None,
        # angle_name=None,
        # # direct addition of bsplines
        # deg=None,
        # **kwdargs,
    # ):
        # """ Add a 2d polar mesh

        # For now only includes radial mesh
        # radius has to be backed-up by:
            # - a radius quantity from a pre-existing rect or tri mesh
            # - a function

        # """

        # # check input data and get input dicts
        # dref, ddata, dmesh = _checks_2d_polar.check(
            # coll=self,
            # # key
            # key=key,
            # # polar
            # radius=radius,
            # angle=angle,
            # radius2d=radius2d,
            # angle2d=angle2d,
            # # parameters
            # radius_dim=radius_dim,
            # radius_quant=radius_quant,
            # radius_name=radius_name,
            # radius_units=radius_units,
            # angle_dim=angle_dim,
            # angle_quant=angle_quant,
            # angle_name=angle_name,
        # )

        # # update dicts
        # self.update(dref=dref, ddata=ddata, dobj=dobj)

        # # special treatment of angle2d
        # if dmesh[key]['angle2d'] is not None:
            # drefa, ddataa, kR, kZ = _compute.angle2d_zone(
                # coll=self,
                # key=dmesh[key]['angle2d'],
                # keyrad2d=dmesh[key]['radius2d'],
                # key_ptsO=dmesh[key]['pts_O'],
                # res=res,
                # keym0=key,
            # )

            # # update dicts
            # self.update(dref=drefa, ddata=ddataa)
            # if 'azone' in self.get_lparam(self._which_mesh):
                # self.set_param(
                    # key=key,
                    # param='azone',
                    # value=(kR, kZ),
                    # which=self._which_mesh,
                # )
            # else:
                # self.add_param(
                    # 'azone', value={key: (kR, kZ)}, which=self._which_mesh,
                # )

        # # optional bspline
        # if deg is not None:
            # self.add_bsplines(key=key, deg=deg)

    # -----------------
    # remove mesh
    # ------------------

    def remove_mesh(self, key=None, propagate=None):
        return _generic_mesh.remove_mesh(
            coll=self,
            key=key,
            propagate=propagate,
        )

    # -----------------
    # crop
    # ------------------

    def crop(self, key=None, crop=None, thresh_in=None, remove_isolated=None):
        """ Crop a mesh / bspline

        Uses:
            - a mask of bool for each mesh elements
            - a 2d (R, Z) closed polygon

        If applied on a bspline, cropping is double-checked to make sure
        all remaining bsplines have full support domain
        """

        crop, key, thresh_in = _cropping.crop(
            coll=self,
            key=key,
            crop=crop,
            thresh_in=thresh_in,
            remove_isolated=remove_isolated,
        )

        # add crop data
        wm = self._which_mesh
        mtype = self.dobj[wm][key]['type']

        keycrop = f'{key}_crop'
        if mtype == 'rect':
            ref = tuple([
                self._ddata[k0]['ref'][0]
                for k0 in self._dobj[wm][key]['cents']
            ])
        else:
            ref = self.ddata[self.dobj[wm][key]['cents'][0]]['ref']

        self.add_data(
            key=keycrop,
            data=crop,
            ref=ref,
            dim='bool',
            quant='bool',
        )

        # update obj
        self._dobj[wm][key]['crop'] = keycrop
        self._dobj[wm][key]['crop_thresh'] = thresh_in

        # safety check
        self.get_mesh_outline(key, plot_debug=True)

    # -----------------
    # indices
    # ------------------

    def select_ind(
        self,
        key=None,
        ind=None,
        elements=None,
        returnas=None,
        crop=None,
    ):
        """ Return ind for selected key (mesh or bspline) as:
                - tuple (default)
                - 'flat'

        Can covert one into the other
        """
        return _select._select_ind(
            coll=self,
            key=key,
            ind=ind,
            elements=elements,
            returnas=returnas,
            crop=crop,
        )

    def select_mesh_elements(
        self,
        key=None,
        ind=None,
        elements=None,
        returnas=None,
        return_neighbours=None,
        crop=None,
    ):
        """ Return indices or values of selected knots / cent

        Can be used to convert tuple (R, Z) indices to flat (RZ,) indices
        Can return values instead of indices
        Can return indices / values of neighbourgs

        """
        # check key
        key, _, _ = _generic_mesh._get_key_mesh_vs_bplines(
            coll=self,
            key=key,
            which=self._which_mesh,
        )

        # get ind
        mtype = self.dobj[self._which_mesh][key]['type']
        if mtype == 'rect':
            returnas_ind = tuple
        else:
            returnas_ind = bool

        ind = self.select_ind(
            key=key,
            ind=ind,
            elements=elements,
            returnas=returnas_ind,
            crop=crop,
        )

        return _select._select_mesh(
            coll=self,
            key=key,
            ind=ind,
            elements=elements,
            returnas=returnas,
            return_neighbours=return_neighbours,
        )

    # -----------------
    # interp tools
    # ------------------

    def get_sample_mesh(
        self,
        key=None,
        res=None,
        grid=None,
        mode=None,
        x0=None,
        x1=None,
        Dx0=None,
        Dx1=None,
        imshow=None,
        in_mesh=None,
        # store
        store=None,
        kx0=None,
        kx1=None,
    ):
        """ Return a sampled version of the chosen mesh """
        return _sample.sample_mesh(
            coll=self,
            key=key,
            res=res,
            grid=grid,
            mode=mode,
            x0=x0,
            x1=x1,
            Dx0=Dx0,
            Dx1=Dx1,
            imshow=imshow,
            in_mesh=in_mesh,
            # store
            store=store,
            kx0=kx0,
            kx1=kx1,
        )

    # -----------------
    # outline
    # ------------------

    def get_mesh_outline(
        self,
        key=None,
        closed=None,
        plot_debug=None,
    ):
        """ Get outline of 2d mesh """
        return _outline._get_outline(
            coll=self,
            key=key,
            closed=closed,
            plot_debug=plot_debug,
        )

    # -----------------
    # plotting
    # ------------------

    def plot_mesh(
        self,
        key=None,
        ind_knot=None,
        ind_cent=None,
        crop=None,
        color=None,
        dax=None,
        dmargin=None,
        fs=None,
        dleg=None,
        connect=None,
    ):

        return _plot.plot_mesh(
            coll=self,
            key=key,
            ind_knot=ind_knot,
            ind_cent=ind_cent,
            crop=crop,
            color=color,
            dax=dax,
            dmargin=dmargin,
            fs=fs,
            dleg=dleg,
            connect=connect,
        )