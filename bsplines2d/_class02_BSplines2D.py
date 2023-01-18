# -*- coding: utf-8 -*-


# Built-in
import copy


# Common
import numpy as np


# local
from ._class01_Mesh2D import Mesh2D as Previous
from . import _class02_checks as _checks
from . import _class02_compute as _compute
from . import _class01_rect_cropping as _cropping
from . import _class02_operators as _operators


__all__ = ['BSplines2D']


# #############################################################################
# #############################################################################
#                       
# #############################################################################


class BSplines2D(Previous):

    _ddef = copy.deepcopy(Previous._ddef)
    _ddef['params']['ddata'].update({
        'bsplines': {'cls': (tuple, str), 'def': ''},
    })
    _ddef['params']['dobj'] = None
    _ddef['params']['dref'] = None

    _dshow = dict(Previous._dshow)
 
    # -----------------
    # bsplines
    # ------------------

    def add_bsplines(self, key=None, deg=None):
        """ Add bspline basis functions on the chosen mesh """

        # --------------
        # check inputs

        keym, keybs, deg = _checks._mesh_bsplines(
            key=key,
            lkeys=list(self.dobj[self._which_mesh].keys()),
            deg=deg,
        )

        # ------------
        # get bsplines

        nd = self.dobj[self._which_mesh][keym]['nd']
        mtype = self.dobj[self._which_mesh][keym]['type']
        if nd == '1d':
            dref, ddata, dobj = _compute._mesh1d_bsplines(
                coll=self, keym=keym, keybs=keybs, deg=deg,
            )
        elif mtype == 'rect':
            dref, ddata, dobj = _compute._mesh2DRect_bsplines(
                coll=self, keym=keym, keybs=keybs, deg=deg,
            )
        elif mtype == 'tri':
            dref, ddata, dobj = _compute._mesh2DTri_bsplines(
                coll=self, keym=keym, keybs=keybs, deg=deg,
            )
        else:
            dref, ddata, dobj = _compute._mesh2Dpolar_bsplines(
                coll=self, keym=keym, keybs=keybs, deg=deg, # angle=angle,
            )

        # --------------
        # update dict and crop if relevant

        self.update(dobj=dobj, ddata=ddata, dref=dref)
        if mtype == 'rect':
            _cropping.add_cropbs_from_crop(
                coll=self,
                keybs=keybs,
                keym=keym,
            )

    # -----------------
    # add_data
    # ------------------

    def update(
        self,
        dobj=None,
        ddata=None,
        dref=None,
        harmonize=None,
    ):
        " Overload datastock update() method "

        # if ddata => check ref for bsplines
        if ddata is not None:
            for k0, v0 in ddata.items():
                (
                    ddata[k0]['ref'], ddata[k0]['data'],
                ) = _checks.add_data_meshbsplines_ref(
                    coll=self,
                    ref=v0['ref'],
                    data=v0['data'],
                )

        # update
        super().update(
            dobj=dobj,
            ddata=ddata,
            dref=dref,
            harmonize=harmonize,
        )

        # assign bsplines
        _checks._set_data_bsplines(coll=self)
        
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
        super().crop(
            key=key, 
            crop=crop, 
            thresh_in=thresh_in,
            remove_isolated=remove_isolated,
        )

        # also crop bsplines
        for k0 in self.dobj.get(self._which_bsplines, {}).keys():
            if self.dobj[self._which_bsplines][k0][self._which_mesh] == key:
                _cropping.add_cropbs_from_crop(coll=self, keybs=k0, keym=key)

    # -----------------
    # get data subset
    # ------------------

    def get_dict_bsplines(self):
        """ Return dict of profiles2d with associated bsplines as values """

        # dict of profiles2d
        dk = {
            k0: v0[self._which_bsplines]
            for k0, v0 in self._ddata.items()
            if v0[self._which_bsplines] not in [None, '']
        }
        dk.update({k0: k0 for k0 in self._dobj[self._which_bsplines].keys()})

        return dk

    # -----------------
    # indices
    # ------------------

    def select_bsplines(
        self,
        key=None,
        ind=None,
        returnas=None,
        return_cents=None,
        return_knots=None,
        crop=None,
    ):
        """ Return indices or values of selected knots / cent

        Can be used to convert tuple (R, Z) indices to flat (RZ,) indices
        Can return values instead of indices
        Can return indices / values of neighbourgs

        """
        return _compute._select_bsplines(
            coll=self,
            key=key,
            ind=ind,
            returnas=returnas,
            return_cents=return_cents,
            return_knots=return_knots,
            crop=crop,
        )

    # -----------------
    # Integration operators
    # ------------------

    def add_bsplines_operator(
        self,
        key=None,
        operator=None,
        geometry=None,
        crop=None,
        store=None,
        returnas=None,
        # specific to deg = 0
        centered=None,
        # to return gradR, gradZ, for D1N2 deg 0, for tomotok
        returnas_element=None,
    ):
        """ Get a matrix operator to compute an integral

        operator specifies the integrand:
            - 'D0N1': integral of the value
            - 'D0N2': integral of the squared value
            - 'D1N2': integral of the squared gradient
            - 'D2N2': integral of the squared laplacian

        geometry specifies in which geometry:
            - 'linear': linear geometry (cross-section = surface)
            - 'toroidal': toroildal geometry (cross-section = volumic slice)

        """

        return _operators.get_bsplines_operator(
            coll=self,
            key=key,
            operator=operator,
            geometry=geometry,
            crop=crop,
            store=store,
            returnas=returnas,
            # specific to deg = 0
            centered=centered,
            # to return gradR, gradZ, for D1N2 deg 0, for tomotok
            returnas_element=returnas_element,
        )

    # -----------------
    # binning tools
    # ------------------

    def binning(
        self,
        keys=None,
        bins=None,
        axis=None,
        key_bs=None,
    ):
        """ Return binned data along the desired axis """
        
        return _binning.binning(
            coll=self,
            keys=keys,
            bins=bins,
            axis=axis,
            key_bs=key_bs,
        )
        

    # -----------------
    # interp tools
    # ------------------

    def _interpolate_bsplines(
        self,
        # interpolation base, 1d or 2d
        key=None,
        # external coefs (optional)
        coefs=None,
        # interpolation points
        x0=None,
        x1=None,
        grid=None,
        # time: t or indt
        t=None,
        indt=None,
        indt_strict=None,
        # bsplines
        indbs=None,
        # parameters
        details=None,
        reshape=None,
        res=None,
        crop=None,
        nan0=None,
        val_out=None,
        imshow=None,
        return_params=None,
        # storing
        store=None,
        inplace=None,
    ):
        
        
        


    """
    def get_sample_bspline(self, key=None, res=None, grid=None, mode=None):
        return _compute.sample_bsplines(
            coll=self,
            key=key,
            res=res,
            grid=grid,
            mode=mode,
        )
    """

    def _check_qr12RPZ(
        self,
        quant=None,
        ref1d=None,
        ref2d=None,
        q2dR=None,
        q2dPhi=None,
        q2dZ=None,
        group1d=None,
        group2d=None,
    ):

        if group1d is None:
            group1d = self._group1d
        if group2d is None:
            group2d = self._group2d

        lc0 = [quant is None, ref1d is None, ref2d is None]
        lc1 = [q2dR is None, q2dPhi is None, q2dZ is None]
        if np.sum([all(lc0), all(lc1)]) != 1:
            msg = (
                "Please provide either (xor):\n"
                + "\t- a scalar field (isotropic emissivity):\n"
                + "\t\tquant : scalar quantity to interpolate\n"
                + "\t\t\tif quant is 1d, intermediate reference\n"
                + "\t\t\tfields are necessary for 2d interpolation\n"
                + "\t\tref1d : 1d reference field on which to interpolate\n"
                + "\t\tref2d : 2d reference field on which to interpolate\n"
                + "\t- a vector (R,Phi,Z) field (anisotropic emissivity):\n"
                + "\t\tq2dR :  R component of the vector field\n"
                + "\t\tq2dPhi: R component of the vector field\n"
                + "\t\tq2dZ :  Z component of the vector field\n"
                + "\t\t=> all components have the same time and mesh!\n"
            )
            raise Exception(msg)

        # Check requested quant is available in 2d or 1d
        if all(lc1):
            (
                idquant, idref1d, idref2d,
            ) = _compute._get_possible_ref12d(
                dd=self._ddata,
                key=quant, ref1d=ref1d, ref2d=ref2d,
                group1d=group1d,
                group2d=group2d,
            )
            idq2dR, idq2dPhi, idq2dZ = None, None, None
            ani = False
        else:
            idq2dR, msg = _compute._get_keyingroup_ddata(
                dd=self._ddata,
                key=q2dR, group=group2d, msgstr='quant', raise_=True,
            )
            idq2dPhi, msg = _compute._get_keyingroup_ddata(
                dd=self._ddata,
                key=q2dPhi, group=group2d, msgstr='quant', raise_=True,
            )
            idq2dZ, msg = _compute._get_keyingroup_ddata(
                dd=self._ddata,
                key=q2dZ, group=group2d, msgstr='quant', raise_=True,
            )
            idquant, idref1d, idref2d = None, None, None
            ani = True
        return idquant, idref1d, idref2d, idq2dR, idq2dPhi, idq2dZ, ani

    def _interp_pts2d_to_quant1d(
        self,
        pts=None,
        vect=None,
        t=None,
        quant=None,
        ref1d=None,
        ref2d=None,
        q2dR=None,
        q2dPhi=None,
        q2dZ=None,
        interp_t=None,
        interp_space=None,
        fill_value=None,
        Type=None,
        group0d=None,
        group1d=None,
        group2d=None,
        return_all=None,
    ):
        """ Return the value of the desired 1d quantity at 2d points

        For the desired inputs points (pts):
            - pts are in (X, Y, Z) coordinates
            - space interpolation is linear on the 1d profiles
        At the desired input times (t):
            - using a nearest-neighbourg approach for time

        """
        # Check inputs
        if group0d is None:
            group0d = self._group0d
        if group1d is None:
            group1d = self._group1d
        if group2d is None:
            group2d = self._group2d
        # msg = "Only 'nearest' available so far for interp_t!"
        # assert interp_t == 'nearest', msg

        # Check requested quant is available in 2d or 1d
        idquant, idref1d, idref2d, idq2dR, idq2dPhi, idq2dZ, ani = \
                self._check_qr12RPZ(
                    quant=quant, ref1d=ref1d, ref2d=ref2d,
                    q2dR=q2dR, q2dPhi=q2dPhi, q2dZ=q2dZ,
                    group1d=group1d, group2d=group2d,
                )

        # Check the pts is (3,...) array of floats
        idmesh = None
        if pts is None:
            # Identify mesh to get default points
            if ani:
                idmesh = [id_ for id_ in self._ddata[idq2dR]['ref']
                          if self._dref[id_]['group'] == group2d][0]
            else:
                if idref1d is None:
                    idmesh = [id_ for id_ in self._ddata[idquant]['ref']
                              if self._dref[id_]['group'] == group2d][0]
                else:
                    idmesh = [id_ for id_ in self._ddata[idref2d]['ref']
                              if self._dref[id_]['group'] == group2d][0]

            # Derive pts
            pts = self._get_pts_from_mesh(key=idmesh)

        pts = np.atleast_2d(pts)
        if pts.shape[0] != 3:
            msg = (
                "pts must be np.ndarray of (X,Y,Z) points coordinates\n"
                + "Can be multi-dimensional, but 1st dimension is (X,Y,Z)\n"
                + "    - Expected shape : (3,...)\n"
                + "    - Provided shape : {}".format(pts.shape)
            )
            raise Exception(msg)

        # Check t
        lc = [t is None, type(t) is str, type(t) is np.ndarray]
        assert any(lc)
        if lc[1]:
            assert t in self._ddata.keys()
            t = self._ddata[t]['data']

        # Interpolation (including time broadcasting)
        # this is the second slowest step (~0.08 s)
        func = self._get_finterp(
            idquant=idquant, idref1d=idref1d, idref2d=idref2d,
            idq2dR=idq2dR, idq2dPhi=idq2dPhi, idq2dZ=idq2dZ,
            idmesh=idmesh,
            interp_t=interp_t, interp_space=interp_space,
            fill_value=fill_value, ani=ani, Type=Type,
            group0d=group0d, group2d=group2d,
        )

        # Check vect of ani
        c0 = (
            ani is True
            and (
                vect is None
                or not (
                    isinstance(vect, np.ndarray)
                    and vect.shape == pts.shape
                )
            )
        )
        if c0:
            msg = (
                "Anisotropic field interpolation needs a field of local vect\n"
                + "  => Please provide vect as (3, npts) np.ndarray!"
            )
            raise Exception(msg)

        # This is the slowest step (~1.8 s)
        val, t = func(pts, vect=vect, t=t)

        # return
        if return_all is None:
            return_all = True
        if return_all is True:
            dout = {
                't': t,
                'pts': pts,
                'ref1d': idref1d,
                'ref2d': idref2d,
                'q2dR': idq2dR,
                'q2dPhi': idq2dPhi,
                'q2dZ': idq2dZ,
                'interp_t': interp_t,
                'interp_space': interp_space,
            }
            return val, dout
        else:
            return val

    def _interp_pts2d_to_quant1d(
        self,
        pts=None,
        vect=None,
        t=None,
        quant=None,
        ref1d=None,
        ref2d=None,
        q2dR=None,
        q2dPhi=None,
        q2dZ=None,
        interp_t=None,
        interp_space=None,
        fill_value=None,
        Type=None,
        group0d=None,
        group1d=None,
        group2d=None,
        return_all=None,
    ):
        """ Return the value of the desired 1d quantity at 2d points

        For the desired inputs points (pts):
            - pts are in (X, Y, Z) coordinates
            - space interpolation is linear on the 1d profiles
        At the desired input times (t):
            - using a nearest-neighbourg approach for time

        """
        # Check inputs
        if group0d is None:
            group0d = self._group0d
        if group1d is None:
            group1d = self._group1d
        if group2d is None:
            group2d = self._group2d
        # msg = "Only 'nearest' available so far for interp_t!"
        # assert interp_t == 'nearest', msg

        # Check requested quant is available in 2d or 1d
        idquant, idref1d, idref2d, idq2dR, idq2dPhi, idq2dZ, ani = \
                self._check_qr12RPZ(
                    quant=quant, ref1d=ref1d, ref2d=ref2d,
                    q2dR=q2dR, q2dPhi=q2dPhi, q2dZ=q2dZ,
                    group1d=group1d, group2d=group2d,
                )

        # Check the pts is (3,...) array of floats
        idmesh = None
        if pts is None:
            # Identify mesh to get default points
            if ani:
                idmesh = [id_ for id_ in self._ddata[idq2dR]['ref']
                          if self._dref[id_]['group'] == group2d][0]
            else:
                if idref1d is None:
                    idmesh = [id_ for id_ in self._ddata[idquant]['ref']
                              if self._dref[id_]['group'] == group2d][0]
                else:
                    idmesh = [id_ for id_ in self._ddata[idref2d]['ref']
                              if self._dref[id_]['group'] == group2d][0]

            # Derive pts
            pts = self._get_pts_from_mesh(key=idmesh)

        pts = np.atleast_2d(pts)
        if pts.shape[0] != 3:
            msg = (
                "pts must be np.ndarray of (X,Y,Z) points coordinates\n"
                + "Can be multi-dimensional, but 1st dimension is (X,Y,Z)\n"
                + "    - Expected shape : (3,...)\n"
                + "    - Provided shape : {}".format(pts.shape)
            )
            raise Exception(msg)

        # Check t
        lc = [t is None, type(t) is str, type(t) is np.ndarray]
        assert any(lc)
        if lc[1]:
            assert t in self._ddata.keys()
            t = self._ddata[t]['data']

        # Interpolation (including time broadcasting)
        # this is the second slowest step (~0.08 s)
        func = self._get_finterp(
            idquant=idquant, idref1d=idref1d, idref2d=idref2d,
            idq2dR=idq2dR, idq2dPhi=idq2dPhi, idq2dZ=idq2dZ,
            idmesh=idmesh,
            interp_t=interp_t, interp_space=interp_space,
            fill_value=fill_value, ani=ani, Type=Type,
            group0d=group0d, group2d=group2d,
        )

        # Check vect of ani
        c0 = (
            ani is True
            and (
                vect is None
                or not (
                    isinstance(vect, np.ndarray)
                    and vect.shape == pts.shape
                )
            )
        )
        if c0:
            msg = (
                "Anisotropic field interpolation needs a field of local vect\n"
                + "  => Please provide vect as (3, npts) np.ndarray!"
            )
            raise Exception(msg)

        # This is the slowest step (~1.8 s)
        val, t = func(pts, vect=vect, t=t)

        # return
        if return_all is None:
            return_all = True
        if return_all is True:
            dout = {
                't': t,
                'pts': pts,
                'ref1d': idref1d,
                'ref2d': idref2d,
                'q2dR': idq2dR,
                'q2dPhi': idq2dPhi,
                'q2dZ': idq2dZ,
                'interp_t': interp_t,
                'interp_space': interp_space,
            }
            return val, dout
        else:
            return val

    def interpolate_profile2d(
        # ressources
        self,
        # interpolation base, 1d or 2d
        key=None,
        # external coefs (optional)
        coefs=None,
        # interpolation points
        R=None,
        Z=None,
        radius=None,
        angle=None,
        grid=None,
        radius_vs_time=None,
        azone=None,
        # time: t or indt
        t=None,
        indt=None,
        indt_strict=None,
        # bsplines
        indbs=None,
        # parameters
        details=None,
        reshape=None,
        res=None,
        crop=None,
        nan0=None,
        val_out=None,
        imshow=None,
        return_params=None,
        # storing
        store=None,
        inplace=None,
    ):
        """ Interpolate desired profile2d (i.e.: data on bsplines)

        Interpolate:
            - key: a data on bsplines
            - coefs: external-provided set of coefs

        coefs can only be provided if:
            - details = False
            - key = keybs
            - coefs is a scalar or has shape = shapebs

        At points:
            - R:  R coordinates (np.ndarray or scalar)
            - Z:  Z coordinates (np.ndarray, same shape as R, or scalar)
            - grid: bool, if True mesh R x Z
            - indt: if provided, only interpolate at desired time indices

        With options:
            - details: bool, if True returns value for each bspline
            - indbs:   optional, select bsplines for which to interpolate
            - reshape: bool,
            - res:  optional, resolution to generate R and Z if they are None
            - crop: bool, whether to use the cropped mesh
            - nan0: value for out-of-mesh points
            - imshow: bool, whether to return as imshow (transpose)
            - return_params: bool, whether to return dict of input params

        """

        return _compute.interp2d(
            # ressources
            coll=self,
            # interpolation base, 1d or 2d
            key=key,
            # external coefs (optional)
            coefs=coefs,
            # interpolation points
            R=R,
            Z=Z,
            radius=radius,
            angle=angle,
            grid=grid,
            radius_vs_time=radius_vs_time,
            azone=azone,
            # time: t or indt
            t=t,
            indt=indt,
            indt_strict=indt_strict,
            # bsplines
            indbs=indbs,
            # parameters
            details=details,
            reshape=reshape,
            res=res,
            crop=crop,
            nan0=nan0,
            val_out=val_out,
            imshow=imshow,
            return_params=return_params,
            # storing
            store=store,
            inplace=inplace,
        )

    # TBF after polar meshes
    def interpolate_2dto1d(
        # resources
        self,
        # interpolation base
        key1d=None,
        key2d=None,
        # interpolation pts
        R=None,
        Z=None,
        grid=None,
        # parameters
        interp_t=None,
        fill_value=None,
        ani=False,
    ):

        return _compute.interp2dto1d(
            coll=self,
            key1d=key1d,
            key2d=key2d,
            R=R,
            Z=Z,
            grid=grid,
            crop=crop,
            nan0=nan0,
            return_params=return_params,
        )

    # -----------------
    # plotting
    # ------------------

    def plot_bsplines(
        self,
        key=None,
        indbs=None,
        indt=None,
        knots=None,
        cents=None,
        res=None,
        plot_mesh=None,
        val_out=None,
        nan0=None,
        cmap=None,
        dax=None,
        dmargin=None,
        fs=None,
        dleg=None,
    ):

        return _plot.plot_bspline(
            coll=self,
            key=key,
            indbs=indbs,
            indt=indt,
            knots=knots,
            cents=cents,
            res=res,
            plot_mesh=plot_mesh,
            val_out=val_out,
            nan0=nan0,
            cmap=cmap,
            dax=dax,
            dmargin=dmargin,
            fs=fs,
            dleg=dleg,
        )

    def plot_profile2d(
        self,
        # inputs
        key=None,
        coefs=None,
        indt=None,
        res=None,
        # plot options
        vmin=None,
        vmax=None,
        cmap=None,
        dax=None,
        dmargin=None,
        fs=None,
        dcolorbar=None,
        dleg=None,
        # interactivity
        dinc=None,
        connect=None,
    ):
        return _plot.plot_profile2d(
            coll=self,
            # inputs
            key=key,
            coefs=coefs,
            indt=indt,
            res=res,
            # plot options
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            dax=dax,
            dmargin=dmargin,
            fs=fs,
            dcolorbar=dcolorbar,
            dleg=dleg,
            # interactivity
            dinc=dinc,
            connect=connect,
        )
