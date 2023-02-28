# -*- coding: utf-8 -*-


# Built-in
import warnings


# Common
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import datastock as ds

# specific
from . import _generic_check
from . import _class01_checks as _checks


# #############################################################################
# #############################################################################
#                           Main
# #############################################################################


def plot_as_profile2d(
    # ressources
    coll=None,
    # inputs
    key=None,
    # parameters
    dres=None,
    # contours
    dlevels=None,
    ref_com=None,
    # figure
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

    # --------------
    # check input

    if connect is None:
        connect = True

    (
        key, keybs, keym, nd, mtype,
        submesh, subbs, submtype,
        dlevels,
        cmap, dcolorbar, dleg,
        connect,
    ) = _check(
        coll=coll,
        key=key,
        dlevels=dlevels,
        # figure
        cmap=cmap,
        dcolorbar=dcolorbar,
        dleg=dleg,
        # interactivity
        connect=connect,
    )

    # --------------
    #  Prepare data

    (
        coll2, dkeys, interp,
        dlevels, lcol,
    ) = _prepare(
        coll=coll,
        key=key,
        keybs=keybs,
        keym=keym,
        dres=dres,
        mtype=mtype,
        # submesh
        submesh=submesh,
        subbs=subbs,
        submtype=submtype,
        # levels
        dlevels=dlevels,
        ref_com=ref_com,
    )

    # -----------------
    # case with submesh

    if submesh is not None:

        dax, dgroup = _plot_submesh(
            coll=coll,
            coll2=coll2,
            key=key,
            keym=keym,
            keybs=keybs,
            # figure
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            dax=dax,
            dmargin=dmargin,
            fs=fs,
            dcolorbar=dcolorbar,
            dleg=dleg,
            interp=interp,
            dkeys=dkeys,
            lcol=lcol,
        )

    # -------------------
    # call right function

    else:
        dax, dgroup = coll2.plot_as_array(
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            dax=dax,
            dmargin=dmargin,
            fs=fs,
            dcolorbar=dcolorbar,
            dleg=dleg,
            interp=interp,
            connect=False,
            **dkeys,
        )

        
    # -----------
    # add levels
    
    if dlevels is not None:
        _add_levels_2d(
            key=key,
            dax=dax,
            dgroup=dgroup,
            coll2=coll2,
            dlevels=dlevels,
            # figure
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            dmargin=dmargin,
            fs=fs,
            dcolorbar=dcolorbar,
            dleg=dleg,
            interp=interp,
            dkeys=dkeys,
        )
    
    
    # -----------
    # connect

    if connect is True:
        dax.setup_interactivity(kinter='inter0', dgroup=dgroup, dinc=dinc)
        dax.disconnect_old()
        dax.connect()

        dax.show_commands()
        return dax
    else:
        return dax, dgroup
    
    
# #############################################################################
# #############################################################################
#                       Check
# #############################################################################


def _check(
    coll=None,
    key=None,
    dlevels=None,
    # figure
    cmap=None,
    dcolorbar=None,
    dleg=None,
    # interactivity
    connect=None,
):

    # ----------
    # keys

    # key
    dk = coll.get_profiles2d()
    key = ds._generic_check._check_var(
        key, 'key',
        types=str,
        allowed=list(dk.keys()),
    )

    wm = coll._which_mesh
    wbs = coll._which_bsplines

    keybs = dk[key]
    # refbs = coll.dobj[wbs][keybs]['ref']

    keym = coll.dobj[wbs][keybs][wm]
    nd = coll.dobj[wm][keym]['nd']
    mtype = coll.dobj[wm][keym]['type']

    submesh = coll.dobj[wm][keym]['submesh']
    if submesh == '':
        submesh = None

    if submesh is None:
        subbs = keybs
        submtype = mtype
    else:
        subbs = coll.dobj[wm][keym]['subbs']
        # subkey = coll.dobj[wm][keym]['subkey']
        submtype = coll.dobj[wm][keym]['type']

    # ----------
    # dlevels
    
    if dlevels is not None:
        
        dp2d = coll.get_profiles2d()
        if isinstance(dlevels, (float, int)):
            dlevels = {key: np.r_[dlevels]}
            
        elif isinstance(dlevels, (np.ndarray, list, tuple)):
            dlevels = np.atleast_1d(dlevels).ravel()
            dlevels = {key: dlevels}
            
        c0 = (
            isinstance(dlevels, dict)
            and all([kk in dp2d.keys() for kk, vv in dlevels.items()])
        )
        if not c0:
            msg = (
                "Arg dlevels must be a dict with:\n"
                "\t- keys: valid keys of 2d profile data\n"
                "\t- values: iterable of level values\n"
                f"\nProvided: {dlevels}"
            )
            raise Exception(msg)
            
        for k0, v0 in dlevels.items():
            
            if isinstance(v0, (np.ndarray, list, tuple)):
                dlevels[k0] = {'levels': np.atleast_1d(v0).ravel()}
                v0 = dlevels[k0]
            
            if not isinstance(v0, dict) or 'levels' not in v0.keys():
                msg = f"dlevels['{k0}'] must have key 'levels'"
                raise Exception(msg)
            
            # check fields
            dlevels[k0]['levels'] = np.atleast_1d(v0['levels']).ravel()
            dlevels[k0]['color'] = v0.get('color', 'k')

    # ----------
    # figure

    # cmap
    if cmap is None:
        cmap = 'viridis'

    # dcolorbar
    defdcolorbar = {
        # 'location': 'right',
        'fraction': 0.15,
        'orientation': 'vertical',
    }
    dcolorbar = ds._generic_check._check_var(
        dcolorbar, 'dcolorbar',
        default=defdcolorbar,
        types=dict,
    )

    # dleg
    defdleg = {
        'bbox_to_anchor': (1.1, 1.),
        'loc': 'upper left',
        'frameon': True,
    }
    dleg = ds._generic_check._check_var(
        dleg, 'dleg',
        default=defdleg,
        types=(bool, dict),
    )

    # ----------
    # interactivity

    connect = ds._generic_check._check_var(
        connect, 'connect',
        types=bool,
        default=True,
    )

    return (
        key, keybs, keym, nd, mtype,
        submesh, subbs, submtype,
        dlevels,
        cmap, dcolorbar, dleg,
        connect,
    )


# #############################################################################
# #############################################################################
#                       Prepare
# #############################################################################


def _prepare(
    coll=None,
    key=None,
    keybs=None,
    keym=None,
    coefs=None,
    indt=None,
    dres=None,
    mtype=None,
    # submesh
    submesh=None,
    subbs=None,
    submtype=None,
    # levels
    dlevels=None,
    ref_com=None,
):

    # ------------
    # misc

    # deg and
    # wm = coll._which_mesh
    wbs = coll._which_bsplines
    deg = coll.dobj[wbs][subbs]['deg']
    if deg == 0:
        interp = 'nearest'
    elif deg == 1:
        interp = 'bilinear'
    elif deg >= 2:
        interp = 'bicubic'

    # -----------------
    # get plotting mesh

    coll2, dbs = coll.interpolate_all_bsplines(
        key=key,
        dres=dres,
        submesh=True,
    )

    bs2d = [k0 for k0, v0 in dbs.items() if len(v0['ref']) == 2][0]
    rX, rY = dbs[bs2d]['ref']
    lr1d = [k0 for k0 in coll2.ddata[key]['ref'] if k0 not in [rX, rY]]
    ndim = coll2.ddata[key]['data'].ndim

    dkeys = {
        'key': key,
        'keyX': coll2.get_ref_vector(ref=rX)[3],
        'keyY': coll2.get_ref_vector(ref=rY)[3],
        'keyZ': None,
        'keyU': None,
    }

    if ndim >= 3:
        dkeys['keyZ'] = coll2.get_ref_vector(ref=lr1d[0])[3]
        # uniform = ds._plot_as_array._check_uniform_lin(
            # k0=keyZ, ddata=coll2.ddata,
        # )
        # if not uniform:
            # keyZ = None
        if ndim == 4:
            dkeys['keyU'] = coll2.get_ref_vector(ref=lr1d[1])[3]

    # -----------------
    # optional contours

    if dlevels is not None:

        for ii, (k0, v0) in enumerate(dlevels.items()):        

            refi = coll.ddata[k0]['ref']            

            # get contours
            dout, dref = coll.get_profile2d_contours(
                key=k0,
                levels=v0['levels'],
                ref_com=ref_com if ref_com in refi else None,
                res=dres if isinstance(dres, (int, float)) else None,
                store=False,
                return_dref=True,
                key_cont0="cont0",
                key_cont1="cont1",
            )
    
            # ref
            dlevels[k0]['dref'] = dref
            
            # axis
            ref = dout['cont0']['ref']
            axis = [
                ii for ii, rr in enumerate(ref)
                if rr not in coll.dref.keys()
            ]
    
            # refZ, refU
            refZ, refU = None, None
            if len(ref) == 3:
                refZ = [
                    rr for ii, rr in enumerate(ref)
                    if ii not in axis
                ][0]
                
            # populate
            for k1 in ['cont0', 'cont1']:
                dlevels[k0][k1] = dout[k1]
                dlevels[k0]['cont0']['key'] = f'{k0}_{k1}'
                
            dlevels[k0]['refZ'] = refZ
            dlevels[k0]['refU'] = refU
            dlevels[k0]['axis'] = axis

    # --------------------
    # optional for submesh

    lcol = ['k', 'r', 'b', 'g', 'm', 'c', 'y']

    return (
        coll2, dkeys, interp,
        dlevels, lcol,
    )


# #############################################################################
# #############################################################################
#                   Plot fixed mesh with levels
# #############################################################################


def _add_levels_2d(
    key=None,
    dax=None,
    dgroup=None,
    coll2=None,
    dlevels=None,
    # figure
    vmin=None,
    vmax=None,
    cmap=None,
    dmargin=None,
    fs=None,
    dcolorbar=None,
    dleg=None,
    interp=None,
    dkeys=None,
):

    # ---------------------------------
    # add make contours as single lines

    ndim = len(dgroup)
    for ii, (k0, v0) in enumerate(dlevels.items()):
        
        for k1 in ['cont0', 'cont1']:
        
            v1 = v0[k1]
            sh = v1['data'].shape
            shnan = [1 if ii == v0['axis'][0] else ss for ii, ss in enumerate(sh)]
    
            dlevels[k0][k1]['data'] = np.append(
                v1['data'],
                np.full(tuple(shnan), np.nan),
                axis=v0['axis'][0],
            )
    
            sh = dlevels[k0][k1]['data'].shape
            newpts = sh[v0['axis'][0]]*sh[v0['axis'][1]]
            sh = tuple(np.r_[
                sh[:v0['axis'][0]],
                newpts,
                sh[v0['axis'][1]+1:]
            ].astype(int))
            
            newref = tuple(np.r_[
                v1['ref'][:v0['axis'][0]],
                [v0['dref']['npts']['key']],
                v1['ref'][v0['axis'][1]+1:],
            ])
    
            dlevels[k0][k1]['data'] = dlevels[k0][k1]['data'].swapaxes(
                v0['axis'][0],
                v0['axis'][1],
            ).reshape(sh)
            dlevels[k0][k1]['ref'] = newref

        dlevels[k0]['dref']['npts']['size'] = newpts
        
        if ii == 0:
            dax.add_ref(**dlevels[k0]['dref']['npts'])
            dax.add_ref(**dlevels[k0]['dref']['levels'])

        for k1 in ['cont0', 'cont1']:
            dax.add_data(**dlevels[k0][k1])

    # -----------
    # add contour

    kax = 'matrix'
    if dax.dax.get(kax) is not None:
        ax = dax.dax[kax]['handle']

        if ndim == 2:

            for ii, (k0, v0) in enumerate(dlevels.items()):
                ax.plot(
                    v0['cont0']['data'],
                    v0['cont1']['data'],
                    ls='-',
                    lw=1.,
                    c=dlevels[k0]['color'],
                )

        elif ndim == 3:

            for ii, (k0, v0) in enumerate(dlevels.items()):

                if v0['refZ'] is None:
                    ax.plot(
                        v0['cont0']['data'],
                        v0['cont1']['data'],
                        ls='-',
                        lw=1.,
                        c=dlevels[k0]['color'],
                    )
                    
                else:
                    # slice
                    sli = [
                        slice(None) if ii == v0['axis'][0] else 0
                        for ii in range(ndim-1)
                    ]
                    
                    # plot
                    l0, = ax.plot(
                        v0['cont0']['data'][sli],
                        v0['cont1']['data'][sli],
                        ls='-',
                        lw=1.,
                        c=dlevels[k0]['color'],
                    )
        
                    # store mobile
                    km = f'{k0}_contours'
                    dax.add_mobile(
                        key=km,
                        handle=l0,
                        # group_vis='Z',
                        # refs=[(refZ, ref_lvls), (refZ, ref_lvls)],
                        refs=[(v0['refZ'],), (v0['refZ'],)],
                        data=[v0['cont0']['key'], v0['cont1']['key']],
                        dtype=['xdata', 'ydata'],
                        axes=kax,
                        ind=0,
                    )

        else:
            raise NotImplementedError()

            # sli = [slice(None)]

            # l0, = ax.plot(
            # )

            # k0 = f'contour{ii}'
            # dax.add_mobile(
            # )
            
    # --------------------
    # add horizontal lines
    
    kax = 'radial'
    if dax.dax.get(kax) is not None and key in dlevels.keys():
        ax = dax.dax[kax]['handle']
            
        for ii, ll in enumerate(dlevels[key]['levels']):
            ax.axhline(ll, c='k', ls='--')
            

# #############################################################################
# #############################################################################
#                       Utilities
# #############################################################################


def _plot_bsplines_get_dx01(coll=None, km=None):
    # Get minimum distances

    wm = coll._which_mesh
    mtype = coll.dobj[wm][km]['type']
    if mtype == 'rect':
        knots0, knots1 = coll.dobj[wm][km]['knots']
        knots0 = coll.ddata[knots0]['data']
        knots1 = coll.ddata[knots1]['data']
        dx0 = np.min(np.diff(knots0))
        dx1 = np.min(np.diff(knots1))

    elif mtype == 'tri':
        indtri = coll.ddata[coll.dobj['mesh'][km]['ind']]['data']
        kknots = coll.dobj['mesh'][km]['knots']
        knots0 = coll.ddata[kknots[0]]['data']
        knots1 = coll.ddata[kknots[1]]['data']
        x0 = knots0[indtri]
        x1 = knots1[indtri]
        dist = np.mean(np.array([
            np.sqrt((x0[:, 1] - x0[:, 0])**2 + (x1[:, 1] - x1[:, 0])**2),
            np.sqrt((x0[:, 2] - x0[:, 1])**2 + (x1[:, 2] - x1[:, 1])**2),
            np.sqrt((x0[:, 2] - x0[:, 0])**2 + (x1[:, 2] - x1[:, 0])**2),
        ]))
        dx0, dx1 = dist, dist

    x0minmax = [knots0.min(), knots0.max()]
    x1minmax = [knots1.min(), knots1.max()]
    return dx0, dx1, x0minmax, x1minmax


# #############################################################################
# #############################################################################
#                       plot with submesh
# #############################################################################


def _plot_submesh(
    coll=None,
    coll2=None,
    key=None,
    keym=None,
    keybs=None,
    # figure
    vmin=None,
    vmax=None,
    cmap=None,
    dax=None,
    dmargin=None,
    fs=None,
    dcolorbar=None,
    dleg=None,
    interp=None,
    dkeys=None,
    lcol=None,
):

    if dax is None:
        dax = _plot_profile2d_submesh_create_axes(
            fs=fs,
            dmargin=dmargin,
        )

    # plot usual parts
    dax, dgroup = coll2.plot_as_array(
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        dax=dax,
        dmargin=dmargin,
        fs=fs,
        dcolorbar=dcolorbar,
        dleg=dleg,
        connect=False,
        interp=interp,
        label=True,
        **dkeys,
    )

    # ------------------
    # add radial profile to dax

    kradius, lkradial, lkdet, reft = _plot_profile2d_polar_add_radial(
        coll=coll,
        key=key,
        keym=keym,
        keybs=keybs,
        dax=dax,
    )

    assert (reft is not None) == ('Z' in dgroup.keys())
    if reft is not None and reft not in dgroup['Z']['ref']:
        dgroup['Z']['ref'].append(reft)
        dgroup['Z']['data'].append('index')

    # ------------------
    # add radial profile

    kax = 'radial'
    if dax.dax.get(kax) is not None:
        ax = dax.dax[kax]['handle']
        for ii in range(len(lkradial)):

            if reft is None:
                l0, = ax.plot(
                    dax.ddata[kradius]['data'],
                    dax.ddata[lkradial[ii]]['data'],
                    c=lcol[ii],
                    ls='-',
                    lw=2,
                )
            else:
                l0, = ax.plot(
                    dax.ddata[kradius]['data'],
                    dax.ddata[lkradial[ii]]['data'][0, :],
                    c=lcol[ii],
                    ls='-',
                    lw=2,
                )

                kl = f"radial{ii}"
                dax.add_mobile(
                    key=kl,
                    handle=l0,
                    refs=(reft,),
                    data=[lkradial[ii]],
                    dtype=['ydata'],
                    axes=kax,
                    ind=0,
                )

        if lkdet is not None:
            for ii in range(len(lkdet)):
                if reft is None:
                    l0, = ax.plot(
                        dax.ddata[kradius]['data'],
                        dax.ddata[lkdet[ii]]['data'],
                        ls='-',
                        lw=1,
                    )
                else:
                    l0, = ax.plot(
                        dax.ddata[kradius]['data'],
                        dax.ddata[lkdet[ii]]['data'][0, :],
                        ls='-',
                        lw=1,
                    )

                    kl = f"radial_det{ii}"
                    dax.add_mobile(
                        key=kl,
                        handle=l0,
                        refs=(reft,),
                        data=[lkdet[ii]],
                        dtype=['ydata'],
                        axes=kax,
                        ind=0,
                    )

        ax.set_xlim(
            dax.ddata[kradius]['data'].min(),
            dax.ddata[kradius]['data'].max(),
        )

        if vmin is not None:
            ax.set_ylim(bottom=vmin)
        if vmax is not None:
            ax.set_ylim(top=vmax)
            
    return dax, dgroup
            
            
def _plot_profile2d_polar_add_radial(
    coll=None,
    key=None,
    keym=None,
    keybs=None,
    dax=None,
):

    # -------------
    # key to radius
    
    wm = coll._which_mesh
    wbs = coll._which_bsplines
    
    kr2d = coll.dobj[wm][keym]['subkey'][0]
    kr = coll.dobj[wm][keym]['knots'][0]
    rr = coll.ddata[kr]['data']
    rad = np.linspace(rr[0], rr[-1], rr.size*10)

    # -----------------
    # get angle if any
    
    clas = coll.dobj[wbs][keybs]['class']
    # if clas.knotsa is None:
    if True:
        angle = None
    elif len(clas.shapebs) == 2:
        ka = coll.dobj[wbs][keybs]['apex'][1]
        angle = coll.ddata[ka]['data']
    elif np.sum(clas.nbs_a_per_r > 1) == 1:
        i0 = (clas.nbs_a_per_r > 1).nonzero()[0][0]
        angle = coll.dobj[wbs][keybs]['class'].apex_per_bs_a[i0]
    else:
        pass

    # angle
    if angle is None:
        radmap = rad
        anglemap = angle
    else:
        radmap = np.repeat(rad[:, None], angle.size, axis=1)
        anglemap = np.repeat(angle[None, :], rad.size, axis=0)

    # ----
    # reft

    refc = [rr for rr in coll.ddata[key]['ref'] if rr in coll.ddata[kr2d]['ref']]
    if len(refc) == 1:
        refc = refc[0]
    else:
        refc = None

    reft, keyt, _, dind = coll.get_ref_vector_common(
        keys=[key, kr2d],
        ref=refc,
        dim=None,
        quant=None,
        name=None,
        units=None,
        values=None,
        indices=None,
        ind_strict=None,
    )[1:]

    # radial total profile
    # radial, t_radial, _ = coll.interpolate(
    dout = coll.interpolate(
        keys=key,
        ref_key=keybs,
        x0=radmap,
        x1=anglemap,
        grid=False,
        # ref_com=keyt,
    )[key]

    radial = dout['data']
    # if reft is not None and radial.ndim == radmap.ndim:
    #     radial = np.repeat(radial[None, ...], t_radial.size, axis=0)

    # -------------------------------
    # details for purely-radial cases
    
    if angle is None:
        # radial_details, t_radial, _ = coll.interpolate(
        dout_details = coll.interpolate(
            ref_key=keybs,
            x0=rad,
            grid=False,
            details=True,
        )[f'{keybs}_details']

        radial_details = dout_details['data']
        if reft is None:
            radial_details = radial_details * coll.ddata[key]['data'][None, :]
            refdet = ('nradius',)
        else:
            refdet = (reft, 'nradius')
            if reft == coll.get_time(key)[2]:
                radial_details = (
                    radial_details[None, :, :]
                    * coll.ddata[key]['data'][:, None, :]
                )
            elif key in dind.keys():
                radial_details = (
                    radial_details[None, :, :]
                    * coll.ddata[key]['data'][dind[key]['ind'], None, :]
                )

        nbs = radial_details.shape[-1]

    # -----------
    # add to dax
    
    dax.add_ref(key='nradius', size=rad.size)
    if angle is not None:
        dax.add_ref(key='nangle', size=angle.size)

    if reft is not None:
        assert radial.ndim > 1 and radial.shape[0] > 1
        # if angle is not None:
            # ref = (reft, 'nangle', 'nradius')
        # else:
        ref = (reft, 'nradius')
    else:
        # if angle is not None:
            # ref = ('nangle', 'nradius')
        # else:
        ref = 'nradius'

    # ------------
    # add to ddata
    
    kradius = 'radius'
    dax.add_data(key=kradius, data=rad, ref='nradius')
    if angle is None:
        lk = ['radial']
        dax.add_data(key=lk[0], data=radial, ref=ref)
        lkdet = [f'radial-detail-{ii}' for ii in range(nbs)]
        for ii in range(nbs):
            dax.add_data(
                key=lkdet[ii], data=radial_details[..., ii], ref=refdet,
            )

    else:
        kangle = 'angle'
        dax.add_data(key=kangle, data=angle, ref='nangle')
        lkdet = None
        lk = [f'radial-{ii}' for ii in range(angle.size)]
        if reft is None:
            for ii in range(angle.size):
                dax.add_data(key=lk[ii], data=radial[:, ii], ref=ref)
        else:
            for ii in range(angle.size):
                dax.add_data(key=lk[ii], data=radial[:, :, ii], ref=ref)

    return kradius, lk, lkdet, reft


# ##############################################################
# ##############################################################
#          create default axes for submesh
# ##############################################################


def _plot_profile2d_submesh_create_axes(
    fs=None,
    dmargin=None,
):

    if fs is None:
        fs = (15, 9)

    if dmargin is None:
        dmargin = {
            'left': 0.05, 'right': 0.95,
            'bottom': 0.05, 'top': 0.95,
            'hspace': 0.4, 'wspace': 0.3,
        }

    fig = plt.figure(figsize=fs)
    gs = gridspec.GridSpec(ncols=6, nrows=6, **dmargin)

    # axes for image
    ax0 = fig.add_subplot(gs[:4, 2:4], aspect='auto')

    # axes for vertical profile
    ax1 = fig.add_subplot(gs[:4, 4], sharey=ax0)

    # axes for horizontal profile
    ax2 = fig.add_subplot(gs[4:, 2:4], sharex=ax0)

    # axes for traces
    ax3 = fig.add_subplot(gs[2:4, :2])

    # axes for traces
    ax7 = fig.add_subplot(gs[:2, :2], sharey=ax2)

    # axes for text
    ax4 = fig.add_subplot(gs[:3, 5], frameon=False)
    ax5 = fig.add_subplot(gs[3:, 5], frameon=False)
    ax6 = fig.add_subplot(gs[4:, :2], frameon=False)

    # dax
    dax = {
        # data
        'matrix': {'handle': ax0, 'type': 'matrix'},
        'vertical': {'handle': ax1, 'type': 'misc'},
        'horizontal': {'handle': ax2, 'type': 'misc'},
        'traces': {'handle': ax3, 'type': 'misc'},
        'radial': {'handle': ax7, 'type': 'misc'},
        # text
        'textX': {'handle': ax4, 'type': 'text'},
        'textY': {'handle': ax5, 'type': 'text'},
        'textZ': {'handle': ax6, 'type': 'text'},
    }
    return dax