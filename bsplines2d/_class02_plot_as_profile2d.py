# -*- coding: utf-8 -*-


# Common
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datastock as ds

# specific


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
    dunique_mesh_2d=None,
    # contours
    dlevels=None,
    ref_com=None,
    # details
    plot_details=None,
    # ref vectors
    dref_vectorZ=None,
    dref_vectorU=None,
    ref_vector_strategy=None,
    uniform=None,
    # interpolation
    val_out=None,
    nan0=None,
    # figure
    dvminmax=None,
    # vmin=None,
    # vmax=None,
    cmap=None,
    dax=None,
    dmargin=None,
    fs=None,
    dcolorbar=None,
    dleg=None,
    # interactivity
    dinc=None,
    connect=None,
    show_commands=None,
):

    # --------------
    # check input

    if connect is None:
        connect = True

    (
        dkeys,
        dlevels,
        cmap, _,
        dvminmax,
        dcolorbar, dleg,
        connect,
    ) = _check(
        coll=coll,
        keys=key,
        dlevels=dlevels,
        # plotting
        dvminmax=dvminmax,
        # vmin=vmin,
        # vmax=vmax,
        cmap=cmap,
        # figure
        dcolorbar=dcolorbar,
        dleg=dleg,
        # interactivity
        connect=connect,
    )

    # --------------
    #  Prepare data

    (
        collax, dkeys,
        dlevels, lcol,
    ) = _prepare(
        coll=coll,
        dkeys=dkeys,
        dres=dres,
        dunique_mesh_2d=dunique_mesh_2d,
        # interpolation
        val_out=val_out,
        nan0=nan0,
        # levels
        dlevels=dlevels,
        ref_com=ref_com,
        # ref vectors
        dref_vectorZ=dref_vectorZ,
        dref_vectorU=dref_vectorU,
    )

    lkeys = ['key', 'keyX', 'keyY', 'keyZ', 'keyU']

    # --------------
    #  Prepare dax

    daxi = None
    if dax is not None:
        if not all([k0 in dax.keys() for k0 in dkeys.keys()]):
            dax = {k0: dax for k0 in dkeys.keys()}

    # nmax, color
    if len(dkeys) == 1:
        nmax = None
        color_dict = {k0: None for k0 in dkeys.keys()}
    else:
        nmax = 1
        color_dict = {k0: lcol[ii] for ii, k0 in enumerate(dkeys.keys())}

    # -----------------
    # loop on profile 2d

    for k0, v0 in dkeys.items():

        # -----------------
        # case with submesh

        if dax is not None:
            daxi = dax[k0]

        if v0['submesh'] is not None:
            collax, dgroup = _plot_submesh(
                coll=coll,
                collax=collax,
                key=k0,
                keym=v0['keym'],
                keybs=v0['keybs'],
                # ref vector
                dref_vector=dref_vectorZ,
                ref_vector_strategy=ref_vector_strategy,
                uniform=uniform,
                # details
                plot_details=plot_details,
                # plotting
                dvminmax=dvminmax,
                # vmin=vmin,
                # vmax=vmax,
                cmap=cmap,
                color_dict=color_dict[k0],
                # figure
                dax=daxi,
                dmargin=dmargin,
                fs=fs,
                dcolorbar=dcolorbar,
                dleg=dleg,
                interp=v0['interp'],
                dkeys={k1: v0[k1] for k1 in lkeys},
                lcol=lcol,
            )

        # -------------------
        # without submesh

        else:

            collax, dgroup = collax.plot_as_array(
                dvminmax=dvminmax,
                # vmin=vmin,
                # vmax=vmax,
                cmap=cmap,
                dax=daxi,
                dmargin=dmargin,
                fs=fs,
                dcolorbar=dcolorbar,
                dleg=dleg,
                uniform=uniform,
                interp=v0['interp'],
                color_dict=color_dict[k0],
                nmax=nmax,
                connect=False,
                inplace=True,
                **{k1: v0[k1] for k1 in lkeys},
            )

        # -----------
        # add levels

        if dlevels is not None:
            _add_levels_2d(
                collax=collax,
                dgroup=dgroup,
                dlevels=dlevels,
            )

    # -----------
    # connect

    if connect is True:
        collax.setup_interactivity(kinter='inter0', dgroup=dgroup, dinc=dinc)
        collax.disconnect_old()
        collax.connect()

        collax.show_commands(verb=show_commands)
        return collax
    else:
        return collax, dgroup


# #############################################################################
# #############################################################################
#                       Check
# #############################################################################


def _check(
    coll=None,
    keys=None,
    dlevels=None,
    # figure
    dvminmax=None,
    # vmin=None,
    # vmax=None,
    cmap=None,
    cmap_err=None,
    dcolorbar=None,
    dleg=None,
    # interactivity
    connect=None,
):

    # ----------
    # keys

    # key
    dk = coll.get_profiles2d()

    if isinstance(keys, str):
        keys = [keys]
    keys = ds._generic_check._check_var_iter(
        keys, 'keys',
        types=(list, tuple),
        types_iter=str,
        allowed=list(dk.keys()),
    )

    # check refs
    dkeys = {}
    wm = coll._which_mesh
    wbs = coll._which_bsplines
    for ii, k0 in enumerate(keys):

        keybs = dk[k0]
        keym = coll.dobj[wbs][keybs][wm]
        nd = coll.dobj[wm][keym]['nd']
        mtype = coll.dobj[wm][keym]['type']

        # submesh
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

        # ref_other
        ref_other = [
            rr for rr in coll.ddata[k0]['ref']
            if rr not in coll.dobj[wbs][keybs]['ref']
        ]
        if ii == 1:
            if ref_other != dkeys[keys[0]]['ref_other']:
                lstr = [f"\t- {kk}: {coll.ddata[kk]['ref']}" for kk in keys]
                msg = (
                    "Provided keys must be profile2d with identical ref!"
                    " (apart from the actual profile2d bsplines)\n"
                    "\t- Provided:"
                    + "\n".join(lstr)
                )
                raise Exception(msg)

        dkeys[k0] = {
            'keybs': keybs,
            'keym': keym,
            'nd': nd,
            'mtype': mtype,
            'submesh': submesh,
            'subbs': subbs,
            'submtype': submtype,
            'ref_other': ref_other,
        }

    # ----------
    # dlevels

    if dlevels is not None:

        dp2d = coll.get_profiles2d()
        if isinstance(dlevels, (float, int)):
            dlevels = {k0: np.r_[dlevels] for k0 in keys}

        elif isinstance(dlevels, (np.ndarray, list, tuple)):
            dlevels = np.atleast_1d(dlevels).ravel()
            dlevels = {k0: dlevels for k0 in keys}

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
            dlevels[k0]['ls'] = v0.get('ls', '-')
            dlevels[k0]['lw'] = v0.get('lw', 1)

    # ----------
    # cmap

    # cmap
    if cmap is None:
        cmap = 'viridis'

    # cmap_err
    if cmap_err is None:
        cmap_err = 'seismic'

    # ----------
    # vmin, vmax

    if dvminmax is None:
        dvminmax = {
            'data': {
                'min': min([0] + [np.nanmin(coll.ddata[kk]['data']) for kk in keys]),
                'max': max([np.nanmax(coll.ddata[kk]['data']) for kk in keys]),
            },
        }

    # ----------
    # figure

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
        dkeys,
        dlevels,
        cmap, cmap_err,
        dvminmax,
        dcolorbar, dleg,
        connect,
    )


# #############################################################################
# #############################################################################
#                       Prepare
# #############################################################################


def _prepare(
    coll=None,
    dkeys=None,
    dres=None,
    dunique_mesh_2d=None,
    # levels
    dlevels=None,
    ref_com=None,
    # interpolation
    val_out=None,
    nan0=None,
    # ref vectors
    dref_vectorZ=None,
    dref_vectorU=None,
):

    # -----------
    # check

    if dref_vectorZ is None:
        dref_vectorZ = {}

    if dref_vectorU is None:
        dref_vectorU = {}

    # ---------------------------
    # get interpolated collection

    coll2, dbs = coll.interpolate_all_bsplines(
        keys=list(dkeys.keys()),
        dres=dres,
        dunique_mesh_2d=dunique_mesh_2d,
        submesh=True,
        # interpolation
        val_out=val_out,
        nan0=nan0,
        # for plotting => uniform
        for_plotting=True,
    )
    lbs2d = [k0 for k0, v0 in dbs.items() if len(v0['ref']) == 2]

    # -----------------
    # get plotting mesh

    for k0, v0 in dkeys.items():

        dkeys[k0].update(_get_dkey(
            coll=coll,
            subbs=v0['subbs'],
            kdata=k0,
            lbs2d=lbs2d,
            dbs=dbs,
            coll2=coll2,
            dref_vectorZ=dref_vectorZ,
            dref_vectorU=dref_vectorU,
        ))

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
                dlevels[k0][k1]['key'] = f'{k0}_{k1}'

            dlevels[k0]['refZ'] = refZ
            dlevels[k0]['refU'] = refU
            dlevels[k0]['axis'] = axis

    # --------------------
    # optional for submesh

    lcol = ['k', 'b', 'g', 'r', 'm', 'c', 'y']

    return (
        coll2, dkeys,
        dlevels, lcol,
    )


def _get_dkey(
    coll=None,
    subbs=None,
    kdata=None,
    lbs2d=None,
    dbs=None,
    coll2=None,
    dref_vectorZ=None,
    dref_vectorU=None,
):

    # deg and interp
    wbs = coll._which_bsplines
    deg = coll.dobj[wbs][subbs]['deg']
    if deg == 0:
        interp = 'nearest'
    elif deg == 1:
        interp = 'bilinear'
    elif deg >= 2:
        interp = 'bicubic'

    # key X, Y, Z, U
    bs2d = [k1 for k1 in lbs2d if k1 in coll.ddata[kdata][wbs]][0]
    rX, rY = dbs[bs2d]['ref']
    lr1d = [k1 for k1 in coll2.ddata[kdata]['ref'] if k1 not in [rX, rY]]
    ndim = coll2.ddata[kdata]['data'].ndim

    # higher dimensions
    keyZ, keyU = None, None
    if ndim >= 3:
        keyZ = coll2.get_ref_vector(
            ref=lr1d[0],
            **dref_vectorZ,
        )[3]
        # uniform = ds._plot_as_array._check_uniform_lin(
            # k0=keyZ, ddata=coll2.ddata,
        # )
        # if not uniform:
            # keyZ = None
        if ndim == 4:
            keyU = coll2.get_ref_vector(
                ref=lr1d[1],
                **dref_vectorU,
            )[3]

    return {
        'deg': deg,
        'interp': interp,
        'key': kdata,
        'keyX': coll2.get_ref_vector(ref=rX)[3],
        'keyY': coll2.get_ref_vector(ref=rY)[3],
        'keyZ': keyZ,
        'keyU': keyU,
    }


# #############################################################################
# #############################################################################
#                   Plot fixed mesh with levels
# #############################################################################


def _add_levels_2d(
    collax=None,
    key=None,
    dgroup=None,
    dlevels=None,
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
            collax.add_ref(**dlevels[k0]['dref']['npts'])
            collax.add_ref(**dlevels[k0]['dref']['levels'])

        for k1 in ['cont0', 'cont1']:
            collax.add_data(**dlevels[k0][k1])

    # -----------
    # add contour

    kax = 'matrix'
    if collax.dax.get(kax) is not None:
        ax = collax.dax[kax]['handle']

        if ndim == 2:

            for ii, (k0, v0) in enumerate(dlevels.items()):
                ax.plot(
                    v0['cont0']['data'],
                    v0['cont1']['data'],
                    lw=dlevels[k0]['lw'],
                    c=dlevels[k0]['color'],
                    ls=dlevels[k0]['ls'],
                )

        elif ndim == 3:

            for ii, (k0, v0) in enumerate(dlevels.items()):

                if v0['refZ'] is None:
                    ax.plot(
                        v0['cont0']['data'],
                        v0['cont1']['data'],
                        lw=dlevels[k0]['lw'],
                        c=dlevels[k0]['color'],
                        ls=dlevels[k0]['ls'],
                    )

                else:
                    # slice
                    sli = tuple([
                        slice(None) if ii == v0['axis'][0] else 0
                        for ii in range(ndim-1)
                    ])

                    # plot
                    l0, = ax.plot(
                        v0['cont0']['data'][sli],
                        v0['cont1']['data'][sli],
                        lw=dlevels[k0]['lw'],
                        c=dlevels[k0]['color'],
                        ls=dlevels[k0]['ls'],
                    )

                    # store mobile
                    km = f'{k0}_contours'
                    collax.add_mobile(
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
            # collax.add_mobile(
            # )

    # --------------------
    # add horizontal lines

    kax = 'radial'
    if collax.dax.get(kax) is not None and key in dlevels.keys():
        ax = collax.dax[kax]['handle']

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
    collax=None,
    key=None,
    keym=None,
    keybs=None,
    # ref vetcor
    dref_vector=None,
    ref_vector_strategy=None,
    uniform=None,
    # plot_details
    plot_details=None,
    # figure
    dvminmax=None,
    # vmin=None,
    # vmax=None,
    cmap=None,
    color_dict=None,
    dax=None,
    dmargin=None,
    fs=None,
    dcolorbar=None,
    dleg=None,
    interp=None,
    dkeys=None,
    lcol=None,
):

    # ------------
    # check inputs

    plot_details = ds._generic_check._check_var(
        plot_details, 'plot_details',
        types=bool,
        default=False,
    )

    if dax is None:
        dax = _plot_profile2d_submesh_create_axes(
            fs=fs,
            dmargin=dmargin,
        )

    # ------------
    # preliminary

    # plot usual parts
    collax, dgroup = collax.plot_as_array(
        dvminmax=dvminmax,
        # vmin=vmin,
        # vmax=vmax,
        cmap=cmap,
        color_dict=color_dict,
        dax=dax,
        dmargin=dmargin,
        fs=fs,
        dcolorbar=dcolorbar,
        dleg=dleg,
        connect=False,
        interp=interp,
        label=True,
        inplace=True,
        uniform=uniform,
        **dkeys,
    )

    # ------------------
    # add radial profile to dax

    kradius, lkradial, lkdet, reft = _plot_profile2d_polar_add_radial(
        coll=coll,
        key=key,
        keym=keym,
        keybs=keybs,
        collax=collax,
        dref_vector=dref_vector,
        ref_vector_strategy=ref_vector_strategy,
        plot_details=plot_details,
    )

    if (reft is not None) != ('Z' in dgroup.keys()):
        msg = (
            "(reft is not None) != ('Z' in dgroup.keys()):\n"
            f"reft = {reft}\n"
            f"dgroup.keys() = {dgroup.keys()}\n"
            f"dref_vector = {dref_vector}"
        )
        raise Exception(msg)

    if reft is not None and reft not in dgroup['Z']['ref']:
        dgroup['Z']['ref'].append(reft)
        dgroup['Z']['data'].append('index')

    # ------------------
    # add radial profile

    kax = 'radial'
    lax = [k0 for k0, v0 in collax.dax.items() if kax in v0['type']]
    if len(lax) == 1:
        kax = lax[0]
        ax = collax.dax[kax]['handle']
        for ii in range(len(lkradial)):

            if reft is None:
                l0, = ax.plot(
                    collax.ddata[kradius]['data'],
                    collax.ddata[lkradial[ii]]['data'],
                    ls='-',
                    lw=2,
                    c=color_dict,
                )

            elif collax.ddata[lkradial[ii]]['data'].ndim == 2:
                l0, = ax.plot(
                    collax.ddata[kradius]['data'],
                    collax.ddata[lkradial[ii]]['data'][0, :],
                    ls='-',
                    lw=2,
                    c=color_dict,
                )

                kl = f"{key}_radial{ii}"
                collax.add_mobile(
                    key=kl,
                    handle=l0,
                    refs=(reft,),
                    data=[lkradial[ii]],
                    dtype=['ydata'],
                    axes=kax,
                    ind=0,
                )
            else:
                msg = (
                    "spectral, radial and time-dependent profile2d"
                    " plotting not implemented yet!"
                )
                raise NotImplementedError(msg)

        if lkdet is not None:
            for ii in range(len(lkdet)):
                if reft is None:
                    l0, = ax.plot(
                        collax.ddata[kradius]['data'],
                        collax.ddata[lkdet[ii]]['data'],
                        ls='-',
                        lw=1,
                    )
                else:
                    l0, = ax.plot(
                        collax.ddata[kradius]['data'],
                        collax.ddata[lkdet[ii]]['data'][0, :],
                        ls='-',
                        lw=1,
                    )

                    kl = f"{key}_radial_det{ii}"
                    collax.add_mobile(
                        key=kl,
                        handle=l0,
                        refs=(reft,),
                        data=[lkdet[ii]],
                        dtype=['ydata'],
                        axes=kax,
                        ind=0,
                    )

        ax.set_xlim(
            collax.ddata[kradius]['data'].min(),
            collax.ddata[kradius]['data'].max(),
        )

        if dvminmax.get('data', {}).get('min') is not None:
            ax.set_ylim(bottom=dvminmax['data']['min'])
        if  dvminmax.get('data', {}).get('min') is not None:
            ax.set_ylim(top=dvminmax['data']['max'])

    return collax, dgroup


def _plot_profile2d_polar_add_radial(
    coll=None,
    key=None,
    keym=None,
    keybs=None,
    collax=None,
    # ref_vector
    dref_vector=None,
    ref_vector_strategy=None,
    # details
    plot_details=None,
):

    # ----------
    # check

    if dref_vector is None:
        dref_vector = {}

    # -------------
    # key to radius

    wm = coll._which_mesh
    wbs = coll._which_bsplines

    kr2d = coll.dobj[wm][keym]['subkey'][0]
    kr = coll.dobj[wm][keym]['knots'][0]

    # special case of deg = 1 => use knots directly
    rr = coll.ddata[kr]['data']
    if coll.dobj[wbs][keybs]['deg'] == 1:
        rad = rr
    else:
        rad = np.linspace(rr[0], rr[-1], 10*rr.size)

    # temporary keys
    refr = f'{key}_nradius'
    kradius = f'{key}_radius'

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

    # find reft
    reft, keyt, _, dind = coll.get_ref_vector_common(
        keys=[key, kr2d],
        ref=refc,
        strategy=ref_vector_strategy,
        **dref_vector,
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

    if angle is None and plot_details is True:
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
            refdet = (refr,)

        else:
            refdet = (reft, refr)
            if reft == coll.get_time(key)[2]:
                radial_details = (
                    radial_details[None, :, :]
                    * coll.ddata[key]['data'][:, None, :]
                )

            elif key in dind.keys():

                if dind[key].get('ind') is None:
                    radial_details = (
                        radial_details[None, :, :]
                        * coll.ddata[key]['data'][:, None, :]
                    )

                else:
                    radial_details = (
                        radial_details[None, :, :]
                        * coll.ddata[key]['data'][dind[key]['ind'], None, :]
                    )


        nbs = radial_details.shape[-1]

    # -------------
    # add to collax

    collax.add_ref(key=refr, size=rad.size)
    if angle is not None:
        collax.add_ref(key='nangle', size=angle.size)

    ref = list(dout['ref'])
    ref[ref.index(None)] = refr
    ref = tuple(ref)

    # ------------
    # add to ddata

    lkdet = None
    collax.add_data(key=kradius, data=rad, ref=refr)
    if angle is None:
        lk = [f'{key}_radial']
        collax.add_data(key=lk[0], data=radial, ref=ref)
        if plot_details is True:
            lkdet = [f'{key}_radial_detail_{ii}' for ii in range(nbs)]
            for ii in range(nbs):
                collax.add_data(
                    key=lkdet[ii], data=radial_details[..., ii], ref=refdet,
                )

    else:
        kangle = 'angle'
        collax.add_data(key=kangle, data=angle, ref='nangle')
        lkdet = None
        lk = [f'radial-{ii}' for ii in range(angle.size)]
        if reft is None:
            for ii in range(angle.size):
                collax.add_data(key=lk[ii], data=radial[:, ii], ref=ref)
        else:
            for ii in range(angle.size):
                collax.add_data(key=lk[ii], data=radial[:, :, ii], ref=ref)

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
        'matrix': {'handle': ax0},
        'vertical': {'handle': ax1},
        'horizontal': {'handle': ax2},
        'tracesZ': {'handle': ax3},
        'radial': {'handle': ax7},
        # text
        'textX': {'handle': ax4},
        'textY': {'handle': ax5},
        'textZ': {'handle': ax6},
    }
    return dax