# -*- coding: utf-8 -*-


# Common
import numpy as np
import datastock as ds
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# ##################################################
# ##################################################
#                       Main
# ##################################################


def main(
    coll=None,
    key=None,
    res=None,
    # slice
    Z=None,
    phi=None,
    # domain
    DR=None,
    DZ=None,
    Dphi=None,
    # option
    reshape_2d=None,
    # plot
    plot=None,
    dax=None,
    color=None,
):

    # --------------
    # check inputs
    # --------------

    (
        res, Z, phi, domain, reshape_2d, plot,
    ) = _check(
        res=res,
        # slice
        Z=Z,
        phi=phi,
        # domain
        DR=DR,
        DZ=DZ,
        Dphi=Dphi,
        # option
        reshape_2d=reshape_2d,
        plot=plot,
    )

    # --------------
    # prepare
    # --------------

    (
        func_RZphi_from_ind,
        func_ind_from_domain,
    ) = coll.get_sample_mesh_3d_func(
        key=key,
        res_RZ=res,
        res_phi=res,
        mode='abs',
    )

    # --------------
    # compute
    # --------------

    if phi is None:
        indr, indz, indphi = _horizontal_slice(
            res=res,
            # sampling
            func_ind_from_domain=func_ind_from_domain,
            func_RZphi_from_ind=func_RZphi_from_ind,
            # slice
            Z=Z,
            # domain
            **domain,
        )

    else:
        indr, indz, indphi = _poloidal_slice(
            res=res,
            # sampling
            func_ind_from_domain=func_ind_from_domain,
            func_RZphi_from_ind=func_RZphi_from_ind,
            # slice
            phi=phi,
            # option
            reshape_2d=reshape_2d,
            # domain
            **domain,
        )

    # ------------
    # get points
    # ------------

    pts_r, pts_z, pts_phi, _ = func_RZphi_from_ind(
        indr=indr,
        indz=indz,
        indphi=indphi,
    )

    # --------------
    # output
    # --------------

    dout = {
        # indices
        'indr': {
            'data': indr,
            'units': 'index',
        },
        'indz': {
            'data': indz,
            'units': 'index',
        },
        'indphi': {
            'data': indphi,
            'units': 'index',
        },
        # points
        'pts_r': {
            'data': pts_r,
            'units': 'm',
        },
        'pts_z': {
            'data': pts_z,
            'units': 'm',
        },
        'pts_phi': {
            'data': pts_phi,
            'units': 'rad',
        },
    }

    # -----------------
    # optional plotting
    # -----------------

    if plot is True:
        _plot(
            dout,
            dax=dax,
            color=color,
        )

    return dout


# ##############################################
# ##############################################
#                   Check
# ##############################################


def _check(
    res=None,
    # slice
    Z=None,
    phi=None,
    # domain
    DR=None,
    DZ=None,
    Dphi=None,
    # option
    reshape_2d=None,
    # plot
    plot=None,
    color=None,
):

    # ---------
    # res
    # ---------

    res = float(ds._generic_check._check_var(
        res, 'res',
        types=(int, float),
        sign=">0",
    ))

    # ---------
    # Z vs phi
    # ---------

    lc = [
        Z is not None,
        phi is not None,
    ]
    if np.sum(lc) != 1:
        msg = (
            "Please provide Z xor phi!\n"
            "\t- Z: height of horizontal slice\n"
            "\t- phi: toroidal angle of poloidal slice\n"
            "Provided:\n"
            f"\t- Z = {Z}\n"
            f"\t- phi = {phi}\n"
        )
        raise Exception(msg)

    if phi is None:
        Z = float(ds._generic_check._check_var(
            Z, 'Z',
            types=(int, float),
        ))
    else:
        phi = float(ds._generic_check._check_var(
            phi, 'phi',
            types=(int, float),
        ))
        phi = np.arctan2(np.sin(phi), np.cos(phi))

    # ---------
    # domain
    # ---------

    domain = {
        'DR': DR,
        'DZ': DZ,
        'Dphi': Dphi,
    }

    for k0, v0 in domain.items():
        if v0 is not None:
            v0 = ds._generic_check._check_flat1darray(
                v0, k0,
                dtype=float,
                unique=True,
                size=2,
            )

            if v0[0] >= v0[1]:
                msg = (
                    f"If provided, arg '{k0}' must be strictly increasing!\n"
                    + "For Dphi, use +2pi if needed\n" if k0 == 'Dphi' else ''
                )
                raise Exception(msg)

    if phi is None:
        del domain['DZ']
    else:
        del domain['Dphi']

    # ---------
    # reshape_2d
    # ---------

    reshape_2d = ds._generic_check._check_var(
        reshape_2d, 'reshape_2d',
        types=bool,
        default=True,
    )

    # ---------
    # plot
    # ---------

    plot = ds._generic_check._check_var(
        plot, 'plot',
        types=bool,
        default=False,
    )

    return res, Z, phi, domain, reshape_2d, plot


# ##############################################
# ##############################################
#               Horizontal slice
# ##############################################


def _horizontal_slice(
    res=None,
    # pre-sampling
    func_ind_from_domain=None,
    func_RZphi_from_ind=None,
    # slice
    Z=None,
    # domain
    DR=None,
    Dphi=None,
):
    # -----------
    # domain Z
    # -----------

    dZ = 1.5*res
    DZ = (Z - dZ, Z + dZ)

    indr, indz, indphi = func_ind_from_domain(
        DR=DR,
        DZ=DZ,
        Dphi=Dphi,
    )

    # ------------
    # select plane
    # ------------

    izu = np.unique(indz)
    _, zz, _, _ = func_RZphi_from_ind(indz=izu)

    iz = np.argmin(np.abs(zz - Z))
    iok = indz == izu[iz]

    return indr[iok], indz[iok], indphi[iok]


# ##############################################
# ##############################################
#               Poloidal slice
# ##############################################


def _poloidal_slice(
    res=None,
    # pre-sampling
    func_ind_from_domain=None,
    func_RZphi_from_ind=None,
    # slice
    phi=None,
    # domain
    DR=None,
    DZ=None,
    # option
    reshape_2d=None,
):
    # -----------
    # domain Z
    # -----------

    dphi = np.pi/8
    Dphi = (phi - dphi, phi + dphi)

    indr, indz, indphi = func_ind_from_domain(
        DR=DR,
        DZ=DZ,
        Dphi=Dphi,
    )

    # ------------
    # select plane
    # ------------

    iru = np.unique(indr)
    iphi = np.zeros(iru.shape)
    indr_new, indz_new, indphi_new = [], [], []
    for ii, iri in enumerate(iru):
        ind = indr == iri
        rr, _, phii, dV = func_RZphi_from_ind(
            indr=indr[ind],
            indphi=indphi[ind]
        )

        iphi[ii] = indphi[ind][np.argmin(np.abs(phii - phi))]

        ind[ind] = indphi[ind] == iphi[ii]
        indr_new.append(indr[ind])
        indz_new.append(indz[ind])
        indphi_new.append(indphi[ind])

    # ---------
    # reshape
    # ---------
    ir = np.concatenate(indr_new)
    iz = np.concatenate(indz_new)
    iphi = np.concatenate(indphi_new)
    if reshape_2d is True:
        i0u = np.unique(ir)
        i1u = np.unique(iz)

        indr = -np.ones((i0u.size, i1u.size), dtype=int)
        indz = -np.ones((i0u.size, i1u.size), dtype=int)
        indphi = -np.ones((i0u.size, i1u.size), dtype=int)

        for ii, iri in enumerate(i0u):
            ind = ir == iri

            indsz = np.argsort(iz[ind])
            iiz = np.searchsorted(i1u, np.sort(iz[ind]))
            sli = (iri, iiz)

            indr[sli] = iri
            indz[sli] = iz[ind][indsz]
            indphi[sli] = iphi[ind][indsz]

    else:
        indr = ir
        indz = iz
        indphi = iphi

    return indr, indz, indphi


# ##############################################
# ##############################################
#               Plot
# ##############################################


def _plot(
    dout=None,
    dax=None,
    color=None,
):
    # --------------
    # prepare data
    # --------------

    iok = dout['indr']['data'] >= 0
    xx = dout['pts_r']['data'][iok] * np.cos(dout['pts_phi']['data'][iok])
    yy = dout['pts_r']['data'][iok] * np.sin(dout['pts_phi']['data'][iok])

    # --------------
    # prepare figure
    # --------------

    if dax is None:

        dmargin = {
            'left': 0.05, 'right': 0.95,
            'bottom': 0.06, 'top': 0.90,
            'hspace': 0.20, 'wspace': 0.20,
        }

        fig = plt.figure(figsize=(14, 8))
        gs = gridspec.GridSpec(ncols=3, nrows=1, **dmargin)

        # --------------
        # prepare axes

        ax0 = fig.add_subplot(gs[:, 0], aspect='equal', adjustable='box')
        ax0.set_ylabel('Z (m)', size=12, fontweight='bold')
        ax0.set_xlabel('R (m)', size=12, fontweight='bold')
        ax0.set_title("cross-section", size=14, fontweight='bold')

        ax1 = fig.add_subplot(gs[:, 1], aspect='equal', adjustable='box')
        ax1.set_xlabel('X (m)', size=12, fontweight='bold')
        ax1.set_ylabel('Y (m)', size=12, fontweight='bold')
        ax1.set_title("horizontal", size=14, fontweight='bold')

        ax2 = fig.add_subplot(
            gs[:, 2],
            aspect='equal',
            adjustable='box',
            projection='3d',
        )
        ax2.set_xlabel('X (m)', size=12, fontweight='bold')
        ax2.set_ylabel('Y (m)', size=12, fontweight='bold')
        ax2.set_zlabel('Z (m)', size=12, fontweight='bold')
        ax2.set_title("3d", size=14, fontweight='bold')

        dax = {
            'cross': {'handle': ax0},
            'hor': {'handle': ax1},
            '3d': {'handle': ax2},
        }

    # --------------
    # plot cross
    # --------------

    kax = 'cross'
    if dax.get(kax) is not None:
        ax = dax[kax]['handle']

        ax.plot(
            dout['pts_r']['data'][iok],
            dout['pts_z']['data'][iok],
            marker='.',
            linestyle='None',
            ms=6,
            color=color,
        )

    # --------------
    # plot hor
    # --------------

    kax = 'hor'
    if dax.get(kax) is not None:
        ax = dax[kax]['handle']

        ax.plot(
            xx,
            yy,
            marker='.',
            linestyle='None',
            ms=6,
            color=color,
        )

    # --------------
    # plot 3d
    # --------------

    kax = '3d'
    if dax.get(kax) is not None:
        ax = dax[kax]['handle']

        ax.plot(
            xx,
            yy,
            dout['pts_z']['data'][iok],
            marker='.',
            linestyle='None',
            ms=6,
            color=color,
        )

    return
