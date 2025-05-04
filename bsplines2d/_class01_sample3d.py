# -*- coding: utf-8 -*-


# Common
import numpy as np
import datastock as ds
import matplotlib.pyplot as plt
from matplotlib.path import Path


# tofu
from . import _generic_mesh


# ##############################################################
# ##############################################################
#                       Main
# ##############################################################


def main(
    coll=None,
    key=None,
    res_RZ=None,
    mode=None,
    res_phi=None,
):

    # -------------
    # check inputs

    (
        key, res_phi,
    ) = _check(
        coll=coll,
        key=key,
        res_phi=res_phi,
    )

    # ---------------------
    # prepare dsamp
    # ---------------------

    dsamp = coll.get_sample_mesh(
        key=key,
        res=res_RZ,
        grid=False,
        mode=mode,
        x0=None,
        x1=None,
        Dx0=None,
        Dx1=None,
        imshow=False,
        in_mesh=True,
        store=False,
    )

    x0u = dsamp['x0']['data']
    x1u = dsamp['x1']['data']
    index = dsamp['ind']['data']

    # ------------------
    # phi from res_phi
    # ------------------

    nphi = np.ceil(3 + x0u * (2*np.pi) / res_phi).astype(int)

    # -----------------
    # func_pts_from_ind
    # -----------------

    func_RZphi_from_ind = _get_func_RZphi_from_ind(
        x0u=x0u,
        x1u=x1u,
        nphi=nphi,
    )

    # -----------------
    # func_pts_from_domain
    # -----------------

    func_ind_from_domain = _get_func_ind_from_domain(
        x0u=x0u,
        x1u=x1u,
        index=index,
        nphi=nphi,
        # debug
        func_RZphi_from_ind=func_RZphi_from_ind,
    )

    return func_RZphi_from_ind, func_ind_from_domain


# ##########################################################
# ##########################################################
#               checks
# ##########################################################


def _check(
    coll=None,
    key=None,
    res_phi=None,
):

    # ---------
    # mesh
    # ---------

    # key
    wm = coll._which_mesh
    key, _, cat = _generic_mesh._get_key_mesh_vs_bplines(
        coll=coll,
        key=key,
        which=wm,
    )

    lok = ['rect', 'tri']
    if coll.dobj[wm][key]['type'] not in lok:
        msg = (
            f"Only accepts mtype in {lok}\n"
            f"Provided: {key}\n"
        )
        raise Exception(msg)

    # ---------
    # res_phi
    # ---------

    res_phi = float(ds._generic_check._check_var(
        res_phi, 'res_phi',
        types=(int, float),
        sign='>0',
    ))

    return key, res_phi


# ##########################################################
# ##########################################################
#           pts_from_ind
# ##########################################################


def _get_func_RZphi_from_ind(
    x0u=None,
    x1u=None,
    nphi=None,
):

    # --------------
    # func
    # -------------

    def func(
        indr=None,
        indz=None,
        indphi=None,
        # resources
        x0u=x0u,
        x1u=x1u,
        nphi=nphi,
    ):

        # ----------
        # check
        # ----------

        if indphi is not None:
            c0 = (indr is not None) and (indr.shape == indphi.shape)
            if not c0:
                indr_str = 'None' if indr is None else indr.shape
                msg = (
                    "Arg indr and indphi must be the same shape!\n"
                    f"Provided:\n"
                    f"\t- indr.shape = {indr_str}\n"
                    f"\t- indphi.shape = {indphi.shape}\n"
                )
                raise Exception(msg)

        # ----------
        # R, Z
        # ----------

        # R
        if indr is not None:
            R = x0u[indr]
        else:
            R = None

        # Z
        if indz is not None:
            Z = x1u[indz]
        else:
            Z = None

        # dS
        dS = (x0u[1] - x0u[0]) * (x1u[1] - x1u[0])

        # -------------
        # phi
        # -------------

        dV = None
        phi = None
        if indr is not None:
            dV = np.full(indr.shape, np.nan)
            iru = np.unique(indr)

            if indphi is None:
                for iri in iru:
                    phii = _get_phi_from_nphi(nphi[iri])

                    iok = indr == iri
                    dV[iok] = dS * x0u[iri] * (phii[1] - phii[0])

            else:
                phi = np.full(indr.shape, np.nan)
                for iri in iru:
                    phii = _get_phi_from_nphi(nphi[iri])

                    iok = indr == iri
                    dV[iok] = dS * x0u[iri] * (phii[1] - phii[0])
                    phi[iok] = phii[indphi[iok]]

        return R, Z, phi, dV

    return func


# ##########################################################
# ##########################################################
#           pts_ind_from_domain
# ##########################################################


def _get_func_ind_from_domain(
    x0u=None,
    x1u=None,
    nphi=None,
    index=None,
    # debug
    func_RZphi_from_ind=None,
):

    x0f = np.repeat(x0u[:, None], x1u.size, axis=1)
    x1f = np.repeat(x1u[None, :], x0u.size, axis=0)
    ptsf = np.array([x0f.ravel(), x1f.ravel()]).T
    shapef = x0f.shape

    # ---------------------
    # get polygons - cross
    # ---------------------

    # ---------------------
    # get cross-section polygon with margin

    def func(
        # from domain
        DR=None,
        DZ=None,
        Dphi=None,
        # from poly
        pcross0=None,
        pcross1=None,
        phor0=None,
        phor1=None,
        # resources
        index=index,
        x0u=x0u,
        x1u=x1u,
        ptsf=ptsf,
        shapef=shapef,
        nphi=nphi,
        # debug
        debug=None,
        debug_msg=None,
        func_RZphi_from_ind=func_RZphi_from_ind,
    ):
        # ------------------
        # check
        # ------------------

        index = np.copy(index)

        DR, DZ, Dphi = _check_domain(
            DR=DR,
            DZ=DZ,
            Dphi=Dphi,
        )

        # ------------------
        # limits from domain
        # ------------------

        if DR is not None:
            index &= (x0u[:, None] >= DR[0]) & (x0u[:, None] <= DR[1])
        if DZ is not None:
            index &= (x1u[None, :] >= DZ[0]) & (x1u[None, :] <= DZ[1])

        # ------------------
        # limits from pcross
        # ------------------

        if pcross0 is not None:
            pcross = Path(np.array([pcross0, pcross1]).T)
            index &= pcross.contains_points(ptsf).reshape(shapef)

        # R and Z indices
        ir, iz = index.nonzero()
        iru = np.unique(ir)

        # ------------
        # indices
        # ------------

        ind = np.empty((3, 0), dtype=int)
        for ii, iri in enumerate(iru):

            # ---
            # iz

            izi = np.unique(iz[ir == iri])
            if izi.size == 0:
                print(f"SKIP: {ii} {iri}")
                continue

            # -----
            # iphi

            if phor0 is None and Dphi is None:
                iphi = np.arange(0, nphi[iri])

            elif Dphi is not None:
                phii = _get_phi_from_nphi(nphi[iri])
                if Dphi[0] < Dphi[1]:
                    iphi = np.nonzero(
                        (phii >= Dphi[0]) & (phii <= Dphi[1])
                    )[0]
                else:
                    iphi = np.nonzero(
                        (phii >= Dphi[0]) | (phii <= Dphi[1])
                    )[0]
            elif phor0 is not None:
                phii = np.pi*np.linspace(-1, 1, nphi[iri])
                pts = np.array([
                    x0u[iri]*np.cos(phii),
                    x0u[iri]*np.sin(phii),
                ]).T
                path = Path(np.array([phor0, phor1]).T)
                iphi = path.contains_points(pts).nonzero()[0]

            # ---------
            # group ind

            indi = np.array(
                [
                    np.full((iphi.size*izi.size,), iri),
                    np.repeat(izi, iphi.size),
                    np.tile(iphi, izi.size),
                ],
                dtype=int,
            )

            ind = np.concatenate((ind, indi), axis=1)

        # ----------
        # debug

        if debug is True:

            # coordinates
            rr, zz, pp, dV = func_RZphi_from_ind(ind)

            # title
            tit = (
                "func_ind_from_domain\n"
                f"pcross0 {pcross0 is not None} - phor0 {phor0 is not None}\n"
                f"DR {DR is not None} - DZ {DZ is not None}"
                f" - Dphi{Dphi is not None}"
            )
            if debug_msg is not None:
                tit += f"\n{debug_msg}"

            # figure and axes
            fig = plt.figure(figsize=(14, 8))
            fig.suptitle(tit, size=14, fontweight='bold')

            ax0 = fig.add_subplot(1, 2, 1, aspect='equal')
            ax0.set_xlabel("R (m)", size=12, fontweight='bold')
            ax0.set_ylabel("Z (m)", size=12, fontweight='bold')

            ax1 = fig.add_subplot(1, 2, 2, aspect='equal')
            ax1.set_xlabel("X (m)", size=12, fontweight='bold')
            ax1.set_ylabel("Y (m)", size=12, fontweight='bold')

            # domain
            if DR is not None:
                ax0.axvline(DR[0], c='k', ls='--')
                ax0.axvline(DR[1], c='k', ls='--')

            if DZ is not None:
                ax0.axhline(DZ[0], c='k', ls='--')
                ax0.axhline(DZ[1], c='k', ls='--')

            if Dphi is not None:
                maxr = np.max(rr)
                ax1.plot(
                    np.r_[0, maxr]*np.cos(Dphi[0]),
                    np.r_[0, maxr]*np.sin(Dphi[0]),
                    c='k',
                    ls='--',
                )
                ax1.plot(
                    np.r_[0, maxr]*np.cos(Dphi[1]),
                    np.r_[0, maxr]*np.sin(Dphi[1]),
                    c='k',
                    ls='--',
                )

            # polygons
            if pcross0 is not None:
                ax0.fill(pcross0, pcross1, fc=(0.5, 0.5, 0.5, 0.5))
            if phor0 is not None:
                ax1.fill(phor0, phor1, fc=(0.5, 0.5, 0.5, 0.5))

            # points
            ax0.plot(rr, zz, '.')
            ax1.plot(rr*np.cos(pp), rr*np.sin(pp), '.')

        return ind[0, :], ind[1, :], ind[2, :]

    return func


# ##########################################################
# ##########################################################
#           unique get phi function
# ##########################################################


def _get_phi_from_nphi(nphi):
    return np.pi * np.linspace(-1, 1, nphi, endpoint=False)


# ##########################################################
# ##########################################################
#           _check_domain
# ##########################################################


def _check_domain(
    DR=None,
    DZ=None,
    Dphi=None,
):
    # -----------
    #   DR
    # -----------

    if DR is not None:
        DR = ds._generic_check._check_flat1darray(
            DR, 'DR',
            size=2,
            dtype=float,
            unique=True,
            sign='>=0',
        )

    # -----------
    #   DZ
    # -----------

    if DZ is not None:
        DZ = ds._generic_check._check_flat1darray(
            DZ, 'DZ',
            size=2,
            dtype=float,
            unique=True,
        )

    # -----------
    #   Dphi
    # -----------

    if Dphi is not None:
        Dphi = ds._generic_check._check_flat1darray(
            Dphi, 'Dphi',
            size=2,
            dtype=float,
        )
        Dphi = [
            np.arctan2(np.sin(Dphi[0]), np.cos(Dphi[0])),
            np.arctan2(np.sin(Dphi[1]), np.cos(Dphi[1])),
        ]

    return DR, DZ, Dphi
