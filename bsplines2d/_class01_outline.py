# -*- coding: utf-8 -*-


# Common
import numpy as np
import datastock as ds


# ###############################################################
# ###############################################################
#               Main
# ###############################################################


def _get_outline(coll=None, key=None, closed=None):

    # ------------
    # check inputs

    # key
    wm = coll._which_mesh
    lok = [
        k0 for k0, v0 in coll.dobj.get(wm, {}).items()
        if v0['nd'] == '2d'
    ]
    key = ds._generic_check._check_var(
        key, 'key',
        types=str,
        allowed=lok,
    )

    mtype = coll.dobj[wm][key]['type']

    # closed
    closed = ds._generic_check._check_var(
        closed, 'closed',
        types=bool,
        default=True,
    )

    # ------------
    # compute

    if mtype == 'rect':
        x0, x1 = _rect_outline(
            coll=coll,
            key=key,
        )

    else:
        x0, x1 = _tri_outline(
            coll=coll,
            key=key,
        )

    # close
    if closed is True:
        x0 = np.append(x0, x0[0])
        x1 = np.append(x1, x1[0])

    # -------------
    # format output

    k0, k1 = coll.dobj[wm][key]['knots']

    dout = {
        'x0': {
            'data': x0,
            'units': coll.ddata[k0]['units'],
            'dim': coll.ddata[k0]['dim'],
            'name': coll.ddata[k0]['name'],
            'quant': coll.ddata[k0]['quant'],
        },
        'x1': {
            'data': x1,
            'units': coll.ddata[k1]['units'],
            'dim': coll.ddata[k1]['dim'],
            'name': coll.ddata[k1]['name'],
            'quant': coll.ddata[k1]['quant'],
        },
    }

    return dout


# ###############################################################
# ###############################################################
#               rectangular
# ###############################################################


def _rect_outline(coll=None, key=None):

    # ----------
    # prepare

    wm = coll._which_mesh

    # knots
    k0, k1 = coll.dobj[wm][key]['knots']
    knots0 = coll.ddata[k0]['data']
    knots1 = coll.ddata[k1]['data']

    # crop
    crop = coll.dobj[wm][key]['crop']

    # -------
    # compute

    if crop is False:
        x0min, x0max = knots0.min(), knots0.max()
        x1min, x1max = knots1.min(), knots1.max()
        x0 = np.r_[x0min, x0max, x0max, x0min]
        x1 = np.r_[x1min, x1min, x1max, x1max]

    else:
        crop = coll.ddata[crop]['data']

        # indices of knots
        (i0, i1), (ic0, ic1) = coll.select_mesh_elements(
            key=key,
            crop=True,
            elements='knots',
            returnas='ind',
            return_neighbours=True,
        )

        # only keep edge knots
        i0min, i1min = i0.min(), i1.min()
        i0max, i1max = i0.max(), i1.max()

        # get indices of edge knots
        iok = (ic0 >= 0) & (ic1 >= 0)
        ind0 = np.array([
            np.any(~crop[ic0[ii, iok[ii, :]], ic1[ii, iok[ii, :]]])
            for ii in range(i0.size)
        ])

        # get subest of extreme coordinates
        ind = (
            ind0
            | (i0 == i0min) | (i0 == i0max)
            | (i1 == i1min) | (i1 == i1max)
        )
        if not np.any(ind):
            import pdb; pdb.set_trace()     # DB
            pass

        # keep only limits
        i0 = i0[ind]
        i1 = i1[ind]
        ic0 = ic0[ind, :]
        ic1 = ic1[ind, :]

        # get starting point
        i0_min = np.min(i0)
        i1_min = np.min(i1[i0 == i0_min])
        ii0 = ((i0 == i0_min) & (i1 == i1_min)).nonzero()[0][0]

        # pall = list(range(0, i0.size))

        p0 = [i0[ii0], i1[ii0]]
        pall = np.array([i0, i1]).T.tolist()

        lp = [p0]
        pall.remove(p0)

        old = None
        # while len(lp) == 1 or lp[-1] != p0:
        while len(pall) > 0:
            pp, old, ii0 = _next_pp(
                p0=lp[-1],
                ii0=ii0,
                pall=pall,
                i0=i0,
                i1=i1,
                ic0=ic0,
                ic1=ic1,
                crop=crop,
                old=old,
            )

            # ------ DEBUG --------
            if pp is None:
                import matplotlib.pyplot as plt
                msg = (
                    f"\n\t- i0: {i0}\n"
                    f"\t- i1: {i1}\n"
                    f"\t- ii0: {ii0}\n"
                    f"\t- old: {old}\n"
                    f"\t- pp: {pp}\n"
                    f"\t- pall: {pall}\n"
                )
                print(msg)
                dax = coll.plot_mesh(key)
                plt.gca().plot(
                    knots0[np.array(lp)[:, 0]],
                    knots1[np.array(lp)[:, 1]],
                    ls='-',
                    marker='o',
                    lw=2,
                )
                plt.gca().plot(
                    knots0[i0],
                    knots1[i1],
                    ls='None',
                    marker='x',
                )
                plt.gca().plot(
                    knots0[np.array(pall)[:, 0]],
                    knots1[np.array(pall)[:, 1]],
                    ls='None',
                    marker='s',
                )

                break
            # -----------------------

            lp.append(pp)

        i0, i1 = np.array(lp).T
        x0 = knots0[i0]
        x1 = knots1[i1]

    return x0, x1


def _next_pp(
    p0=None,
    ii0=None,
    pall=None,
    i0=None,
    i1=None,
    ic0=None,
    ic1=None,
    crop=None,
    old=None,
):

    # inc0, inc1
    p1 = np.copy(p0)
    found = False
    ldir = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    if old is not None:
        ldir = [pp for pp in ldir if pp != old]

    for (ip0, ip1) in ldir:

        stop = False
        while stop is False:

            p2 = [p1[0] + ip0, p1[1] + ip1]

            if p2 in pall:

                ii2 = ((i0 == p2[0]) & (i1 == p2[1])).nonzero()[0][0]
                icc0 = np.intersect1d(
                    ic0[ii0, :],
                    ic0[ii2, :],
                    assume_unique=False,
                    return_indices=False,
                )
                icc1 = np.intersect1d(
                    ic1[ii0, :],
                    ic1[ii2, :],
                    assume_unique=False,
                    return_indices=False,
                )

                i00 = np.array([
                    (ic0[ii0, :] == cc0).nonzero()[0] for cc0 in icc0
                ]).ravel()
                i11 = np.array([
                    (ic1[ii0, :] == cc1).nonzero()[0] for cc1 in icc1
                ]).ravel()

                ic = np.intersect1d(
                    i00,
                    i11,
                    assume_unique=False,
                    return_indices=False,
                )

                c0 = np.any(crop[ic0[ii0, ic], ic1[ii0, ic]])

                if c0:
                    p1 = p2
                    ii0 = ii2
                    found = True
                    pall.remove(p2)

                else:
                    stop = True

            else:
                stop = True

        if found is True:
            new = (-ip0, -ip1)
            break

    if found is False:
        return None, None, None

    return p1, new, ii0


# ###############################################################
# ###############################################################
#               triangular
# ###############################################################


def _tri_outline(coll=None, key=None):

    msg = "outline not implemented yet for triangular meshes"
    raise NotImplementedError(msg)