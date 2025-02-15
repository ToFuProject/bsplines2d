# -*- coding: utf-8 -*-


import itertools as itt


# Common
import numpy as np
import datastock as ds


# ###############################################################
# ###############################################################
#               Main
# ###############################################################


def _get_outline(
    coll=None,
    key=None,
    closed=None,
    plot_debug=None,
):

    # ------------
    # check inputs
    # ------------

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

    # plot_debug
    plot_debug = ds._generic_check._check_var(
        plot_debug, 'plot_debug',
        types=bool,
        default=False,
    )


    # ------------
    # compute
    # ------------

    if mtype == 'rect':
        i0, i1, knots0, knots1 = _rect_outline(
            coll=coll,
            key=key,
            plot_debug=plot_debug,
        )

        # close
        if closed is True:
            i0 = np.append(i0, i0[0])
            i1 = np.append(i1, i1[0])

        x0 = knots0[i0]
        x1 = knots1[i1]

        dout = {
            'i0': {
                'data': i0,
                'units': None,
                'dim': 'index',
            },
            'i1': {
                'data': i1,
                'units': None,
                'dim': 'index',
            },
        }

    else:
        ind, knots0, knots1 = _tri_outline(
            coll=coll,
            key=key,
            plot_debug=plot_debug,
        )

        # close
        if closed is True:
            ind = np.append(ind, ind[0])

        x0 = knots0[ind]
        x1 = knots1[ind]

        dout = {
            'ind': {
                'data': ind,
                'units': None,
                'dim': 'index',
            },
        }

    # -------------
    # format output
    # -------------

    k0, k1 = coll.dobj[wm][key]['knots']

    dout.update({
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
    })

    return dout


# ###############################################################
# ###############################################################
#               rectangular
# ###############################################################


def _rect_outline(coll=None, key=None, plot_debug=None):

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
        i0min, i0max = np.argmin(knots0), np.argmax(knots0)
        i1min, i1max = np.argmin(knots1), np.argmax(knots1)

        i0 = np.r_[i0min, i0max, i0max, i0min]
        i1 = np.r_[i1min, i1min, i1max, i1max]

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

        # find indices of all knots adjacent to at least one excluded cent
        ind0 = np.zeros(i0.shape, dtype=bool)
        linds = [range(ss) for ss in i0.shape]
        for ind in itt.product(*linds):
            sli = tuple(list(ind) + [slice(None)])
            iok = (ic0[sli] >= 0) & (ic1[sli] >= 0)
            sli = tuple(list(ind) + [iok])
            ind0[ind] = np.any(~crop[ic0[sli], ic1[sli]])

        # get subset of extreme coordinates
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
            if pp is None and plot_debug:
                import matplotlib.pyplot as plt
                msg = (
                    f"Isolated points detected in mesh '{key}':\n"
                    f"\n\t- i0: {i0}\n"
                    f"\t- i1: {i1}\n"
                    f"\t- ii0: {ii0}\n"
                    f"\t- old: {old}\n"
                    f"\t- pp: {pp}\n"
                    f"\t- pall: {pall}\n"
                )
                _ = coll.plot_mesh(key)
                plt.gca().plot(
                    knots0[np.array(lp)[:, 0]],
                    knots1[np.array(lp)[:, 1]],
                    ls='-',
                    marker='o',
                    lw=2,
                    label='optimized outline',
                )
                plt.gca().plot(
                    knots0[i0],
                    knots1[i1],
                    ls='None',
                    marker='x',
                    label='full outline',
                )
                plt.gca().plot(
                    knots0[np.array(pall)[:, 0]],
                    knots1[np.array(pall)[:, 1]],
                    ls='None',
                    marker='s',
                    label='isolated',
                )
                plt.gca().legend(loc='center left', bbox_to_anchor=(1., 0.5))
                raise Exception(msg)
            # -----------------------

            lp.append(pp)

        i0, i1 = np.array(lp).T

    return i0, i1, knots0, knots1


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


def _tri_outline(coll=None, key=None, plot_debug=None):

    # ----------
    # prepare
    # ----------

    wm = coll._which_mesh

    # knots
    k0, k1 = coll.dobj[wm][key]['knots']
    knots0 = coll.ddata[k0]['data']
    knots1 = coll.ddata[k1]['data']

    # indices
    ind = coll.ddata[coll.dobj[wm][key]['ind']]['data']

    # crop (not implemented for triangulare meshes)
    # crop = coll.dobj[wm][key]['crop']

    # -----------------------
    # compute
    # method: check if each segment belongs to 1 or 2 triangles
    # -----------------------

    # ---------------------
    # get all segments

    seg = np.array([
        np.sort(ind[:, :-1], axis=1),
        np.sort(ind[:, ::2], axis=1),
        np.sort(ind[:, 1:], axis=1),
    ])

    # ---------------------
    # list unique segments

    useg = np.unique(np.concatenate(seg, axis=0), axis=0)

    # ---------------------
    # nb of triangles

    ntri = np.array([
        np.any(np.all(seg == ss[None, None, :], axis=-1), axis=0).sum()
        for ss in useg
    ])

    # safety checks
    iout = (ntri != 1) & (ntri != 2)
    if np.any(iout):
        msg = (
            f"Outline of triangular mesh '{key}':\n"
            "Some segments in ntri are use not 1 or 2:\n"
            f"{iout.nonzero()[0]}"
        )
        raise Exception(msg)

    # --------------
    # seg_outline

    seg_out = useg[ntri == 1]

    # each node used exacty twice
    nodu, counts = np.unique(seg_out.ravel(), return_counts=True)
    iout = counts != 2
    if np.any(iout):
        lstr = [f"\t- {nodu[ii]}" for ii in iout.nonzero()[0]]
        msg = (
            f"Outline of triangular mesh '{key}':\n"
            "Some nodes in seg_out are not used exactly twice!\n"
            + "\n".join(lstr)
        )
        raise Exception(msg)

    # check nb of outline segment == number of unique nodes
    nnod = nodu.size
    if seg_out.shape[0] != nnod:
        msg = (
            f"Outline of triangular mesh '{key}':\n"
            "Nb. of outline segments does not match nb. of unique nodes\n"
            f"\t- nb outline segments: {seg_out.shape[0]}\n"
            f"\t- nb of unique nodes:  {nodu.size}\n"
        )
        raise Exception(msg)

    # --------------
    # order

    iout = np.r_[seg_out[0, :], np.zeros((nodu.size-2,), dtype=int)]
    iseg = np.arange(1, nnod)
    for ii in range(2, nnod):

        iok = seg_out[iseg, :] == iout[ii-1]
        iokn = np.any(iok, axis=1).nonzero()[0][0]
        isegi = iseg[iokn]

        iout[ii] = seg_out[isegi, :][~iok[iokn, :]][0]
        iseg = np.delete(iseg, iokn)

    return iout, knots0, knots1