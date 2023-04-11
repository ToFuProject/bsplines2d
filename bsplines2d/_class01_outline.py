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
        crop0, crop1 = (~crop).nonzero()
        crop = np.array([crop0, crop1]).T.tolist()

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

        ind0 = np.array([
            np.any([[ic0[ii, jj], ic1[ii, jj]] in crop for jj in [0, 1, 2, 3]])
            for ii in range(i0.size)
        ])
        ind = (
            ind0
            | (i0 == i0min) | (i0 == i0max)
            | (i1 == i1min) | (i1 == i1max)
        )
        if not np.any(ind):
            import pdb; pdb.set_trace()     # DB
            pass

        i00 = np.min(i0[ind])
        p0 = [i00, np.min(i1[ind][i0[ind] == i00])]
        pall = np.array([i0[ind], i1[ind]]).T.tolist()
        lp = [p0]

        old = None
        while len(lp) == 1 or lp[-1] != p0:
            pp, old = _next_pp(
                p0=lp[-1],
                pall=pall,
                old=old,
            )
            lp.append(pp)

        i0, i1 = np.array(lp).T
        x0 = knots0[i0]
        x1 = knots1[i1]

    return x0, x1


def _next_pp(p0=None, pall=None, old=None):

    # inc0, inc1

    p1 = np.copy(p0)
    found = False
    ldir = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    if old is not None:
        ldir = [pp for pp in ldir if pp != old]

    for (ic0, ic1) in ldir:

        stop = False
        while stop is False:

            p2 = [p1[0] + ic0, p1[1] + ic1]
            if p2 in pall:
                p1 = p2
                found = True
            else:
                stop = True

        if found is True:
            new = (-ic0, -ic1)
            break

    if found is False:
        raise Exception("Next point not found!")

    return p1, new


# ###############################################################
# ###############################################################
#               triangular
# ###############################################################


def _tri_outline(coll=None, key=None):

    msg = "outline not implemented yet for triangular meshes"
    raise NotImplementedError(msg)
