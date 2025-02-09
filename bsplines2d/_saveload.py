# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 08:25:06 2023

@author: dvezinet
"""


import datastock as ds


# ############################################################
# ############################################################
#                 BSplines - saving
# ############################################################


def prepare_bsplines(coll=None):

    # -----------------
    # Remove classes

    wbs = coll._which_bsplines

    dclas = {}
    lkbs = list(coll.dobj.get(wbs, {}).keys())
    for kbs in lkbs:
        dclas[kbs] = coll._dobj[wbs][kbs]['class']
        coll._dobj[wbs][kbs]['class'] = None

    return dclas


def restore_bsplines(coll=None, dclas=None):

    # -----------------
    # Remove classes

    wbs = coll._which_bsplines

    lkbs = list(coll.dobj.get(wbs, {}).keys())
    for kbs in lkbs:
        coll._dobj[wbs][kbs]['class'] = dclas[kbs]

    return


# ############################################################
# ############################################################
#                 BSplines - loading
# ############################################################


def load(
    pfe=None,
    cls=None,
    coll=None,
    allow_pickle=None,
    sep=None,
    verb=None,
):

    # --------------------
    # use datastock.load()

    from . import _class02_compute as _compute

    if cls is None:
        from ._class03_Bins import Bins
        cls = Bins

    coll = ds.load(
        pfe=pfe,
        cls=cls,
        coll=coll,
        allow_pickle=allow_pickle,
        sep=sep,
        verb=verb,
    )

    # ----------------
    # re-build classes

    wm = coll._which_mesh
    wbs = coll._which_bsplines
    lkbs = list(coll.dobj.get(wbs, {}).keys())
    for kbs in lkbs:

        keym = coll.dobj[wbs][kbs][wm]
        nd = coll.dobj[wm][keym]['nd']
        mtype = coll.dobj[wm][keym]['type']
        deg = coll.dobj[wbs][kbs]['deg']
        if nd == '1d':
            dref, ddata, dobj = _compute._mesh1d_bsplines(
                coll=coll, keym=keym, keybs=kbs, deg=deg,
            )
        elif mtype == 'rect':
            dref, ddata, dobj = _compute._mesh2DRect_bsplines(
                coll=coll, keym=keym, keybs=kbs, deg=deg,
            )
        elif mtype == 'tri':
            dref, ddata, dobj = _compute._mesh2DTri_bsplines(
                coll=coll, keym=keym, keybs=kbs, deg=deg,
            )
        else:
            dref, ddata, dobj = _compute._mesh2Dpolar_bsplines(
                coll=coll, keym=keym, keybs=kbs, deg=deg, # angle=angle,
            )

        coll._dobj[wbs][kbs]['class'] = dobj[wbs][kbs]['class']

    return coll