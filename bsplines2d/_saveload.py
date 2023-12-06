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

    lkbs = list(coll.dobj.get(wbs, {}).keys())
    for kbs in lkbs:
        coll._dobj[wbs][kbs]['class'] = None



# ############################################################
# ############################################################
#                 BSplines - loading
# ############################################################


def load(
    pfe=None,
    cls=None,
    allow_pickle=None,
    sep=None,
    verb=None,
):

    # --------------------
    # use datastock.load()

    from ._class02_BSplines2D import BSplines2D

    coll = ds.load(
        pfe=pfe,
        cls=BSplines2D,
        allow_pickle=allow_pickle,
        sep=sep,
        verb=verb,
    )

    # ----------------
    # re-build classes

    wbs = coll._which_bsplines
    lkbs = list(coll.dobj.get(wbs, {}).keys())
    for kbs in lkbs:
        pass

    return coll