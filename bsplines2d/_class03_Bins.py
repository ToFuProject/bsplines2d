# -*- coding: utf-8 -*-


# Built-in
import copy


# Common
import numpy as np


# local
from ._class02_BSplines2D import BSplines2D as Previous
from . import _class03_checks as _checks
from . import _class03_binning as _binning


__all__ = ['Bins']


# #############################################################################
# #############################################################################
#
# #############################################################################


class Bins(Previous):

    _which_bins = 'bins'
    _ddef = copy.deepcopy(Previous._ddef)
    _dshow = dict(Previous._dshow)

    _dshow.update({
        _which_bins: [
            'nd',
            'cents',
            'shape',
            'ref',
        ],
    })

    # -----------------
    # bsplines
    # ------------------

    def add_bins(
        self,
        key=None,
        edges=None,
        # custom names
        key_ref=None,
        key_cents=None,
        key_res=None,
        # attributes
        **kwdargs,
    ):
        """ Add bin """

        # --------------
        # check inputs

        key, dref, ddata, dobj = _checks.check(
            coll=self,
            key=key,
            edges=edges,
            # custom names
            key_cents=key_cents,
            key_ref=key_ref,
            # attributes
            **kwdargs,
        )

        # --------------
        # update dict and crop if relevant

        self.update(dobj=dobj, ddata=ddata, dref=dref)

    # -----------------
    # binning tools
    # ------------------

    def binning(
        self,
        keys=None,
        ref_key=None,
        bins=None,
        # optional storing
        verb=None,
        store=None,
        returnas=None,
        key_store=None,
    ):
        """ Bin data along ref_key

        Binning is treated here as an integral
        Hence, if:
            - the data has units [ph/eV]
            - the ref_key has units [eV]
            - the binned data has units [ph]

        return a dict with data and units per key

        """

        return _binning.binning(
            coll=self,
            keys=keys,
            ref_key=ref_key,
            bins=bins,
            # optional storing
            verb=verb,
            returnas=returnas,
            store=store,
            key_store=key_store,
        )
