

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datastock as ds


# ###########################################
# ###########################################
#               Main
# ###########################################


def main(
    coll=None,
    key=None,
    operator=None,
    geometry=None,
    crop=None,
    # store vs return
    return_param=None,
    # specific to deg = 0
    centered=None,
    # to return gradR, gradZ, for D1N2 deg 0, for tomotok
    returnas_element=None,
    # plotting
    dax=None,
    fs=None,
    dmargin=None,
    fontsize=None,
):

    # ------------
    # check inputs
    # ------------

    key, operator = _check(**locals())

    # ------------
    # get operator
    # ------------

    dout = coll.add_bsplines_operator(
        key=key,
        operator=operator,
        geometry=geometry,
        crop=crop,
        store=False,
        returnas=True,
        return_param=return_param,
        # specific to deg = 0
        centered=centered,
        # to return gradR, gradZ, for D1N2 deg 0, for tomotok
        returnas_element=returnas_element,
    )

    # ------------
    # prepare axes
    # ------------

    if dax is None:
        dax = _get_dax(
            coll=coll,
            key=key,
            operator=operator,
            dmargin=dmargin,
            fs=fs,
            fontsize=fontsize,
        )

    dax = ds._generic_check._check_dax(dax)

    # ------------
    # plot
    # ------------

    if operator == 'D1N2':

        _plot_D1N2(dax=dax, dout=dout)

    return dax


# ###########################################
# ###########################################
#               Check
# ###########################################


def _check(
    coll=None,
    key=None,
    operator=None,
    # unused
    **kwdargs,
):

    # -------------
    # key
    # -------------

    wbs = coll._which_bsplines
    lok = list(coll.dobj.get(wbs, {}).keys())
    key = ds._generic_check._check_var(
        key, 'key',
        types=str,
        allowed=lok,
    )

    # -------------
    # operator
    # -------------

    lok = ['D1N2']
    operator = ds._generic_check._check_var(
        operator, 'operator',
        types=str,
        allowed=lok,
        extra_msg="Not implemented for operator '{operator}'!\n"
    )

    return key, operator


# ###########################################
# ###########################################
#               get_dax
# ###########################################


def _get_dax(
    coll=None,
    key=None,
    operator=None,
    fs=None,
    fontsize=None,
    dmargin=None,
):

    # ------------
    # check inputs
    # ------------

    # dmargin
    if dmargin is None:
        dmargin = {
            'left': 0.08, 'right': 0.95,
            'bottom': 0.08, 'top': 0.9,
            'hspace': 0.1, 'wspace': 0.2,
        }

    # fontsize
    fontsize = float(ds._generic_check._check_var(
        fontsize, 'fontsize',
        types=(int, float),
        sign='>0',
        default=14,
    ))

    if fs is None:
        fs = (15, 10)

    # ------------
    # figure
    # ------------

    fig = plt.figure(figsize=fs)

    wbs = coll._which_bsplines
    tit = (
        f"Operator '{operator}' for {wbs}'{key}'"
    )
    fig.suptitle(tit, size=fontsize + 2, fontweight='bold')

    dax = {}

    # ------------
    # D1N2
    # ------------

    if operator == 'D1N2':

        gs = gridspec.GridSpec(ncols=2, nrows=1, **dmargin)

        # gradR
        ax0 = fig.add_subplot(gs[0, 0])
        ax0.set_title(
            r'$tMM0 = \nabla^T_R\nabla_R$',
            fontsize=fontsize,
            fontweight='bold',
        )
        ax0.set_xlabel(
            'bspline index',
            fontsize=fontsize,
            fontweight='bold',
        )
        ax0.set_ylabel(
            'bspline index',
            fontsize=fontsize,
            fontweight='bold',
        )

        dax['gradR'] = {'handle': ax0}

        # gradZ
        ax = fig.add_subplot(gs[0, 1], sharex=ax0, sharey=ax0)
        ax.set_title(
            r'$tMM1 = \nabla^T_Z\nabla_Z$',
            fontsize=fontsize,
            fontweight='bold',
        )
        ax.set_xlabel(
            'bspline index',
            fontsize=fontsize,
            fontweight='bold',
        )
        ax.set_ylabel(
            'bspline index',
            fontsize=fontsize,
            fontweight='bold',
        )

        dax['gradZ'] = {'handle': ax}

    return dax


# ###########################################
# ###########################################
#               D1N2
# ###########################################


def _plot_D1N2(
    dax=None,
    dout=None,
):

    # ------------
    # gradR
    # ------------

    for kax, kdata in [('gradR', 'tMM0'), ('gradZ', 'tMM1')]:

        if dax.get(kax) is not None:
            ax = dax[kax]['handle']

            # max value
            vmax = np.max(np.abs(dout[kdata]['data']))

            # imshow
            ax.imshow(
                dout[kdata]['data'].todense(),
                origin='upper',
                interpolation='nearest',
                cmap=plt.cm.seismic,
                vmax=vmax,
                vmin=-vmax,
            )

    return
