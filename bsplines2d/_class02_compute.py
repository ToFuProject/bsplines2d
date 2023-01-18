# -*- coding: utf-8 -*-


# Built-in
import warnings

# Common
import numpy as np
from scipy.spatial import ConvexHull
from matplotlib.path import Path
from contourpy import contour_generator
import datastock as ds


# tofu
from . import _generic_mesh
from . import _utils_bsplines
# from . import _class02_checks as _checks
from . import _class02_bsplines_rect
from . import _class02_bsplines_tri
from . import _class02_bsplines_polar
from . import _class02_bsplines_1d


# #############################################################################
# #############################################################################
#                           Mesh2 - Tri - bsplines
# #############################################################################


def _mesh2DTri_bsplines(coll=None, keym=None, keybs=None, deg=None):

    # --------------
    # create bsplines

    kknots = coll.dobj[coll._which_mesh][keym]['knots']
    func_details, func_sum, clas = _class1_bsplines_tri.get_bs2d_func(
        deg=deg,
        knotsR=coll.ddata[kknots[0]]['data'],
        knotsZ=coll.ddata[kknots[1]]['data'],
        cents=coll.ddata[coll.dobj[coll._which_mesh][keym]['ind']]['data'],
        trifind=coll.dobj[coll._which_mesh][keym]['func_trifind'],
    )
    keybsr = f'{keybs}-nbs'
    kbscr = f'{keybs}-apR'
    kbscz = f'{keybs}-apZ'

    bs_cents = clas._get_bs_cents()

    # ----------------
    # format into dict

    dref = {
        # bs index
        keybsr: {
            'size': clas.nbs,
        },
    }

    ddata = {
        kbscr: {
            'data': bs_cents[0, :],
            'units': 'm',
            'dim': 'distance',
            'quant': 'R',
            'name': 'R',
            'ref': (keybsr,),
        },
        kbscz: {
            'data': bs_cents[1, :],
            'units': 'm',
            'dim': 'distance',
            'quant': 'Z',
            'name': 'Z',
            'ref': (keybsr,),
        },
    }

    dobj = {
        'bsplines': {
            keybs: {
                'deg': deg,
                'mesh': keym,
                'ref': (keybsr,),
                'ref-bs': (keybsr,),
                'apex': (kbscr, kbscz),
                'shape': (clas.nbs,),
                'crop': False,
                'func_details': func_details,
                'func_sum': func_sum,
                'class': clas,
            }
        },
    }

    return dref, ddata, dobj


# #############################################################################
# #############################################################################
#                           Mesh1D - bsplines
# #############################################################################


def _mesh1d_bsplines(
    coll=None,
    keym=None,
    keybs=None,
    deg=None,
    which_mesh=None,
    which_bsplines='bsplines',
):
    
    kknots = coll.dobj[which_mesh][keym]['knots'][0]
    knots = coll.ddata[kknots]['data']

    kbsn = f'{keybs}-nbs'
    kbsap = f'{keybs}-ap'

    func_details, func_sum, clas = _class1_bsplines_1d.get_bs2d_func(
        deg=deg,
        knots=knots,
        coll=coll,
    )

    # ------------
    # refs

    ref = (kbsn,)
    apex = (kbsap,)

    # ----------------
    # format into dict

    # dref
    dref = {
        # bs index
        kbsn: {'size': clas.nbs},
    }

    # ddata
    ddata = {
        kbsap: {
            'data': clas.apex_per_bs,
            'units': '',
            'dim': '',
            'quant': '',
            'name': '',
            'ref': (kbsn,),
        },
    }

    # dobj
    dobj = {
        which_bsplines: {
            keybs: {
                'deg': deg,
                which_mesh: keym,
                'ref': ref,
                'ref-bs': (kbsn,),
                'apex': apex,
                'shape': clas.shapebs,
                'func_details': func_details,
                'func_sum': func_sum,
                'class': clas,
                'crop': None,
            }
        },
    }

    return dref, ddata, dobj


# #############################################################################
# #############################################################################
#                           Mesh2DRect - bsplines
# #############################################################################


def _mesh2DRect_bsplines(coll=None, keym=None, keybs=None, deg=None):

    # --------------
    # create bsplines

    kR, kZ = coll.dobj[coll._which_mesh][keym]['knots']
    Rknots = coll.ddata[kR]['data']
    Zknots = coll.ddata[kZ]['data']

    keybsr = f'{keybs}-nbs'
    kRbsapn = f'{keybs}-nR'
    kZbsapn = f'{keybs}-nZ'
    kRbsap = f'{keybs}-apR'
    kZbsap = f'{keybs}-apZ'

    (
        shapebs, Rbs_apex, Zbs_apex,
        knots_per_bs_R, knots_per_bs_Z,
    ) = _class1_bsplines_rect.get_bs2d_RZ(
        deg=deg, Rknots=Rknots, Zknots=Zknots,
    )
    nbs = int(np.prod(shapebs))

    func_details, func_sum, clas = _class1_bsplines_rect.get_bs2d_func(
        deg=deg,
        Rknots=Rknots,
        Zknots=Zknots,
        shapebs=shapebs,
        # knots_per_bs_R=knots_per_bs_R,
        # knots_per_bs_Z=knots_per_bs_Z,
    )

    # ----------------
    # format into dict

    dref = {
        kRbsapn: {
            'size': Rbs_apex.size,
        },
        kZbsapn: {
            'size': Zbs_apex.size,
        },
        keybsr: {
            'size': nbs,
        },
    }

    ddata = {
        kRbsap: {
            'data': Rbs_apex,
            'units': 'm',
            'dim': 'distance',
            'quant': 'R',
            'name': 'R',
            'ref': kRbsapn,
        },
        kZbsap: {
            'data': Zbs_apex,
            'units': 'm',
            'dim': 'distance',
            'quant': 'Z',
            'name': 'Z',
            'ref': kZbsapn,
        },
    }

    dobj = {
        'bsplines': {
            keybs: {
                'deg': deg,
                'mesh': keym,
                'ref': (kRbsapn, kZbsapn),
                'ref-bs': (keybsr,),
                'apex': (kRbsap, kZbsap),
                'shape': shapebs,
                'crop': False,
                'func_details': func_details,
                'func_sum': func_sum,
                'class': clas,
            }
        },
    }

    return dref, ddata, dobj


def add_cropbs_from_crop(coll=None, keybs=None, keym=None):

    # ----------------
    # get

    kcropbs = False
    if coll.dobj[coll._which_mesh][keym]['crop'] is not False:
        kcropm = coll.dobj[coll._which_mesh][keym]['crop']
        cropbs = _get_cropbs_from_crop(
            coll=coll,
            crop=coll.ddata[kcropm]['data'],
            keybs=keybs,
        )
        kcropbs = f'{keybs}-crop'
        kcroppedbs = f'{keybs}-nbs-crop'

    # ----------------
    # optional crop

    if kcropbs is not False:

        # add cropped flat reference
        coll.add_ref(
            key=kcroppedbs,
            size=int(cropbs.sum()),
        )

        coll.add_data(
            key=kcropbs,
            data=cropbs,
            ref=coll._dobj['bsplines'][keybs]['ref'],
            dim='bool',
            quant='bool',
        )
        coll._dobj['bsplines'][keybs]['crop'] = kcropbs


def _mesh2DRect_bsplines_knotscents(
    returnas=None,
    return_knots=None,
    return_cents=None,
    ind=None,
    deg=None,
    Rknots=None,
    Zknots=None,
    Rcents=None,
    Zcents=None,
):

    # -------------
    # check inputs

    return_knots = ds._generic_check._check_var(
        return_knots, 'return_knots',
        types=bool,
        default=True,
    )
    return_cents = ds._generic_check._check_var(
        return_cents, 'return_cents',
        types=bool,
        default=True,
    )
    if return_knots is False and return_cents is False:
        return

    # -------------
    # compute

    if return_knots is True:

        knots_per_bs_R = _utils_bsplines._get_knots_per_bs(
            Rknots, deg=deg, returnas=returnas,
        )
        knots_per_bs_Z = _utils_bsplines._get_knots_per_bs(
            Zknots, deg=deg, returnas=returnas,
        )
        if ind is not None:
            knots_per_bs_R = knots_per_bs_R[:, ind[0]]
            knots_per_bs_Z = knots_per_bs_Z[:, ind[1]]

        nknots = knots_per_bs_R.shape[0]
        knots_per_bs_R = np.tile(knots_per_bs_R, (nknots, 1))
        knots_per_bs_Z = np.repeat(knots_per_bs_Z, nknots, axis=0)

    if return_cents is True:

        cents_per_bs_R = _utils_bsplines._get_cents_per_bs(
            Rcents, deg=deg, returnas=returnas,
        )
        cents_per_bs_Z = _utils_bsplines._get_cents_per_bs(
            Zcents, deg=deg, returnas=returnas,
        )
        if ind is not None:
            cents_per_bs_R = cents_per_bs_R[:, ind[0]]
            cents_per_bs_Z = cents_per_bs_Z[:, ind[1]]

        ncents = cents_per_bs_R.shape[0]
        cents_per_bs_R = np.tile(cents_per_bs_R, (ncents, 1))
        cents_per_bs_Z = np.repeat(cents_per_bs_Z, ncents, axis=0)

    # -------------
    # return

    if return_knots is True and return_cents is True:
        out = (
            (knots_per_bs_R, knots_per_bs_Z), (cents_per_bs_R, cents_per_bs_Z)
        )
    elif return_knots is True:
        out = (knots_per_bs_R, knots_per_bs_Z)
    else:
        out = (cents_per_bs_R, cents_per_bs_Z)
    return out


# #############################################################################
# #############################################################################
#                           Mesh2D - polar - bsplines
# #############################################################################


def _mesh2Dpolar_bsplines(
    coll=None,
    keym=None,
    keybs=None,
    angle=None,
    deg=None,
):

    # ---------------
    # create bsplines

    kknots = coll.dobj[coll._which_mesh][keym]['knots']
    knotsr = coll.ddata[kknots[0]]['data']
    if len(kknots) == 2:
        angle = coll.ddata[kknots[1]]['data']

    func_details, func_sum, clas = _class1_bsplines_polar.get_bs2d_func(
        deg=deg,
        knotsr=knotsr,
        angle=angle,
        coll=coll,
    )

    keybsnr = f'{keybs}-nr'
    keybsn = f'{keybs}-nbs'
    keybsapr = f'{keybs}-apr'

    # ------------
    # refs

    if clas.knotsa is None:
        ref = (keybsnr,)
        apex = (keybsapr,)
    elif len(clas.shapebs) == 2:
        keybsna = f'{keybs}-na'
        keybsapa = f'{keybs}-apa'
        ref = (keybsnr, keybsna)
        apex = (keybsapr, keybsapa)
    else:
        ref = (keybsn,)
        apex = (keybsapr,)

        # check angle vs angle2d
        mesh = coll._which_mesh
        angle2d = coll.dobj[mesh][keym]['angle2d']
        if angle2d is None:
            msg = (
                "Poloidal bsplines require mesh with angle2d!\n"
                f"\t- self.dobj['{mesh}']['{keym}']['angle2d'] = {angle2d}"
            )
            raise Exception(msg)

    # bs_cents = clas._get_bs_cents()

    # ----------------
    # format into dict

    # dref
    dref = {
        # bs index
        keybsnr: {'size': clas.nbs_r},
        keybsn: {'size': clas.nbs},
    }
    if len(clas.shapebs) == 2:
        dref[keybsna] = {'size': clas.nbs_a_per_r[0]}

    # ddata
    ddata = {
        keybsapr: {
            'data': clas.apex_per_bs_r,
            'units': '',
            'dim': '',
            'quant': '',
            'name': '',
            'ref': (keybsnr,),
        },
    }
    if len(clas.shapebs) == 2:
        ddata[keybsapa] = {
            'data': clas.apex_per_bs_a[0],
            'units': 'rad',
            'dim': 'angle',
            'quant': '',
            'name': '',
            'ref': (keybsna,),
        }

    # dobj
    dobj = {
        'bsplines': {
            keybs: {
                'deg': deg,
                'mesh': keym,
                'ref': ref,
                'ref-bs': (keybsn,),
                'apex': apex,
                'shape': clas.shapebs,
                'func_details': func_details,
                'func_sum': func_sum,
                'class': clas,
                'crop': coll.dobj[coll._which_mesh][keym]['crop'],
            }
        },
    }

    return dref, ddata, dobj


def _mesh2DPolar_bsplines_knotscents(
    returnas=None,
    return_knots=None,
    return_cents=None,
    ind=None,
    deg=None,
    # resources
    clas=None,
    rknots=None,
    aknots=None,
    rcents=None,
    acents=None,
):

    # -------------
    # check inputs

    return_knots = ds._generic_check._check_var(
        return_knots, 'return_knots',
        types=bool,
        default=True,
    )
    return_cents = ds._generic_check._check_var(
        return_cents, 'return_cents',
        types=bool,
        default=True,
    )
    if return_knots is False and return_cents is False:
        return

    # -------------
    # compute

    if return_knots is True:

        knots_per_bs_r = _utils_bsplines._get_knots_per_bs(
            rknots, deg=deg, returnas=returnas,
        )
        knots_per_bs_Z = _utils_bsplines._get_knots_per_bs(
            Zknots, deg=deg, returnas=returnas,
        )
        if ind is not None:
            knots_per_bs_R = knots_per_bs_R[:, ind[0]]
            knots_per_bs_Z = knots_per_bs_Z[:, ind[1]]

        nknots = knots_per_bs_R.shape[0]
        knots_per_bs_R = np.tile(knots_per_bs_R, (nknots, 1))
        knots_per_bs_Z = np.repeat(knots_per_bs_Z, nknots, axis=0)

    if return_cents is True:

        cents_per_bs_R = _utils_bsplines._get_cents_per_bs(
            Rcents, deg=deg, returnas=returnas,
        )
        cents_per_bs_Z = _utils_bsplines._get_cents_per_bs(
            Zcents, deg=deg, returnas=returnas,
        )
        if ind is not None:
            cents_per_bs_R = cents_per_bs_R[:, ind[0]]
            cents_per_bs_Z = cents_per_bs_Z[:, ind[1]]

        ncents = cents_per_bs_R.shape[0]
        cents_per_bs_R = np.tile(cents_per_bs_R, (ncents, 1))
        cents_per_bs_Z = np.repeat(cents_per_bs_Z, ncents, axis=0)

    # -------------
    # return

    if return_knots is True and return_cents is True:
        out = (
            (knots_per_bs_R, knots_per_bs_Z), (cents_per_bs_R, cents_per_bs_Z)
        )
    elif return_knots is True:
        out = (knots_per_bs_R, knots_per_bs_Z)
    else:
        out = (cents_per_bs_R, cents_per_bs_Z)
    return out


# #############################################################################
# #############################################################################
#                           Mesh2DRect - crop
# #############################################################################


def _get_cropbs_from_crop(coll=None, crop=None, keybs=None):

    if isinstance(crop, str) and crop in coll.ddata.keys():
        crop = coll.ddata[crop]['data']

    shref = coll.dobj[coll._which_mesh][coll.dobj['bsplines'][keybs]['mesh']]['shape-c']
    if crop.shape != shref:
        msg = "Arg crop seems to have the wrong shape!"
        raise Exception(msg)

    keym = coll.dobj['bsplines'][keybs][coll._which_mesh]
    kRk, kZk = coll.dobj['mesh'][keym]['knots']
    kRc, kZc = coll.dobj['mesh'][keym]['cents']

    cents_per_bs_R, cents_per_bs_Z = _mesh2DRect_bsplines_knotscents(
        returnas='ind',
        return_knots=False,
        return_cents=True,
        ind=None,
        deg=coll.dobj['bsplines'][keybs]['deg'],
        Rknots=coll.ddata[kRk]['data'],
        Zknots=coll.ddata[kZk]['data'],
        Rcents=coll.ddata[kRc]['data'],
        Zcents=coll.ddata[kZc]['data'],
    )

    shapebs = coll.dobj['bsplines'][keybs]['shape']
    cropbs = np.array([
        [
            np.all(crop[cents_per_bs_R[:, ii], cents_per_bs_Z[:, jj]])
            for jj in range(shapebs[1])
        ]
        for ii in range(shapebs[0])
    ], dtype=bool)

    return cropbs


# #############################################################################
# #############################################################################
#                           Mesh2DRect - interp utility
# #############################################################################


def _get_keyingroup_ddata(
    dd=None, dd_name='data',
    key=None, monot=None,
    msgstr=None, raise_=False,
):
    """ Return the unique data key matching key

    Here, key can be interpreted as name / source / units / quant...
    All are tested using select() and a unique match is returned
    If not unique match an error message is either returned or raised

    """

    # ------------------------
    # Trivial case: key is actually a ddata key

    if key in dd.keys():
        return key, None

    # ------------------------
    # Non-trivial: check for a unique match on other params

    dind = _select(
        dd=dd, dd_name=dd_name,
        dim=key, quant=key, name=key, units=key, source=key,
        monot=monot,
        log='raw',
        returnas=bool,
    )
    ind = np.array([ind for kk, ind in dind.items()])

    # Any perfect match ?
    nind = np.sum(ind, axis=1)
    sol = (nind == 1).nonzero()[0]
    key_out, msg = None, None
    if sol.size > 0:
        if np.unique(sol).size == 1:
            indkey = ind[sol[0], :].nonzero()[0]
            key_out = list(dd.keys())[indkey]
        else:
            lstr = "[dim, quant, name, units, source]"
            msg = "Several possible matches in {} for {}".format(lstr, key)
    else:
        lstr = "[dim, quant, name, units, source]"
        msg = "No match in {} for {}".format(lstr, key)

    # Complement error msg and optionally raise
    if msg is not None:
        lk = ['dim', 'quant', 'name', 'units', 'source']
        dk = {
            kk: (
                dind[kk].sum(),
                sorted(set([vv[kk] for vv in dd.values()]))
            ) for kk in lk
        }
        msg += (
            "\n\nRequested {} could not be identified!\n".format(msgstr)
            + "Please provide a valid (unique) key/name/dim/quant/units:\n\n"
            + '\n'.join([
                '\t- {} ({} matches): {}'.format(kk, dk[kk][0], dk[kk][1])
                for kk in lk
            ])
            + "\nProvided:\n\t'{}'".format(key)
        )
        if raise_:
            raise Exception(msg)
    return key_out, msg


def _get_possible_ref12d(
    dd=None,
    key=None, ref1d=None, ref2d=None,
    group1d='radius',
    group2d='mesh2d',
):

    # Get relevant lists
    kq, msg = _get_keyingroup_ddata(
        dd=dd,
        key=key, group=group2d, msgstr='quant', raise_=False,
    )

    if kq is not None:
        # The desired quantity is already 2d
        k1d, k2d = None, None

    else:
        # Check if the desired quantity is 1d
        kq, msg = _get_keyingroup_ddata(
            dd=dd,
            key=key, group=group1d,
            msgstr='quant', raise_=True,
        )

        # Get dict of possible {ref1d: lref2d}
        ref = [rr for rr in dd[kq]['ref'] if dd[rr]['group'] == (group1d,)][0]
        lref1d = [
            k0 for k0, v0 in dd.items()
            if ref in v0['ref'] and v0['monot'][v0['ref'].index(ref)] is True
        ]

        # Get matching ref2d with same quant and good group
        lquant = list(set([dd[kk]['quant'] for kk in lref1d]))
        dref2d = {
            k0: [
                kk for kk in _select(
                    dd=dd, quant=dd[k0]['quant'],
                    log='all', returnas=str,
                )
                if group2d in dd[kk]['group']
                and not isinstance(dd[kk]['data'], dict)
            ]
            for k0 in lref1d
        }
        dref2d = {k0: v0 for k0, v0 in dref2d.items() if len(v0) > 0}

        if len(dref2d) == 0:
            msg = (
                "No match for (ref1d, ref2d) for ddata['{}']".format(kq)
            )
            raise Exception(msg)

        # check ref1d
        if ref1d is None:
            if ref2d is not None:
                lk = [k0 for k0, v0 in dref2d.items() if ref2d in v0]
                if len(lk) == 0:
                    msg = (
                        "\nNon-valid interpolation intermediate\n"
                        + "\t- provided:\n"
                        + "\t\t- ref1d = {}, ref2d = {}\n".format(ref1d, ref2d)
                        + "\t- valid:\n{}".format(
                            '\n'.join([
                                '\t\t- ref1d = {}  =>  ref2d in {}'.format(
                                    k0, v0
                                )
                                for k0, v0 in dref2d.items()
                            ])
                        )
                    )
                    raise Exception(msg)
                if kq in lk:
                    ref1d = kq
                else:
                    ref1d = lk[0]
            else:
                if kq in dref2d.keys():
                    ref1d = kq
                else:
                    ref1d = list(dref2d.keys())[0]
        else:
            ref1d, msg = _get_keyingroup_ddata(
                dd=dd,
                key=ref1d, group=group1d,
                msgstr='ref1d', raise_=False,
            )
        if ref1d not in dref2d.keys():
            msg = (
                "\nNon-valid interpolation intermediate\n"
                + "\t- provided:\n"
                + "\t\t- ref1d = {}, ref2d = {}\n".format(ref1d, ref2d)
                + "\t- valid:\n{}".format(
                    '\n'.join([
                        '\t\t- ref1d = {}  =>  ref2d in {}'.format(
                            k0, v0
                        )
                        for k0, v0 in dref2d.items()
                    ])
                )
            )
            raise Exception(msg)

        # check ref2d
        if ref2d is None:
            ref2d = dref2d[ref1d][0]
        else:
            ref2d, msg = _get_keyingroup_ddata(
                dd=dd,
                key=ref2d, group=group2d,
                msgstr='ref2d', raise_=False,
            )
        if ref2d not in dref2d[ref1d]:
            msg = (
                "\nNon-valid interpolation intermediate\n"
                + "\t- provided:\n"
                + "\t\t- ref1d = {}, ref2d = {}\n".format(ref1d, ref2d)
                + "\t- valid:\n{}".format(
                    '\n'.join([
                        '\t\t- ref1d = {}  =>  ref2d in {}'.format(
                            k0, v0
                        )
                        for k0, v0 in dref2d.items()
                    ])
                )
            )
            raise Exception(msg)

    return kq, ref1d, ref2d
    coll=None,
    key=None,
    R=None,
    Z=None,
    indbs=None,
    indt=None,
    grid=None,
    details=None,
    reshape=None,
    res=None,
    crop=None,
    nan0=None,
    imshow=None,
    return_params=None,
):

    # ---------------
    # check inputs

    # TBD
    pass

    # ---------------
    # post-treatment

    if nan0 is True:
        val[val == 0] = np.nan

    # ------
    # return

    if return_params is True:
        return val, dparams
    else:
        return val


# #############################################################################
# #############################################################################
#                           Mesh2DRect - operators
# #############################################################################


def get_bsplines_operator(
    coll,
    key=None,
    operator=None,
    geometry=None,
    crop=None,
    store=None,
    returnas=None,
    # specific to deg = 0
    centered=None,
    # to return gradR, gradZ, for D1N2 deg 0, for tomotok
    returnas_element=None,
):

    # check inputs
    lk = list(coll.dobj.get('bsplines', {}).keys())
    key = ds._generic_check._check_var(
        key, 'key',
        types=str,
        allowed=lk,
    )

    store = ds._generic_check._check_var(
        store, 'store',
        default=True,
        types=bool,
    )

    returnas = ds._generic_check._check_var(
        returnas, 'returnas',
        default=store is False,
        types=bool,
    )

    crop = ds._generic_check._check_var(
        crop, 'crop',
        default=True,
        types=bool,
    )

    # cropbs
    cropbs = coll.dobj['bsplines'][key]['crop']
    keycropped = coll.dobj['bsplines'][key]['ref-bs'][0]
    if cropbs not in [None, False] and crop is True:
        cropbs_flat = coll.ddata[cropbs]['data'].ravel(order='F')
        if coll.dobj['bsplines'][key]['deg'] == 0:
            cropbs = coll.ddata[cropbs]['data']
        keycropped = f"{keycropped}-crop"
    else:
        cropbs = False
        cropbs_flat = False

    # compute and return
    (
        opmat, operator, geometry, dim,
    ) = coll.dobj['bsplines'][key]['class'].get_operator(
        operator=operator,
        geometry=geometry,
        cropbs_flat=cropbs_flat,
        # specific to deg=0
        cropbs=cropbs,
        centered=centered,
        # to return gradR, gradZ, for D1N2 deg 0, for tomotok
        returnas_element=returnas_element,
    )

    # cropping
    if operator == 'D1':
        ref = (keycropped, keycropped)
    elif operator == 'D0N1':
        ref = (keycropped,)
    elif 'N2' in operator:
        ref = (keycropped, keycropped)

    return opmat, operator, geometry, dim, ref, crop, store, returnas, key


# #################################################################
# #################################################################
#               Contour computation
# #################################################################


def _get_contours(
    RR=None,
    ZZ=None,
    val=None,
    levels=None,
    largest=None,
    uniform=None,
):
    """ Return R, Z coordinates of contours (time-dependent)

    For contourpy algorithm, the dimensions shoud be (ny, nx), from meshgrid

    RR = (nz, nr)
    ZZ = (nz, nr)
    val = (nt, nz, nr)
    levels = (nlevels,)

    cR = (nt, nlevels, nmax) array of R coordinates
    cZ = (nt, nlevels, nmax) array of Z coordinates

    The contour coordinates are uniformzied to always have the same nb of pts

    """

    # -------------
    # check inputs

    if largest is None:
        largest = False

    if uniform is None:
        uniform = True

    # val.shape = (nt, nR, nZ)
    lc = [
        val.shape == RR.shape,
        val.ndim == RR.ndim + 1 and val.shape[1:] == RR.shape,
    ]
    if lc[0]:
        val = val[None, ...]
    elif lc[1]:
        pass
    else:
        msg = "Incompatible val.shape!"
        raise Exception(msg)

    nt, nR, nZ = val.shape

    # ------------------------
    # Compute list of contours

    # compute contours at rknots
    # see https://github.com/matplotlib/matplotlib/blob/main/src/_contour.h

    contR = [[] for ii in range(nt)]
    contZ = [[] for ii in range(nt)]
    for ii in range(nt):

        # define map
        contgen = contour_generator(
            x=RR,
            y=ZZ,
            z=val[ii, ...],
            name='serial',
            corner_mask=None,
            line_type='Separate',
            fill_type=None,
            chunk_size=None,
            chunk_count=None,
            total_chunk_count=None,
            quad_as_tri=True,       # for sub-mesh precision
            # z_interp=<ZInterp.Linear: 1>,
            thread_count=0,
        )

        for jj in range(len(levels)):

            # compute concatenated contour
            no_cont = False
            cj = contgen.lines(levels[jj])

            c0 = (
                isinstance(cj, list)
                and all([
                    isinstance(cjj, np.ndarray)
                    and cjj.ndim == 2
                    and cjj.shape[1] == 2
                    for cjj in cj
                ])
            )
            if not c0:
                msg = f"Wrong output from contourpy!\n{cj}"
                raise Exception(msg)

            if len(cj) > 0:
                cj = [
                    cc[np.all(np.isfinite(cc), axis=1), :]
                    for cc in cj
                    if np.sum(np.all(np.isfinite(cc), axis=1)) >= 3
                ]

                if len(cj) == 0:
                    no_cont = True
                elif len(cj) == 1:
                    cj = cj[0]
                elif len(cj) > 1:
                    if largest:
                        nj = [
                            0.5*np.abs(np.sum(
                                (cc[1:, 0] + cc[:-1, 0])
                                *(cc[1:, 1] - cc[:-1, 1])
                            ))
                            for cc in cj
                        ]
                        cj = cj[np.argmax(nj)]
                    else:
                        ij = np.cumsum([cc.shape[0] for cc in cj])
                        cj = np.concatenate(cj, axis=0)
                        cj = np.insert(cj, ij, np.nan, axis=0)

                elif np.sum(np.all(~np.isnan(cc), axis=1)) < 3:
                    no_cont = True
            else:
                no_cont = True

            if no_cont is True:
                cj = np.full((3, 2), np.nan)

            contR[ii].append(cj[:, 0])
            contZ[ii].append(cj[:, 1])

    # ------------------------------------------------
    # Interpolate / concatenate to uniformize as array

    if uniform:
        ln = [[pp.size for pp in cc] for cc in contR]
        nmax = np.max(ln)
        cR = np.full((nt, len(levels), nmax), np.nan)
        cZ = np.full((nt, len(levels), nmax), np.nan)

        for ii in range(nt):
            for jj in range(len(levels)):
                cR[ii, jj, :] = np.interp(
                    np.linspace(0, ln[ii][jj], nmax),
                    np.arange(0, ln[ii][jj]),
                    contR[ii][jj],
                )
                cZ[ii, jj, :] = np.interp(
                    np.linspace(0, ln[ii][jj], nmax),
                    np.arange(0, ln[ii][jj]),
                    contZ[ii][jj],
                )

        return cR, cZ
    else:
        return contR, contZ


# #############################################################################
# #############################################################################
#                   Polygon simplification
# #############################################################################


def _simplify_polygon(pR=None, pZ=None, res=None):
    """ Use convex hull with a constraint on the maximum discrepancy """

    # ----------
    # preliminary 1: check there is non redundant point

    dp = np.sqrt((pR[1:] - pR[:-1])**2 + (pZ[1:] - pZ[:-1])**2)
    ind = (dp > 1.e-6).nonzero()[0]
    pR = pR[ind]
    pZ = pZ[ind]

    # check new poly is closed
    if (pR[0] != pR[-1]) or (pZ[0] != pZ[-1]):
        pR = np.append(pR, pR[0])
        pZ = np.append(pZ, pZ[0])

    # check it is counter-clockwise
    clock = np.nansum((pR[1:] - pR[:-1]) * (pZ[1:] + pZ[:-1]))
    if clock > 0:
        pR = pR[::-1]
        pZ = pZ[::-1]

    # threshold = diagonal of resolution + 10%
    thresh = res * np.sqrt(2) * 1.1

    # ----------
    # preliminary 2: get convex hull and copy

    poly = np.array([pR, pZ]).T
    iconv = ConvexHull(poly, incremental=False).vertices

    # close convex hull to iterate on edges
    pR_conv = np.append(pR[iconv], pR[iconv[0]])
    pZ_conv = np.append(pZ[iconv], pZ[iconv[0]])

    # copy to create new polygon that will serve as buffer
    pR_bis, pZ_bis = np.copy(pR), np.copy(pZ)

    # -------------------------
    # loop on convex hull edges

    for ii in range(pR_conv.size - 1):

        pR1, pR2 = pR_conv[ii], pR_conv[ii+1]
        pZ1, pZ2 = pZ_conv[ii], pZ_conv[ii+1]
        i0 = np.argmin(np.hypot(pR_bis - pR1, pZ_bis - pZ1))

        # make sure it starts from p1
        pR_bis = np.append(pR_bis[i0:], pR_bis[:i0])
        pZ_bis = np.append(pZ_bis[i0:], pZ_bis[:i0])

        # get indices of closest points to p1, p2
        i1 = np.argmin(np.hypot(pR_bis - pR1, pZ_bis - pZ1))
        i2 = np.argmin(np.hypot(pR_bis - pR2, pZ_bis - pZ2))

        # get corresponding indices of poly points to be included
        if i2 == i1 + 1:
            itemp = [i1, i2]

        else:
            # several points in-between
            # => check they are all within distance before exclusing them

            # get unit vector of segment
            norm12 = np.hypot(pR2 - pR1, pZ2 - pZ1)
            u12R = (pR2 - pR1) / norm12
            u12Z = (pZ2 - pZ1) / norm12

            # get points standing between p1 nd p2
            lpR = pR_bis[i1 + 1:i2]
            lpZ = pZ_bis[i1 + 1:i2]

            # indices of points standing too far from edge (use cross-product)
            iout = np.abs(u12R*(lpZ - pZ1) - u12Z*(lpR - pR1)) > thresh

            # if any pts too far => include all pts
            if np.any(iout):
                itemp = np.arange(i1, i2 + 1)
            else:
                itemp = [i1, i2]

        # build pts_in
        pR_in = pR_bis[itemp]
        pZ_in = pZ_bis[itemp]

        # concatenate to add to new polygon
        pR_bis = np.append(pR_in, pR_bis[i2 + 1:])
        pZ_bis = np.append(pZ_in, pZ_bis[i2 + 1:])

    # check new poly is closed
    if (pR_bis[0] != pR_bis[-1]) or (pZ_bis[0] != pZ_bis[-1]):
        pR_bis = np.append(pR_bis, pR_bis[0])
        pZ_bis = np.append(pZ_bis, pZ_bis[0])

    return pR_bis, pZ_bis


# #############################################################################
# #############################################################################
#                   radius2d special points handling
# #############################################################################


def radius2d_special_points(
    coll=None,
    key=None,
    keym0=None,
    res=None,
):

    keybs = coll.ddata[key]['bsplines']
    keym = coll.dobj['bsplines'][keybs]['mesh']
    mtype = coll.dobj[coll._which_mesh][keym]['type']
    assert mtype in ['rect', 'tri']

    # get map sampling
    RR, ZZ = coll.get_sample_mesh(
        key=keym,
        res=res,
        grid=True,
    )

    # get map
    val, t, _ = coll.interpolate_profile2d(
        key=key,
        R=RR,
        Z=ZZ,
        grid=False,
        imshow=True,        # for contour
    )

    # get min max values
    rmin = np.nanmin(val)
    rmax = np.nanmax(val)

    # get contour of 0
    cR, cZ = _get_contours(
        RR=RR,
        ZZ=ZZ,
        val=val,
        levels=[rmin + 0.05*(rmax-rmin)],
    )

    # dref
    ref_O = f'{keym0}-pts-O-n'
    dref = {
        ref_O: {'size': 1},
    }

    # get barycenter 
    if val.ndim == 3:
        assert cR.shape[1] == 1
        ax_R = np.nanmean(cR[:, 0, :], axis=-1)[:, None]
        ax_Z = np.nanmean(cZ[:, 0, :], axis=-1)[:, None]
        reft = coll.ddata[key]['ref'][0]
        ref = (reft, ref_O)
    else:
        ax_R = np.r_[np.nanmean(cR)]
        ax_Z = np.r_[np.nanmean(cZ)]
        ref = (ref_O,)

    kR = f'{keym0}-pts-O-R'
    kZ = f'{keym0}-pts-O-Z'
    ddata = {
        kR: {
            'ref': ref,
            'data': ax_R,
            'dim': 'distance',
            'quant': 'R',
            'name': 'O-points_R',
            'units': 'm',
        },
        kZ: {
            'ref': ref,
            'data': ax_Z,
            'dim': 'distance',
            'quant': 'Z',
            'name': 'O-points_Z',
            'units': 'm',
        },
    }

    return dref, ddata, kR, kZ


# #############################################################################
# #############################################################################
#                   angle2d discontinuity handling
# #############################################################################


def angle2d_zone(
    coll=None,
    key=None,
    keyrad2d=None,
    key_ptsO=None,
    res=None,
    keym0=None,
):

    keybs = coll.ddata[key]['bsplines']
    keym = coll.dobj['bsplines'][keybs]['mesh']
    mtype = coll.dobj[coll._which_mesh][keym]['type']
    assert mtype in ['rect', 'tri']

    # --------------
    # prepare

    hastime, hasvect, reft, keyt = coll.get_time(key=key)[:4]
    if hastime:
        nt = coll.dref[reft]['size']
    else:
        msg = (
            "Non time-dependent angle2d not implemented yet\n"
            "=> ping @Didou09 on Github to open an issue"
        )
        raise NotImplementedError(msg)

    if res is None:
        res = _get_sample_mesh_res(
            coll=coll,
            keym=keym,
            mtype=mtype,
        )

    # get map sampling
    RR, ZZ = coll.get_sample_mesh(
        key=keym,
        res=res/2.,
        grid=True,
        imshow=True,    # for contour
    )

    # get map
    val, t, _ = coll.interpolate_profile2d(
        key=key,
        R=RR,
        Z=ZZ,
        grid=False,
        azone=False,
    )
    val[np.isnan(val)] = 0.
    amin = np.nanmin(val)
    amax = np.nanmax(val)

    # get contours of absolute value
    cRmin, cZmin = _get_contours(
        RR=RR,
        ZZ=ZZ,
        val=val,
        levels=[amin + 0.10*(amax - amin)],
        largest=True,
        uniform=True,
    )
    cRmax, cZmax = _get_contours(
        RR=RR,
        ZZ=ZZ,
        val=val,
        levels=[amax - 0.10*(amax - amin)],
        largest=True,
        uniform=True,
    )

    cRmin, cZmin = cRmin[:, 0, :], cZmin[:, 0, :]
    cRmax, cZmax = cRmax[:, 0, :], cZmax[:, 0, :]

    rmin = np.full(cRmin.shape, np.nan)
    rmax = np.full(cRmax.shape, np.nan)

    # get points inside contour 
    for ii in range(nt):
        rmin[ii, :], _, _ = coll.interpolate_profile2d(
            key=keyrad2d,
            R=cRmin[ii, :],
            Z=cZmin[ii, :],
            grid=False,
            indt=ii,
        )
        rmax[ii, :], _, _ = coll.interpolate_profile2d(
            key=keyrad2d,
            R=cRmax[ii, :],
            Z=cZmax[ii, :],
            grid=False,
            indt=ii,
        )

    # get magnetic axis
    kR, kZ = key_ptsO
    axR = coll.ddata[kR]['data']
    axZ = coll.ddata[kZ]['data']
    assert coll.ddata[kR]['ref'][0] == coll.ddata[key]['ref'][0]

    start_min = np.nanargmin(rmin, axis=-1)
    start_max = np.nanargmin(rmax, axis=-1)

    # re-order from start_min, start_max
    lpR, lpZ = [], []
    for ii in range(rmin.shape[0]):
        imin = np.r_[
            np.arange(start_min[ii], rmin.shape[1]),
            np.arange(0, start_min[ii]),
        ]

        cRmin[ii] = cRmin[ii, imin]
        cZmin[ii] = cZmin[ii, imin]
        rmin[ii] = rmin[ii, imin]
        # check it is counter-clockwise
        clock = np.nansum(
            (cRmin[ii, 1:] - cRmin[ii, :-1])
            *(cZmin[ii, 1:] + cZmin[ii, :-1])
        )
        if clock > 0:
            cRmin[ii, :] = cRmin[ii, ::-1]
            cZmin[ii, :] = cZmin[ii, ::-1]
            rmin[ii, :] = rmin[ii, ::-1]

        imax = np.r_[
            np.arange(start_max[ii], rmax.shape[1]),
            np.arange(0, start_max[ii])
        ]
        cRmax[ii] = cRmax[ii, imax]
        cZmax[ii] = cZmax[ii, imax]
        rmax[ii] = rmax[ii, imax]
        # check it is clockwise
        clock = np.nansum(
            (cRmax[ii, 1:] - cRmax[ii, :-1])
            *(cZmax[ii, 1:] + cZmax[ii, :-1])
        )
        if clock < 0:
            cRmax[ii, :] = cRmax[ii, ::-1]
            cZmax[ii, :] = cZmax[ii, ::-1]
            rmax[ii, :] = rmax[ii, ::-1]

        # i0
        dr = np.diff(rmin[ii, :])
        i0 = (np.isnan(dr) | (dr < 0)).nonzero()[0][0]
        # rmin[ii, i0-1:] = np.nan
        dr = np.diff(rmax[ii, :])
        i1 = (np.isnan(dr) | (dr < 0)).nonzero()[0][0]
        # rmax[ii, i1-1:] = np.nan

        # polygon
        pR = np.r_[axR[ii], cRmin[ii, :i0-1], cRmax[ii, :i1-1][::-1]]
        pZ = np.r_[axZ[ii], cZmin[ii, :i0-1], cZmax[ii, :i1-1][::-1]]

        pR, pZ = _simplify_polygon(pR=pR, pZ=pZ, res=res)

        lpR.append(pR)
        lpZ.append(pZ)

    # Ajust sizes
    nb = np.array([pR.size for pR in lpR])

    # 
    nmax = np.max(nb)
    pR = np.full((nt, nmax), np.nan)
    pZ = np.full((nt, nmax), np.nan)

    for ii in range(nt):
        pR[ii, :] = np.interp(
            np.linspace(0, nb[ii], nmax),
            np.arange(0, nb[ii]),
            lpR[ii],
        )
        pZ[ii, :] = np.interp(
            np.linspace(0, nb[ii], nmax),
            np.arange(0, nb[ii]),
            lpZ[ii],
        )

    # ----------------
    # prepare output dict

    # ref
    kref = f'{keym0}-azone-npt'
    dref = {
        kref: {'size': nmax}
    }

    # data
    kR = f'{keym0}-azone-R'
    kZ = f'{keym0}-azone-Z'
    ddata = {
        kR: {
            'data': pR,
            'ref': (reft, kref),
            'units': 'm',
            'dim': 'distance',
            'quant': 'R',
            'name': None,
        },
        kZ: {
            'data': pZ,
            'ref': (reft, kref),
            'units': 'm',
            'dim': 'distance',
            'quant': 'R',
            'name': None,
        },
    }

    return dref, ddata, kR, kZ


def angle2d_inzone(
    coll=None,
    keym0=None,
    keya2d=None,
    R=None,
    Z=None,
    t=None,
    indt=None,
):


    # ------------
    # prepare points

    if R.ndim == 1:
        shape0 = None
        pts = np.array([R, Z]).T
    else:
        shape0 = R.shape
        pts = np.array([R.ravel(), Z.ravel()]).T

    # ------------
    # prepare path

    kazR, kazZ = coll.dobj[coll._which_mesh][keym0]['azone']
    pR = coll.ddata[kazR]['data']
    pZ = coll.ddata[kazZ]['data']

    hastime, hasvect, reft, keyt, tnew, dind = coll.get_time(
        key=kazR,
        t=t,
        indt=indt,
    )

    # ------------
    # test points

    if hastime:
        if dind is None:
            nt = coll.dref[reft]['size']
            ind = np.zeros((nt, R.size), dtype=bool)
            for ii in range(nt):
                path = Path(np.array([pR[ii, :], pZ[ii, :]]).T)
                ind[ii, :] = path.contains_points(pts)
        else:
            import pdb; pdb.set_trace()     # DB
            raise NotImplementedError()
            # TBC / TBF
            nt = None
            ind = np.zeros((nt, R.size), dtype=bool)
            for ii in range(nt):
                path = Path(np.array([pR[ii, :], pZ[ii, :]]).T)
                ind[ii, :] = path.contains_points(pts)

    else:
        path = Path(np.array([pR, pZ]).T)
        ind = path.contains_points(pts)

    # -------------------------
    # fromat output and return

    if shape0 is not None:
        if hastime:
            ind = ind.reshape(tuple(np.r_[nt, shape0]))
        else:
            ind = ind.reshape(shape0)

    return ind
