#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Interface to the C SATSI library.

:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2014
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from future.utils import native_str

import ctypes as C
import glob
import inspect
import itertools
import numpy as np
import multiprocessing
import os


# Don't use more than 8 cores to avoid funny situations on clusters.
MAX_NUMBER_OF_CORES = min(multiprocessing.cpu_count(), 8)


LIB_DIR = os.path.join(os.path.dirname(os.path.abspath(
    inspect.getfile(inspect.currentframe()))), "lib")

cache = []


def load_lib():
    if cache:
        return cache[0]
    else:
        # Enable a couple of different library naming schemes.
        possible_files = glob.glob(os.path.join(LIB_DIR, "satsi*.so"))
        if not possible_files:
            raise ValueError("Could not find suitable satsi shared "
                             "library.")
        filename = possible_files[0]
        lib = C.CDLL(filename)
        cache.append(lib)
        return lib


clibsatsi = load_lib()


# Results from the 2D satsi tradeoff calculation.
class Satsi2DTradeoffResult(C.Structure):
    _fields_ = [('mech_misfit', C.c_double),
                ('mvar', C.c_double)]


clibsatsi.satsi_2D_tradeoff.argtypes = [
    # x_in
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1,
                           flags=native_str('C_CONTIGUOUS')),
    # y_in
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1,
                           flags=native_str('C_CONTIGUOUS')),
    # ddir_in
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,
                           flags=native_str('C_CONTIGUOUS')),
    # dip_in
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,
                           flags=native_str('C_CONTIGUOUS')),
    # rake_in
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,
                           flags=native_str('C_CONTIGUOUS')),
    # input_length
    C.c_int,
    # cwt (damping_parameter)
    C.c_double]


clibsatsi.satsi_2D_tradeoff.restype = Satsi2DTradeoffResult


def _apply_calculate_2D_tradeoff(args):
    fault_plane_solutions, damping = args
    result = clibsatsi.satsi_2D_tradeoff(
        fault_plane_solutions["x"].values,
        fault_plane_solutions["y"].values,
        fault_plane_solutions["dip"].values,
        fault_plane_solutions["dip_angle"].values,
        fault_plane_solutions["rake"].values,
        len(fault_plane_solutions),
        float(damping))
    return (result.mech_misfit, result.mvar)


def calculate_2D_tradeoff(fault_plane_solutions):
    """
    Calculates the tradeoff between misfit and model length in parallel.

    It scales almost linear with the amount of processors, the cost of
    forking is definitely worth it. Using hyperthreading is also worth it.
    It thus always uses the amount of physical cores on a machine up to a
    maximum of 8 to avoid funny situations on clusters.

    :param fault_plane_solutions:
    """
    damp_parameters = [0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8,
                       2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 5.0, 6.0]
    tradeoffs = []
    model_variances = []

    # Apply in parallel as potentially pretty slow.
    with multiprocessing.Pool(
            processes=min(len(damp_parameters), MAX_NUMBER_OF_CORES)) as pool:
        result = pool.map(_apply_calculate_2D_tradeoff,
                          list(zip(itertools.repeat(fault_plane_solutions),
                               damp_parameters)))

    for misfit, variance in result:
        tradeoffs.append(misfit)
        model_variances.append(variance)

    return {
        "damping_value": np.float64(damp_parameters),
        "data_misfit": np.float64(tradeoffs),
        "model_length": np.float64(model_variances),
    }
