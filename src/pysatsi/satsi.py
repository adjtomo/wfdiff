#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2014
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

import numpy as np
import pandas.io.parsers

from . import PysatsiError


def read_fault_plane_solutions(filename):
    # Read first line.
    with open(filename, "rt") as fh:
        first_line = fh.readline().strip().split()
    # The fault plane solutions file is either 5 or 7 columns long. 5
    # columns results in them being interpreted as a 2D list of time
    # independent values. 7 columns are time dependent 3D values.
    if len(first_line) == 5:
        headers = ["x", "y", "dip", "dip_angle", "rake"]
    elif len(first_line) == 7:
        headers = ["x", "y", "z", "t", "dip", "dip_angle", "rake"]
    else:
        raise PysatsiError("fault plane solutions file has neither 5 nor 7 "
                           "columns.")
    # All the C routines use double precision so values will be double
    # precision from the beginning.
    df = pandas.io.parsers.read_csv(filename, sep=" ", skipinitialspace=True,
                                    names=headers, dtype=np.float64)
    return df


def invert_stress(fault_plane_solutions, damping=True,
                  confidence_level=0.95, fraction_valid_fault_planes=0.5,
                  min_event_per_node=20, boostrap_resamplings=2000,
                  time_space_damping_ratio=1.0):
    """
    Invert for the stress tensor.

    :param fault_plane_solutions: The fault plane solutions.
    :param damping: Determines if damping is used in the LSQR. If ``True``
        the damping parameter will be estimated based on
    :param confidence_level:
    :param fraction_valid_fault_planes:
    :param min_event_per_node:
    :param boostrap_resamplings:
    :param time_space_damping_ratio:
    """
    pass