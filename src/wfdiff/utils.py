#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utility functionality for wfdiff.

:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2015
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

import numpy as np

def rightmost_threshold_crossing(
        periods, misfit_values,
        threshold, threshold_is_upper_limit):
    """
    Function finding the rightmost threshold crossing in a fairly stable
    manner.

    Returns the period and misfit value of the crossing.

    :param periods: The periods where the misfits have been measured.
    :param misfit_values: The associated misfit values.
    :param threshold: The threshold.
    :type threshold: float
    :param threshold_is_upper_limit: Flag to determine if the threshold is
        an upper limit or not. Should be ``True`` when a small misfit means a
        better fit to the data. ``False`` otherwise (e.g. for cross
        correlations).
    :type threshold_is_upper_limit: bool
    """
    periods = np.array(periods)
    misfit_values = np.array(misfit_values)

    # Mirror on threshold. Then we only have to deal
    # with one case in the following.
    if not threshold_is_upper_limit:
        misfit_values -= threshold
        misfit_values *= -1.0
        misfit_values += threshold

    # Step 1: Check rightmost value.
    if misfit_values[-1] >= threshold:
        return periods[-1], misfit_values[-1]

    # Step 2: Check all other values.
    if np.all(misfit_values <= threshold):
        return periods[0], misfit_values[0]

    # Step 3: Find the value that first crosses the threshold.
    c_idx = np.where(misfit_values >= threshold)[0][-1]

    # Step 4: Interpolate the two adjacent values to get a fairly
    # granular threshold intersection measurement.
    ip = np.linspace(periods[c_idx], periods[c_idx + 1], 100)
    im = np.interp(ip,
                   periods[c_idx: c_idx + 2],
                   misfit_values[c_idx: c_idx + 2])
    idx = np.argmin(np.abs(im - threshold))

    return ip[idx], im[idx]
