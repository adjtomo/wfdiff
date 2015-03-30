#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Various misfit measurements.

All of them take two traces as input files. For the purpose of
normalization, the first trace is intended to be the low resolution one,
and the second one the high resolution one.

All of them have to either return a dictionary or a list of dictionaries.
Each dictionary is a separate misfit measurement and has to have the
following five keys:

* ``"name"``: Internally used name. Usually the function name. Snake case
    please.
* ``"pretty_name"``: Verbose pretty name for the misfit measurement. Used
    for plotting and similar purposes.
* ``"value"``: Single float denoting the value of the misfit measurement.
* ``"logarithmic_plot"``: Boolean flag determining if the values should be
    plotted with a logarithmic scale or not.
* ``"minimizing_misfit"``: Boolean flag determining if a misfit is a
    classical minimizing misfit, e.g. smaller values are better or not.
    Cross correlations for examples are not, most others are.

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


# Add all misfits function here! Otherwise they will not be discovered.
__all__ = ["l2_norm", "l1_norm", "cross_correlation"]


def l2_norm(tr1, tr2):
    """
    L2 norm normalized by the integrated energy of the second trace.
    """
    return {
        "name": "l2_norm",
        "pretty_name": "Normalized L2 Norm",
        "value": np.log10(np.sum((tr1.data - tr2.data) ** 2) / np.sum(
            tr2.data ** 2)),
        "logarithmic_plot": False,
        "minimizing_misfit": True
    }


def l1_norm(tr1, tr2):
    """
    Normalized L1 norm between two traces. Reference is the L1 norm of the
    second trace.
    """
    return {
        "name": "l1_norm",
        "pretty_name": "Normalized L1 Norm",
        "logarithmic_plot": False,
        "value": np.abs(tr1.data - tr2.data).sum() / np.sum(np.abs(tr2.data)),
        "minimizing_misfit": True
    }


def cross_correlation(tr1, tr2):
    """
    Normalized cross correlation between two traces.
    """
    d = tr1.data
    s = tr2.data
    cc = np.correlate(d, s, mode="full")
    # Normalized cross correlation.
    max_cc_value = cc.max() / np.sqrt((s ** 2).sum() * (d ** 2).sum())
    return {
            "name": "cross_correlation",
            "pretty_name": "Cross Correlation Coefficient",
            "logarithmic_plot": False,
            "value": max_cc_value,
            # The larger the correlation, the better.
            "minimizing_misfit": False
    }
