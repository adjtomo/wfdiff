#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Various misfit measurements.

All of them take two traces as input files.

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


def l2_norm(tr1, tr2):
    return np.sum((np.sum(tr1.data ** 2) - np.sum(tr2.data ** 2)) ** 2)


def l1_norm(tr1, tr2):
    return np.sum(np.abs(np.sum(np.abs(tr1.data)) - np.sum(np.abs(tr2.data))))


def x_corr(tr1, tr2):
    d = tr1.data
    s = tr2.data
    cc = np.correlate(d, s, mode="full")
    time_shift = cc.argmax() - len(d) + 1
    # Normalized cross correlation.
    max_cc_value = cc.max() / np.sqrt((s ** 2).sum() * (d ** 2).sum())
    return max_cc_value, time_shift
