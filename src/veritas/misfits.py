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


def _x_corr(tr1, tr2):
    d = tr1.data
    s = tr2.data
    cc = np.correlate(d, s, mode="full")
    time_shift = cc.argmax() - len(d) + 1
    # Normalized cross correlation.
    max_cc_value = cc.max() / np.sqrt((s ** 2).sum() * (d ** 2).sum())
    return max_cc_value, time_shift


def x_corr_value(tr1, tr2):
    return _x_corr(tr1, tr2)[0]


def x_corr_time_shift(tr1, tr2):
    return _x_corr(tr1, tr2)[1]


def preprocess_traces(tr_a, tr_b):
    """
    Makes sure both traces are sampled at the same points in time.

    Right now it is fairly simple and just interpolates the lower sampled
    trace to the higher sampled one. This is simple but should be fairly
    stable.

    The traces will be changed in place.

    :param tr_a: Trace A.
    :type tr_a: :class:`~obspy.core.trace.Trace`
    :param tr_b: Trace B.
    :type tr_b: :class:`~obspy.core.trace.Trace`
    """
    # Make life easy and always interpolate the lower resolution
    # one to the higher resolution one.
    low_res, high_res = sorted([tr_a, tr_b], key=lambda x: x.stats.npts)

    # Assert the high res one is completely contained in the low
    # res one.
    high_res.trim(low_res.stats.starttime, low_res.stats.endtime,
                  nearest_sample=False)
    low_res.interpolate(sampling_rate=high_res.stats.sampling_rate,
                        method="cubic", starttime=high_res.stats.starttime,
                        npts=high_res.stats.npts)
