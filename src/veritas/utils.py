#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utility functionality.

:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2015
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA


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

    # Cubic spline interpolation.
    # XXX: Maybe use Lanczos resampling?
    low_res.interpolate(sampling_rate=high_res.stats.sampling_rate,
                        method="cubic", starttime=high_res.stats.starttime,
                        npts=high_res.stats.npts)
