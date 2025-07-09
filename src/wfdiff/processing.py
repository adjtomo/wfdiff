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

import obspy
import numpy as np


UNITS_MAP = {
    "displacement": 0,
    "velocity": 1,
    "acceleration": 2
}


def preprocess_traces(tr_a, tr_b, data_units, desired_units, starttime=None,
                      endtime=None):
    """
    Makes sure both traces are sampled at the same points in time.

    Right now it is fairly simple and just interpolates the lower sampled
    trace to the higher sampled one. This is simple but should be fairly
    stable.

    It furthermore cuts the data to the desired start and end times in
    seconds, converts to the desired units and detrends, demeans, and tapers
    to stabilize the filters afterwards.

    The traces will be changed in place.

    :param tr_a: Trace A.
    :type tr_a: :class:`~obspy.core.trace.Trace`
    :param tr_b: Trace B.
    :type tr_b: :class:`~obspy.core.trace.Trace`
    :param data_units: The units of the data. One of ``"displacement"``,
        ``"velocity"``, or ``"acceleration"``.
    :type data_units: str
    :param desired_units: The desired output units of this function. One of
        ``"displacement"``, ``"velocity"``, or ``"acceleration"``.
    :type desired_units: str
    :param starttime: Time in seconds since the first common sample of both
        traces after the interpolation. If given, data will be cut to that
        time.
    :type starttime: float
    :param endtime: Time in seconds since the first common sample of both
        traces after the interpolation. If given, data will be cut to that
        time.
    :type endtime: float
    """
    # Make life easy and always interpolate the lower resolution
    # one to the higher resolution one.
    low_res, high_res = sorted([tr_a, tr_b], key=lambda x: x.stats.npts)

    # Assert the high res one is completely contained in the low
    # res one.
    high_res.trim(low_res.stats.starttime, low_res.stats.endtime,
                  nearest_sample=False)

    # Cubic spline interpolation. After this points, both traces will be
    # sampled at exactly the same points in time.
    # XXX: Maybe use Lanczos resampling?
    low_res.interpolate(sampling_rate=high_res.stats.sampling_rate,
                        method="cubic", starttime=high_res.stats.starttime,
                        npts=high_res.stats.npts)

    # Assemble to dummy stream to facilitate the further processing.
    st = obspy.Stream(traces=[low_res, high_res])

    # Transform to the desired units. Differentiation is a second order
    # finite difference scheme and integration is performed by the
    # trapezoidal rule.
    # This is done before the trimming as it potentially reduces boundary
    # effects and the operations are very cheap in any case.
    # XXX: think about higher order / spectral domain operators for more
    # accuracy but this likely has little influence on what this package is
    # supposed to calculate.

    # Positive means a derivative, negative an anti-derivative. Perform as
    # many times as needed.
    derivatives = UNITS_MAP[desired_units] - UNITS_MAP[data_units]
    for _ in range(derivatives):
        if np.sign(derivatives) == 1:
            st.differentiate()
        else:
            st.integrate()

    # Detrend, demean, taper to stabilize the filters that will follow.
    st.detrend("demean")
    st.detrend("linear")
    st.taper(type="hann", max_percentage=0.03)

    if starttime is not None:
        starttime = st[0].stats.starttime + starttime
    if endtime is not None:
        endtime = st[0].stats.starttime + endtime
    st.trim(starttime=starttime, endtime=endtime, pad=True, fill_value=0)
