#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
I/O Utilities for wfdiff.

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
import obspy
import pandas


def read_specfem_stations_file(filename):
    """
    Reads a SPECFEM stations file to a pandas data frame.

    :param filename: The file to read.
    """
    data = pandas.io.parsers.read_table(
        filename, sep=r"\s+", header=None,
        names=["station", "network", "latitude", "longitude", "elevation",
               "depth"])
    return data


def read_specfem_ascii_waveform_file(filename, network, station, channel):
    """
    Reads SPECFEM ASCII files to a :class:`~obspy.core.stream.Stream` object.

    :param filename: The filename.
    :type filename: str
    :param network: The network id of the data.
    :type network: str
    :param station: The station id of the data.
    :type station: str
    :param channel: The channel id of the data.
    :type channel: str
    """
    time_array, data = np.loadtxt(filename).T
    # Try to get a reasonably accurate sample spacing.
    dt = np.diff(time_array).mean()

    tr = obspy.Trace(data=data)
    tr.stats.network = network
    tr.stats.station = station
    tr.stats.channel = channel
    tr.stats.delta = dt
    tr.stats.starttime += time_array[0]

    return obspy.Stream(traces=[tr])
