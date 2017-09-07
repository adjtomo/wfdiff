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
from obspy.core import UTCDateTime
from obspy.core.event import Event, Origin, Magnitude, FocalMechanism, MomentTensor, Tensor
import pandas

def read_specfem_cmtsolution_file(filename):
    """
    Read a SPECFEM CMTSOLUTION file
    """
    with open(filename) as f:
        content = f.readlines()
    c = [x.strip() for x in content]
    print(c)

    # create event object
    t = Tensor()
    t.m_rr, t.m_tt, t.m_pp, t.m_rt, t.m_rp, t.m_tp = float(c[7].split()[1]), float(c[8].split()[1]), float(c[9].split()[1]),float(c[10].split()[1]), float(c[11].split()[1]),float(c[12].split()[1])
    mt = MomentTensor()
    mt.tensor = t
    org = Origin()
    org.time = UTCDateTime(int(c[0].split()[1]), int(c[0].split()[2]), 
                           int(c[0].split()[3]), int(c[0].split()[4]), 
                           int(c[0].split()[5]), float(c[0].split()[6]))   # Moment tensor
    org.latitude = c[4].split()[1]     # event latitude
    org.longitude = c[5].split()[1]    # event longitude
    org.depth = c[6].split()[1]        # event depth 
    mag = Magnitude()
    mag.mag = c[0].split()[10]         # magnitude
    mag.magnitude_type = "Mw"
    ev = Event()
    ev.origins.append(org)
    ev.magnitudes.append(mag)
    ev.focal_mechanisms.append(mt)
    
    # Not attributes of the event object
    hdur = c[3].split()[2]     # half-duration
    tshift = c[2].split()[2]   # time-shift
    evname = c[1].split()[2]   # event name used in the simulation

    return ev

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
