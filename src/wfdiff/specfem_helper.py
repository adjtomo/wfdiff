#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Helper utilites for SPECFEM
1. Most of these are also used by wfdiff

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
import glob
import os

def read_specfem_cmtsolution_file(filename):
    """
    Read a SPECFEM CMTSOLUTION file
    """
    with open(filename) as f:
        content = f.readlines()
    c = [x.strip() for x in content]

    # create event object
    t = Tensor()
    t.m_rr, t.m_tt, t.m_pp, t.m_rt, t.m_rp, t.m_tp = float(c[7].split()[1]), float(c[8].split()[1]), float(c[9].split()[1]),float(c[10].split()[1]), float(c[11].split()[1]),float(c[12].split()[1])    # Moment tensor
    mt = MomentTensor()
    mt.tensor = t
    org = Origin()
    org.time = UTCDateTime(int(c[0].split()[1]), int(c[0].split()[2]), 
                           int(c[0].split()[3]), int(c[0].split()[4]), 
                           int(c[0].split()[5]), float(c[0].split()[6]))    # origin time
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
    evname = c[1].split()[2]   # event name used in the simulation
    tshift = c[2].split()[2]   # time-shift
    hdur = c[3].split()[2]     # half-duration

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


# SPECFEM encoded information in the filenames and the way this
# is done differs depending on the SPECFEM version. This function
# will be used to extract network id, station id, and the component
# from any waveform file encountered in case SPECFEM ASCII output
# is used.
def get_net_sta_comp(filename):
    net, sta, chan, _ = filename.split(".")
    return net, sta, chan[-1]


#  rotate to RTZ
def add_event_info(ev, stn, waveforms):
    '''
    Rotate SPECFEM files to RTZ
    '''
    
    ddir = os.path.dirname(waveforms)

    st = obspy.Stream()
    # Read waveforms
    for indx, row in stn.iterrows():        
        stnm = row['station']
        net = row['network'] 
        # Add backazimuth and inclination information
        res = obspy.geodetics.base.gps2dist_azimuth(
            ev.origins[0].latitude, ev.origins[0].longitude,
            row['latitude'], row['longitude'])

        baz = res[2] # backazimuth
        inc = 0      # inclination at all synthetic sites

        # Read specfem files for which station info is present
        fname = ddir + '/' + net + '.' + stnm + '.' + '*'
        print(fname)
        for f in glob.glob(fname):
            filename = os.path.basename(f)
            _, _, cha, _ = filename.split('.')
            cha = cha[0:2]
            net, sta, comp = get_net_sta_comp(filename)

            chan =  cha + comp
            tr_E = read_specfem_ascii_waveform_file(ddir + '/' + filename, net, stnm, chan)[0]
            tr_E.stats.starttime = ev.origins[0].time
            tr_E.stats.back_azimuth = baz
            tr_E.stats.inclination = inc
            st.append(tr_E)

    eid = otime2eid(ev.origins[0].time)
    return(st)


def save_as_sac(st, dir_name):
    '''
    Save SPECFEM files as SAC files
    '''
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    
    # write sac waveforms
    print("Saving waveforms in %s" % dir_name)

    for tr in st.traces:
        filename = dir_name + '/' \
                + tr.stats.network + '.' + tr.stats.station + '.' \
                + tr.stats.channel + '.sac'
        tr.write(filename, format='SAC')
    

def otime2eid(otime):
    """
    Convert origin time to origin ID. The origin ID has format: YYYYMMDDHHMMSSsss
    See also eid2otime.

    Example
        otime = "2009-04-07T20:12:55"
        eid = otime2eid(otime)
        print(otime, eid)
    """

    yy = UTCDateTime(otime).year
    mo = UTCDateTime(otime).month
    dd = UTCDateTime(otime).day
    hh = UTCDateTime(otime).hour
    mm = UTCDateTime(otime).minute
    ss = UTCDateTime(otime).second
    ms = int(UTCDateTime(otime).microsecond / 1000.0) # mili?
    eid = '%04d%02d%02d%02d%02d%02d%03d' % (yy, mo, dd, hh, mm, ss, ms)
    return eid
