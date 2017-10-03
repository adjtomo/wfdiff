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
#from tqdm import tqdm
try:
    import pyasdf     # For reading asdf waveform files
except ImportError:
    pass


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

def read_specfem_files(specfem_files, new_format=True):
    '''
    Read specfem ascii files and return a stream
    '''
    ddir = os.path.dirname(specfem_files)
    st = obspy.Stream()
    
    #print('Reading specfem ASCII files')
    for f in glob.glob(specfem_files):
        filename = os.path.basename(f)
        net, sta, chan = get_net_sta_comp(filename, new_format)
        tr = read_specfem_ascii_waveform_file(ddir + '/' + filename, net, sta, chan)[0]
        st.append(tr)

    return st

# SPECFEM encoded information in the filenames and the way this
# is done differs depending on the SPECFEM version. This function
# will be used to extract network id, station id, and the component
# from any waveform file encountered in case SPECFEM ASCII output
# is used.
def get_net_sta_comp(filename, new_format=True):
    '''
    :param filename: specfem filename
    :param new_format: `True` for filename in NET.STA.CHA format
            `False` for filename in STA.NET.CHA format
    '''
    if new_format:
        net, sta, chan, _ = filename.split(".")
    else:
        sta, net, chan, _ = filename.split(".")
    
    return net, sta, chan


def add_event_station_info(st, event, stations):
    '''
    Add backaimuth and inclination to traces. Needed for rotation
    Will add info to traces which are also in the staion file

    :param event: event object
    :param stations: specfem station file read as pandas
    :param st: obspy stream
    '''
    
    # Read waveforms
    for indxs, row in stations.iterrows():        
        sta = row['station']
        net = row['network'] 
        # Add backazimuth and inclination information
        _, _, baz = obspy.geodetics.base.gps2dist_azimuth(
            event.origins[0].latitude, event.origins[0].longitude,
            row['latitude'], row['longitude'])

        # Select stream for this station
        try:
            st_net_sta = st.select(network=row['network'], 
                                   station=row['station'])
        except:
            print('No traces found for', row['network'], row['station'])

        for tr in st_net_sta:
            tr.stats.back_azimuth = baz
            tr.stats.starttime = event.origins[0].time

    return(st)


def save_as_sac(st, dir_name):
    '''
    Save SPECFEM files as SAC files

    :param st: obspy stream
    :param dir_name: output sac directory name
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


def specfem_to_asdf(asdf_filename, folder, stations_file,
                    event_file, wf_tag, new_format = True):
    '''
    convert specfem files into asdf

    :param asdf_filename: name for asdf dataset
    :param folder: folder containing semd files
    :param stations_file: specfem stations file in pandas dataframe format
    :wf_tag: tag for waveforms in asdf dataset
    :new_format: `True` (default) for new specfem name 
            format NET.STA.CHA; `False` for old name STA.NET.CHA
    '''

    def populate_channels(net, specfem_files):
        # Add channel to the station xml
        for filename in list(specfem_files):
            _, _, cha, _ = os.path.basename(filename).split(".")
            if cha[-1] is 'E':
                azimuth=90.0
                dip=0.0
            elif cha[-1] is 'N':
                azimuth=0.0
                dip=0.0
            elif cha[-1] is 'Z':
                azimuth=0.0
                dip=-90.0

            chan = obspy.core.inventory.Channel(
                code=cha,
                location_code="",
                latitude=s["latitude"],
                longitude=s["longitude"],
                elevation=s["elevation"],
                depth=s["depth"],
                azimuth=azimuth,
                dip=dip)
            net.stations[0].channels.append(chan)
        return(net)

    
    files = glob.glob(os.path.join(folder, "*.semd"))
    assert files

    cat = obspy.read_events(event_file)
    assert len(cat) == 1
    event = cat[0]

    with pyasdf.ASDFDataSet(asdf_filename) as ds:
        stations = read_specfem_stations_file(stations_file)

        # Add event info
        try:
            ds.add_quakeml(event)
        except:
            print('Event already present, skipping')

        # Add station info
        for s in stations.iterrows():
            s = s[1]
            net = obspy.core.inventory.Network(code=s["network"], stations=[
                obspy.core.inventory.Station(
                    creation_date=obspy.UTCDateTime(),
                    code=s["station"],
                    latitude=s["latitude"],
                    longitude=s["longitude"],
                    elevation=s["elevation"],
                    site=obspy.core.inventory.Site(name=""))])

            # Add channel info
            try:
                if new_format is True:
                    specfem_files = glob.glob(os.path.join(folder, s["network"] +
                                                          '.' + s["station"] + '.' + '*.semd'))
                else:
                    specfem_files = glob.glob(os.path.join(folder + "/" + s["station"] + 
                                                          '.' + s["network"] + '.' + "*.semd"))
                populate_channels(net, specfem_files)
            except:
                print(s["network"] +'.' + s["station"] + ': waveform file absent')
            
            ds.add_stationxml(obspy.core.inventory.Inventory(
                    networks=[net], source=""))

        # Add waveforms 
        st = read_specfem_files(os.path.join(folder, "*.semd"), 
                                new_format=new_format)
        st = add_event_station_info(st, event, stations)
        ds.add_waveforms(st, tag=wf_tag, event_id=event)


def get_stream_from_asdf(ds, wf_tag):
    '''
    Read asdf file and convert into obspy stream

    :param ds: asdf dataset 
    :param wf_tag: tag for the waveforms, ds.waveforms.wf_tag 
            (example: synthetic, data, gll5, gll7, etc)
    '''

    st = obspy.Stream()
    
    for sta in ds.waveforms:
        tag = 'sta.' + wf_tag
        st = st + eval(tag)

    return(st)


def get_stations_from_asdf(ds):
    '''
    Get staions as pandas dataframe from asdf file

    :param ds: asdf dataset
    '''

    net, sta, lat, lon, ele, dep = [], [], [], [], [], []

    for i, tag in enumerate(ds.waveforms.list()):
        n, s = tag.split('.')
        net.append(n)
        sta.append(s)

        ds_tag = n + '_' + s
        lat.append(ds.waveforms[ds_tag].StationXML.networks[0].stations[0].latitude)
        lon.append(ds.waveforms[ds_tag].StationXML.networks[0].stations[0].longitude)
        ele.append(ds.waveforms[ds_tag].StationXML.networks[0].stations[0].elevation)
        dep.append(ds.waveforms[ds_tag].StationXML.networks[0].stations[0].channels[0].depth)
            
    # Create pandas dataframe
    df = pandas.DataFrame(
        {"network": net,
         "station": sta,
         "latitude": lat,
         "longitude": lon,
         "elevation": ele,
         "depth": dep,
         })
        
    return df


def get_station_info(input_stations):
    if type(input_stations) is pyasdf.asdf_data_set.ASDFDataSet:
        df = get_stations_from_asdf(input_stations)
    elif type(input_stations) is str:
        df = read_specfem_stations_file(input_stations)

    return df
