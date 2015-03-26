#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Main veritas interfaces.

:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2015
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
# The following two lines are the main reason why the code works in Python 2
# and Python 3.
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

from collections import namedtuple
import glob
import os

import numpy as np

from . import logger
from .io import read_specfem_stations_file


Channel = namedtuple("Channel", ["network", "station", "component",
                                 "latitude", "longitude", "filename_a",
                                 "filename_b"])


class WaveformDataSet(object):
    """
    Simple object helping to connect two waveform data sets and associated
    station meta information.

    By convention, dataset A is the one which is assumed to be more wrong
    than dataset B, e.g. the one calculated with the lower polynomial order,
    or no topography, ...
    """
    def __init__(self):
        self.dataset_a = {}
        self.dataset_b = {}
        self._stations = None

    def __iter__(self):
        self._iter_cur_index = 0
        self._iter_channels = self.all_channels
        return self

    def __next__(self):
        if self._iter_cur_index >= len(self._iter_channels):
            raise StopIteration
        cur_channel = self._iter_channels[self._iter_cur_index]
        self._iter_cur_index += 1
        return self.get(*cur_channel)

    def get_coordinates(self, network, station):
        row = self._stations[(self._stations["station"] == station) &
                             (self._stations["network"] == network)]
        return float(row.latitude), float(row.longitude)

    def get(self, network, station, component):
        ds_a = self.dataset_a[(network, station, component)]
        ds_b = self.dataset_b[(network, station, component)]
        latitude, longitude = self.get_coordinates(network, station)
        return Channel(network=network, station=station, component=component,
                       latitude=latitude,  longitude=longitude,
                       filename_a=ds_a, filename_b=ds_b)

    @property
    def all_channels(self):
        """
        Returns a sorted list of all channels that are part of both waveform
        datasets and have station information.
        """
        # These are the channels common in both waveform datasets.
        common_channels = self.common_channels

        all_channels = []
        # Loop over all stations that also have station information.
        for station in self.stations:
            channels = [_i for _i in common_channels if _i[:2] == station[:2]]
            all_channels.extend(channels)

        all_channels = sorted(all_channels)
        return all_channels

    @property
    def waveform_stations(self):
        """
        Get a set of all stations that are part of both datasets.
        """
        return set([(_i[0], _i[1]) for _i in self.common_channels])

    @property
    def station_stations(self):
        """
        Get a set of all stations that are part of the station information.
        """
        return set([(_i[1].network, _i[1].station) for _i in
                    self._stations.iterrows()])

    @property
    def stations(self):
        """
        Get a list of all stations that have station information, as well
        waveform data in both data sets.
        """
        return self.waveform_stations.intersection(self.station_stations)

    @property
    def common_channels(self):
        """
        Returns a set of common channel ids.
        """
        return self.channels_in_a.intersection(self.channels_in_b)

    @property
    def channels_in_a(self):
        """
        Returns a set of channels in dataset A.
        """
        return set(self.dataset_a.keys())

    @property
    def channels_in_b(self):
        """
        Returns a set of channels in dataset B.
        """
        return set(self.dataset_b.keys())

    @property
    def channels_only_in_set_a(self):
        """
        Returns a set of channels only in dataset A.
        """
        return self.channels_in_a.difference(self.common_channels)

    @property
    def channels_only_in_set_b(self):
        """
        Returns a set of channels only in dataset B.
        """
        return self.channels_in_b.difference(self.common_channels)

    def add_waveform_to_dataset_a(self, net_sta_comp, filename):
        """
        Add waveforms to dataset A.

        :param net_sta_comp: A tuple of network id, station id, and the
            component.
        :param filename: The filename of the waveform data.
        """
        # Convert to upper case to be a bit defensive.
        self.dataset_a[tuple([_i.upper() for _i in net_sta_comp])] = filename

    def add_waveform_to_dataset_b(self, net_sta_comp, filename):
        """
        Add waveforms to dataset B.

        :param net_sta_comp: A tuple of network id, station id, and the
            component.
        :param filename: The filename of the waveform data.
        """
        # Convert to upper case to be a bit defensive.
        self.dataset_b[tuple([_i.upper() for _i in net_sta_comp])] = filename

    def set_stations_dataframe(self, df):
        """
        Set the stations dataframe of the data set object.

        :param df: A pandas dataframe object with at least the following
            columns: ``"network"``, ``"station"``, ``"latitude"``,
            ``"longitude"``
        """
        self._stations = df


class Config(object):
    def __init__(self, low_res_seismos, high_res_seismos, station_info,
                 t_min, t_max, dt, get_net_sta_comp_fct=None):
        """

        :param low_res_seismos: A UNIX style wildcard pattern to find the
            low resolution seismograms.
        :type low_res_seismos: str
        :param high_res_seismos: A UNIX style wildcard pattern to find the
            high resolution seismograms.
        :type high_res_seismos: str
        :param station_info: A UNIX style wildcard pattern to find the
            station information.
        :type station_info: str
        :param t_min: The minimum period to test.
        :type t_min: float
        :param t_max: The maximum period to test.
        :type t_max: float
        :param dt: The delta for the tested frequencies.
        :type dt: float
        :param get_net_sta_comp_fct: A function that takes a seismogram
            filename and returns network, station, and component of that
            seismogram. Needed if the filename encodes that information.
            SPECFEM recently changed their naming scheme so no default is
            provided as it would be too dangerous to be wrong.
        :type get_net_sta_comp_fct: function
        """
        self.low_res_seismos = low_res_seismos
        self.high_res_seismos = high_res_seismos
        self.get_net_sta_comp_fct = get_net_sta_comp_fct
        self.wf_dataset = WaveformDataSet()

        # Get a list of frequencies to test. Make sure t_max is included.
        assert t_min < t_max
        self.frequencies = np.arange(t_min, t_max + dt * 0.1, dt)

        self._find_waveform_files()
        self.wf_dataset.set_stations_dataframe(read_specfem_stations_file(
            station_info))

        avail_stations = self.wf_dataset.stations
        wf_s = self.wf_dataset.waveform_stations
        logger.info("%i stations are part of both datasets and also have "
                    "available station information." %
                    len(avail_stations))
        if len(wf_s) > len(avail_stations):
            logger.info("%i stations are part of both datasets but no "
                        "station information exists for them." % (
                len(wf_s) - len(avail_stations)))

        for _i in self.wf_dataset:
            print(_i)

    def _find_waveform_files(self):
        """
        Finds the waveform files for the low and high resolution seismograms.
        """
        low_res = glob.glob(self.low_res_seismos)
        high_res = glob.glob(self.high_res_seismos)

        for filename in low_res:
            self.wf_dataset.add_waveform_to_dataset_a(
                self.get_net_sta_comp_fct(os.path.basename(filename)),
                filename)

        for filename in high_res:
            self.wf_dataset.add_waveform_to_dataset_b(
                self.get_net_sta_comp_fct(os.path.basename(filename)),
                filename)


        c_chan = self.wf_dataset.common_channels
        a_chan = self.wf_dataset.channels_only_in_set_a
        b_chan = self.wf_dataset.channels_only_in_set_b

        logger.info("Found %i waveforms files that are part of both "
                    "data sets." % len(c_chan))

        # Give some additional information in case something is potentially
        # wrong.
        if a_chan:
            logger.warning("%i waveforms of the low resolution data set will "
                           "not be used as the high resolution data set does "
                           "not contain them." % (len(a_chan)))
        if b_chan:
            logger.warning("%i waveforms of the high resolution data set will "
                           "not be used as the low resolution data set does "
                           "not contain them." % (len(b_chan)))
