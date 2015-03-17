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
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

import glob
import os

import numpy as np

from . import logger


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
        self.dataset_a[tuple(net_sta_comp)] = filename

    def add_waveform_to_dataset_b(self, net_sta_comp, filename):
        self.dataset_b[tuple(net_sta_comp)] = filename


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
