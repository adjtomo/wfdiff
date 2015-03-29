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
import json
import os

import obspy
from mpi4py import MPI
import numpy as np

from . import logger
from . import misfits
from .io import read_specfem_stations_file, read_specfem_ascii_waveform_file

COMM = MPI.COMM_WORLD

Channel = namedtuple("Channel", ["network", "station", "component",
                                 "latitude", "longitude", "filename_high",
                                 "filename_low"])

class NumPyJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return json.dumps(list(map(float, obj)))
        return json.JSONEncoder.default(self, obj)


class Results(object):
    def __init__(self):
        self.__results = {}

    def add_result(self, result):
        self.__results[(result["network"], result["station"],
                        result["component"])] = result

    def write(self, filename):
        _results = {}
        for _i in self.__results.values():
            _results["{network}.{station}.{component}".format(**_i)] = _i

        with open(filename, "w") as fh:
            json.dump(_results, fh, sort_keys=True, indent=4,
                      separators=(",", ": "), cls=NumPyJSONEncoder)


MISFIT_MAP = {
    "l2_norm": misfits.l2_norm,
    "l1_norm": misfits.l1_norm,
    "x_corr_time_shift": misfits.x_corr_time_shift,
    "x_corr_value": misfits.x_corr_value
}


class WaveformDataSet(object):
    """
    Simple object helping to connect two waveform data sets and associated
    station meta information.

    `dataset_high` is the one with the "truer" value, be it higher
    polynomial degrees during the simulation or more accuracy due to other
    reasons.
    """
    def __init__(self):
        self.dataset_high = {}
        self.dataset_low = {}
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
        ds_high = self.dataset_high[(network, station, component)]
        ds_low = self.dataset_low[(network, station, component)]
        latitude, longitude = self.get_coordinates(network, station)
        return Channel(network=network, station=station, component=component,
                       latitude=latitude,  longitude=longitude,
                       filename_high=ds_high, filename_low=ds_low)

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
        return self.channels_in_high.intersection(self.channels_in_low)

    @property
    def channels_in_high(self):
        """
        Returns a set of channels in the high resolution data set.
        """
        return set(self.dataset_high.keys())

    @property
    def channels_in_low(self):
        """
        Returns a set of channels in the low resolution data set.
        """
        return set(self.dataset_low.keys())

    @property
    def channels_only_in_high_set(self):
        """
        Returns a set of channels only in the high resolution data set.
        """
        return self.channels_in_high.difference(self.common_channels)

    @property
    def channels_only_in_low_set(self):
        """
        Returns a set of channels only in the low resolution dataset.
        """
        return self.channels_in_low.difference(self.common_channels)

    def add_waveform_to_dataset_high(self, net_sta_comp, filename):
        """
        Add waveforms to the high resolution dataset.

        :param net_sta_comp: A tuple of network id, station id, and the
            component.
        :param filename: The filename of the waveform data.
        """
        # Convert to upper case to be a bit defensive.
        self.dataset_high[tuple([_i.upper() for _i in net_sta_comp])] = \
            filename

    def add_waveform_to_dataset_low(self, net_sta_comp, filename):
        """
        Add waveforms to the low resolution data set.

        :param net_sta_comp: A tuple of network id, station id, and the
            component.
        :param filename: The filename of the waveform data.
        """
        # Convert to upper case to be a bit defensive.
        self.dataset_low[tuple([_i.upper() for _i in net_sta_comp])] = filename

    def set_stations_dataframe(self, df):
        """
        Set the stations dataframe of the data set object.

        :param df: A pandas dataframe object with at least the following
            columns: ``"network"``, ``"station"``, ``"latitude"``,
            ``"longitude"``
        """
        self._stations = df


class WFDiff(object):
    def __init__(self, low_res_seismos, high_res_seismos, station_info,
                 t_min, t_max, dt, starttime=None,
                 endtime=None, get_net_sta_comp_fct=None,
                 is_specfem_ascii=False):
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
        :param is_specfem_ascii: If `True`, the waveform files will be
            assumed to be the ASCII files from SPECFEM.
        :type is_specfem_ascii: bool
        """
        self.low_res_seismos = low_res_seismos
        self.high_res_seismos = high_res_seismos
        self.get_net_sta_comp_fct = get_net_sta_comp_fct
        self.wf_dataset = WaveformDataSet()
        self.is_specfem_ascii = is_specfem_ascii

        self.starttime = starttime
        self.endtime = endtime

        # Get a list of frequencies to test. Make sure t_max is included.
        assert t_min < t_max
        self.periods = np.arange(t_min, t_max + dt * 0.1, dt)

        if COMM.rank == 0:
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
        COMM.barrier()

    def run(self, misfit_type, threshold, output_directory):
        misfit_fct = MISFIT_MAP[misfit_type]


        def split(container, count):
            """
            Simple and elegant function splitting a container into count
            equal chunks.

            Order is not preserved but for the use case at hand this is
            potentially an advantage as data sitting in the same folder thus
            have a higher at being processed at the same time thus the disc
            head does not have to jump around so much. Of course very
            architecture dependent.
            """
            return [container[_i::count] for _i in range(count)]

        if os.path.exists(output_directory):
            raise ValueError("Directory '%s' already exists." %
                             output_directory)

        COMM.barrier()

        if COMM.rank == 0:
            os.makedirs(output_directory)

            jobs = [_i for _i in self.wf_dataset]
            total_length = len(jobs)
            # Collect all jobs on rank 0 and distribute.
            jobs = split(jobs, COMM.size)

            logger.info("Distributing %i jobs across %i cores." % (
                total_length, COMM.size))
        else:
            jobs = None

        jobs = COMM.scatter(jobs, root=0)
        results = []

        for _i, job in enumerate(jobs):
            if self.is_specfem_ascii:
                tr_high = read_specfem_ascii_waveform_file(
                    job.filename_high, network=job.network,
                    station=job.station, channel=job.component)[0]
                tr_low = read_specfem_ascii_waveform_file(
                    job.filename_low, network=job.network, station=job.station,
                    channel=job.component)[0]
            else:
                tr_high = obspy.read(job.filename_high)[0]
                tr_low = obspy.read(job.filename_low)[0]

            misfits.preprocess_traces(tr_high, tr_low)

            tr_high.trim(starttime=tr_high.stats.starttime + self.starttime,
                         endtime=tr_high.stats.starttime + self.endtime)
            tr_low.trim(starttime=tr_low.stats.starttime + self.starttime,
                         endtime=tr_low.stats.starttime + self.endtime)

            # Taper to stabilize the filter.
            tr_high.detrend("demean")
            tr_high.detrend("linear")
            tr_low.detrend("demean")
            tr_low.detrend("linear")
            tr_high.taper(type="hann", max_percentage=0.03)
            tr_low.taper(type="hann", max_percentage=0.03)
            tr_high.differentiate()
            tr_low.differentiate()

            # Calculate the misfit for a number of periods.
            periods = []
            misfit_values = []

            import time
            a = time.time()

            # import matplotlib.pylab as plt
            # plt.figure()
            # count = len(self.periods) + 1
            # plt.subplot(count, 1, 1)
            # plt.plot(tr_high.data, color="red")
            # plt.plot(tr_low.data, color="blue")

            for _i, period in enumerate(self.periods):
                l_tr = tr_low.copy()
                h_tr = tr_high.copy()
                l_tr.filter("lowpass", freq=1.0 / period, corners=3)
                h_tr.filter("lowpass", freq=1.0 / period, corners=3)


                misfit = misfit_fct(l_tr, h_tr)
                periods.append(period)
                misfit_values.append(misfit)

                # plt.subplot(count, 1, _i + 2)
                # plt.plot(h_tr.data, color="red")
                # plt.plot(l_tr.data, color="blue")
                # plt.title("Period %.1f, misfit: %g" % (period, misfit))

            # plt.show()

            periods = np.array(periods)
            misfit_values = np.array(misfit_values)

            # We go from low to high periods. Thus we want to find the first
            # value below the chosen threshold.
            this_threshold = misfit_values.min() + (misfit_values.ptp() *
                                                    threshold)
            idx = np.argmax(misfit_values[::-1] > this_threshold) - 1

            value = periods[-idx]

            import matplotlib.pylab as plt
            plt.style.use("ggplot")
            plt.close()
            plt.figure()
            plt.plot(periods, np.log10(misfit_values))
            plt.xlabel("Minimum period [s]")
            plt.ylabel("log10(wdiff)")
            plt.scatter([periods[-idx]], [misfit_values[-idx]], marker="o")
            plt.xlim(periods[0], periods[-1])
            # plt.hlines(-1, periods[0], periods[-1])
            plt.savefig(os.path.join(
                output_directory,
                "%s_%s_%s.png" % (job.network, job.station,job.component)))


            r = {
                "network": job.network,
                "station": job.station,
                "component": job.component,
                "latitude": job.latitude,
                "longitude": job.longitude,
                "periods": self.periods,
                "misfit_values": misfit_values,
                "misfit_type": misfit_type
            }
            results.append(r)

            b = time.time()
            print("Time taken: ", b - a)

            if COMM.rank == 0:
                print("Approximately %i of %i items have been processed." % (
                    min((_i + 1) * MPI.COMM_WORLD.size, total_length),
                    total_length))

        gathered_results = MPI.COMM_WORLD.gather(results, root=0)

        results = Results()

        if COMM.rank == 0:
            for _i in gathered_results:
                for _j in _i:
                    results.add_result(_j)

            results.write(os.path.join(output_directory, "results.json"))



    def _find_waveform_files(self):
        """
        Finds the waveform files for the low and high resolution seismograms.
        """
        low_res = glob.glob(self.low_res_seismos)
        high_res = glob.glob(self.high_res_seismos)

        for filename in low_res:
            self.wf_dataset.add_waveform_to_dataset_low(
                self.get_net_sta_comp_fct(os.path.basename(filename)),
                filename)

        for filename in high_res:
            self.wf_dataset.add_waveform_to_dataset_high(
                self.get_net_sta_comp_fct(os.path.basename(filename)),
                filename)


        c_chan = self.wf_dataset.common_channels
        a_chan = self.wf_dataset.channels_only_in_high_set
        b_chan = self.wf_dataset.channels_only_in_low_set

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
