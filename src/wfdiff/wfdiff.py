#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Main wfdiff interfaces.

:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2015
    Julien Thurin (jthurin@alaska.edu), 2022 - Cartopy implementation
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
# The following two lines are the main reason why the code works in Python 2
# and Python 3.
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

from collections import defaultdict, namedtuple
import copy
import glob
import json
import os
import sys

import obspy

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
from mpi4py import MPI
import numpy as np

from . import logger, misfits, processing, visualization, watermark
from .specfem_helper import *

try:
    import pyasdf     # For reading asdf waveform files
except ImportError:
    pass

plt.style.use("ggplot")

COMM = MPI.COMM_WORLD

Channel = namedtuple("Channel", ["network", "station", "component",
                                 "latitude", "longitude", "filename_high",
                                 "filename_low"])


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


class Results(object):
    def __init__(self):
        self.__misfit_measurements = {}

    @staticmethod
    def load(filename):
        """
        Load class instance from a JSON file.

        :param filename: Filename to read.
        """
        with open(filename, "r") as fh:
            _results = json.load(fh)
        results = Results()
        del _results["_watermark"]
        results.__misfit_measurements = _results
        return results

    def add_result(self, result):
        name = result["misfit_name"]
        if name not in self.__misfit_measurements:
            self.__misfit_measurements[name] = {
                "misfit_name": result["misfit_name"],
                "misfit_pretty_name": result["misfit_pretty_name"],
                "misfit_logarithmic_plot": result["misfit_logarithmic_plot"],
                "minimizing_misfit": result["minimizing_misfit"],
                "measurements": {}
            }

        self.__misfit_measurements[name]["measurements"][
            "{network}.{station}.{component}".format(**result)] = {
            "network": str(result["network"]),
            "station": str(result["station"]),
            "component": str(result["component"]),
            "latitude": result["latitude"],
            "longitude": result["longitude"],
            "misfit_values": result["misfit_values"],
            "periods": result["periods"]
        }

    def dump(self, filename):
        """
        Serialize the object to disc in form of a JSON file.

        :param filename: The filename to store.
        :type filename: str
        """
        measurements = copy.deepcopy(self.__misfit_measurements)
        measurements["_watermark"] = watermark.get_watermark()

        if sys.version_info[0] < 3:
            mode = "wb"
        else:
            mode = "wt"

        with open(filename, mode) as fh:
            json.dump(measurements, fh, sort_keys=True, indent=4,
                      separators=(u",", u": "))

    @property
    def available_misfits(self):
        """
        Set of all available misfits.
        """
        return set(self.__misfit_measurements.keys())

    def get_available_components_for_misfit(self, misfit_type):
        """
        Set of components that have results.
        """
        return set([
            _i["component"] for _i in
            self.__misfit_measurements[misfit_type]["measurements"].values()])

    def filter(self, misfit, component):
        # Sort to gain some consistency.
        return sorted([
            _i for _i in
            self.__misfit_measurements[misfit]["measurements"].values()
            if _i["component"] == component], key=lambda x: (
            x["network"], x["station"], x["component"]))

    def plot_misfits(self, thresholds, output_directory, output_format = 'pdf'):
        # Make sure all thresholds are available.
        if set(thresholds.keys()) != self.available_misfits:
            raise ValueError("Must specify thresholds for all available "
                             "misfits: '%s'" % self.available_misfits)

        if COMM.rank == 0:
            jobs = []
            for misfit in self.available_misfits:
                for component in \
                        self.get_available_components_for_misfit(misfit):
                    jobs.append((misfit, component))
            jobs = split(jobs, COMM.size)
        else:
            jobs = None

        jobs = COMM.scatter(jobs, root=0)

        for misfit, component in jobs:
            visualization.plot_misfit_curves(
                items=self.filter(misfit, component),
                threshold=thresholds[misfit],
                threshold_is_upper_limit=self.__misfit_measurements[
                    misfit]["minimizing_misfit"],
                logarithmic=self.__misfit_measurements[misfit][
                    "misfit_logarithmic_plot"],
                component=component,
                pretty_misfit_name=self.__misfit_measurements[misfit][
                    "misfit_pretty_name"],
                filename=os.path.join(
                    output_directory,
                    "%s_misfit_curves_%s.%s" % (misfit,
                                                           component, output_format)))

    def plot_misfits_hist(self, thresholds, output_directory, output_format='pdf'):
        # Make sure all thresholds are available.
        if set(thresholds.keys()) != self.available_misfits:
            raise ValueError("Must specify thresholds for all available "
                             "misfits: '%s'" % self.available_misfits)

        if COMM.rank == 0:
            jobs = []
            for misfit in self.available_misfits:
                for component in \
                        self.get_available_components_for_misfit(misfit):
                    jobs.append((misfit, component))
            jobs = split(jobs, COMM.size)
        else:
            jobs = None

        jobs = COMM.scatter(jobs, root=0)

        for misfit, component in jobs:
            visualization.plot_misfit_hist(
                items=self.filter(misfit, component),
                component=component,
                pretty_misfit_name=self.__misfit_measurements[misfit][
                    "misfit_pretty_name"],
                filename=os.path.join(
                    output_directory,
                    "%s_misfit_hist_%s.%s" % (misfit,
                                                           component, output_format)))


    def plot_histograms(self, thresholds, output_directory, output_format='pdf'):
        # Make sure all thresholds are available.
        if set(thresholds.keys()) != self.available_misfits:
            raise ValueError("Must specify thresholds for all available "
                             "misfits: '%s'" % self.available_misfits)

        if COMM.rank == 0:
            jobs = []
            for misfit in self.available_misfits:
                for component in \
                        self.get_available_components_for_misfit(misfit):
                    jobs.append((misfit, component))
            jobs = split(jobs, COMM.size)
        else:
            jobs = None

        jobs = COMM.scatter(jobs, root=0)

        for misfit, component in jobs:
            visualization.plot_histogram(
                items=self.filter(misfit, component),
                threshold=thresholds[misfit],
                threshold_is_upper_limit=self.__misfit_measurements[
                    misfit]["minimizing_misfit"],
                component=component,
                pretty_misfit_name=self.__misfit_measurements[misfit][
                    "misfit_pretty_name"],
                filename=os.path.join(
                    output_directory,
                    "%s_minres_histogram_%s.%s" % (misfit, component, output_format)))

    def plot_maps(self, thresholds, output_directory, event=None, output_format='pdf'):
        # Make sure all thresholds are available.
        if set(thresholds.keys()) != self.available_misfits:
            raise ValueError("Must specify thresholds for all available "
                             "misfits: '%s'" % self.available_misfits)

        if COMM.rank == 0:
            jobs = []
            for misfit in self.available_misfits:
                for component in \
                        self.get_available_components_for_misfit(misfit):
                    jobs.append((misfit, component))
            jobs = split(jobs, COMM.size)
        else:
            jobs = None

        jobs = COMM.scatter(jobs, root=0)

        for misfit, component in jobs:
            visualization.plot_map(
                items=self.filter(misfit, component),
                threshold=thresholds[misfit],
                threshold_is_upper_limit=self.__misfit_measurements[
                    misfit]["minimizing_misfit"],
                component=component,
                pretty_misfit_name=self.__misfit_measurements[misfit][
                    "misfit_pretty_name"],
                filename=os.path.join(
                    output_directory,
                    "%s_minres_peiord_map_%s.%s" % (misfit, component, output_format)),
                event=event)

    def plot_misfit_maps(self, thresholds, output_directory, event=None, output_format='pdf'):
        # Make sure all thresholds are available.
        if set(thresholds.keys()) != self.available_misfits:
            raise ValueError("Must specify thresholds for all available "
                             "misfits: '%s'" % self.available_misfits)

        if COMM.rank == 0:
            jobs = []
            for misfit in self.available_misfits:
                for component in \
                        self.get_available_components_for_misfit(misfit):
                    jobs.append((misfit, component))
            jobs = split(jobs, COMM.size)
        else:
            jobs = None

        jobs = COMM.scatter(jobs, root=0)

        for misfit, component in jobs:
            visualization.plot_misfit_map(
                items=self.filter(misfit, component),
                component=component,
                pretty_misfit_name=self.__misfit_measurements[misfit][
                    "misfit_pretty_name"],
                filename=os.path.join(
                    output_directory,
                    "%s_subplots_map_%s.%s" % (misfit, component, output_format)),
                event=event)

    def plot_all(self, thresholds, output_directory, event_file = None,
                 output_format='pdf'):
        # Make sure all thresholds are available.
        if set(thresholds.keys()) != self.available_misfits:
            raise ValueError("Must specify thresholds for all available "
                             "misfits: '%s'" % self.available_misfits)

        # Read EVENT file if available
        if event_file is None:
            event = None
        else:
            try:
                event = obspy.read_events(event_file)[0]
            except:
                try:
                    asdf_ds = pyasdf.ASDFDataSet(event_file,
                                                 mpi=False, mode="r")
                    event = asdf_ds.events[0]
                except:
                    print('Unrecognized event file - will not plot beachball')
                    event = None

        # Plot results
        self.plot_misfits(thresholds, output_directory, output_format = output_format)
        self.plot_misfits_hist(thresholds, output_directory, output_format = output_format)
        self.plot_histograms(thresholds, output_directory, output_format = output_format)
        self.plot_maps(thresholds, output_directory, event=event, output_format = output_format)
        self.plot_misfit_maps(thresholds, output_directory, event=event, output_format = output_format)


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
    def __init__(self, low_res_seismos, high_res_seismos,
                 stations_file, event_file,
                 t_min, t_max, dt,
                 data_units, desired_analysis_units,
                 rotate_RTZ = False,
                 starttime=None, endtime=None,
                 new_specfem_name=True,
                 trace_tags=['low_res','high_res'],
                 asdf_tags=['ngll5','ngll7'],
                 wf_format='asdf'):
        """

        :param low_res_seismos: A UNIX style wildcard pattern to find the
            low resolution seismograms.
        :type low_res_seismos: str
        :param high_res_seismos: A UNIX style wildcard pattern to find the
            high resolution seismograms.
        :type high_res_seismos: str
        :param stations_file: A UNIX style wildcard pattern to find the
            station information.
        :type stations_file: str
        :param t_min: The minimum period to test.
        :type t_min: float
        :param t_max: The maximum period to test.
        :type t_max: float
        :param dt: The delta for the tested frequencies.
        :type dt: float
        :param data_units: The units of the data, one of ``"displacement"``,
            ``"velocity"``, ``"acceleration"``.
        :type data_units: str
        :param desired_analysis_units: The units in which the analysis
            should be performed. One one of ``"displacement"``,
                ``"velocity"``, ``"acceleration"``.
        :type desired_analysis_units: str
        :param new_specfem_name_format: ``True`` if files are in NET.STA.CHAN format
             ``False`` if files are in STA.NET.CHAN (old naming convention)
        :type new_specfem_name_format: boolean
        :param wf_format: If `specfem`, the waveform files will be
            assumed to be the ASCII files from SPECFEM.
            Other options: `asdf`
        :type is_specfem_ascii: str
        """
        self.low_res_seismos = low_res_seismos
        self.high_res_seismos = high_res_seismos
        self.new_specfem_name = new_specfem_name
        self.wf_dataset = WaveformDataSet()
        self.wf_format = wf_format
        self.trace_tags = trace_tags
        self.asdf_tags = asdf_tags
        self.bundle_jobs_by_channels = True
        self.stations_file = stations_file
        self.rotate_RTZ = rotate_RTZ

        # Read asdf data
        if self.wf_format == 'asdf':
            self.asdf_low = pyasdf.ASDFDataSet(self.low_res_seismos,
                                               mpi=False, mode="r")
            self.asdf_high = pyasdf.ASDFDataSet(self.high_res_seismos,
                                                mpi=False, mode="r")
            # Make sure its the same event in
            assert self.asdf_high.events[0] == self.asdf_low.events[0]
            self.event = self.asdf_high.events[0]
        elif self.wf_format == 'specfem':
            self.event = obspy.read_events(event_file)[0]

        # Check units here.
        acceptable_units = ["displacement", "velocity", "acceleration"]
        data_units = data_units.lower()
        desired_analysis_units = desired_analysis_units.lower()
        if data_units not in acceptable_units:
            raise ValueError(
                "Data unit '%s' is not allowed. Allowed units: "
                "'displacement', 'velocity', or 'acceleration'.")
        if desired_analysis_units not in acceptable_units:
            raise ValueError(
                "Data unit '%s' is not allowed. Allowed units: "
                "'displacement', 'velocity', or 'acceleration'.")
        self.data_units = data_units
        self.desired_analysis_units = desired_analysis_units

        self.starttime = starttime
        self.endtime = endtime

        # Get a list of frequencies to test. Make sure t_max is included.
        assert t_min < t_max
        self.periods = np.arange(t_min, t_max + dt * 0.1, dt)

        if COMM.rank == 0:
            self._find_waveform_files()

            # Get stations dataframe
            if self.wf_format == 'asdf':
                self.stations = get_stations_from_asdf(self.asdf_low)
            else:
                self.stations = read_specfem_stations_file(self.stations_file)

            self.wf_dataset.set_stations_dataframe(self.stations)

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

    def run(self, misfit_types, output_directory,
            save_debug_plots=False, output_format='pdf'):

        misfit_functions = {}

        if self.wf_format == 'asdf':
            self.stations = get_stations_from_asdf(self.asdf_low)
        else:
            self.stations = read_specfem_stations_file(self.stations_file)

        # Check if all the misfit types also have corresponding functions.
        for m_type in misfit_types:
            try:
                fct = getattr(misfits, m_type)
            except AttributeError:
                raise ValueError("Misfit '%s' not known. Known misfits: %s" %
                                 (m_type, ", ".join(["'%s'" % _i for _i in
                                                     misfits.__all__])))
            misfit_functions[m_type] = fct

        if os.path.exists(output_directory):
            raise ValueError("Directory '%s' already exists." %
                             output_directory)

        COMM.barrier()

        debug_folder = os.path.join(output_directory, "debug_plots")

        # Rank zero figures out what to do and distributes it.
        if COMM.rank == 0:
            os.makedirs(output_directory)

            if save_debug_plots:
                if not os.path.exists(debug_folder):
                    os.makedirs(debug_folder)

            jobs = [_i for _i in self.wf_dataset]

            total_length = len(jobs)
            # Collect all jobs on rank 0 and distribute.
            jobs = split(jobs, COMM.size)

            logger.info("Distributing %i jobs across %i cores." % (
                total_length, COMM.size))
        else:
            jobs = None

        jobs = COMM.scatter(jobs, root=0)

        # Each rank now calculates things. No load balancing but as data
        # traces very likely all have the same length, this is not needed.
        results = []
        for jobnum, job in enumerate(jobs):
            # Read the waveform traces. Fork depending on SPECFEM ASCII
            # output or not. If not data is read with ObsPy which should be
            # able to read almost anything.
            if self.wf_format == 'specfem':
                st_high = read_specfem_files(job.filename_high)
                st_low = read_specfem_files(job.filename_low)
            elif self.wf_format == 'asdf':
                st_high = self.asdf_high.waveforms[ \
                    job.network + '_' + job.station][self.asdf_tags[1]]
                st_low = self.asdf_low.waveforms[ \
                    job.network + '_' + job.station][self.asdf_tags[0]]
            else:
                st_high = obspy.read(job.filename_high)[0]
                st_low = obspy.read(job.filename_low)[0]

            # Rotate
            if self.rotate_RTZ is True:
                # Compute back_azimuth using event info and add it to trace.stats
                # And then rotate
                st_high = add_event_station_info(st_high, self.event, self.stations)
                st_high.rotate('NE->RT')
                st_low = add_event_station_info(st_low, self.event, self.stations)
                st_low.rotate('NE->RT')

            # Sort to make sure they are in same order
            st_high.sort()
            st_low.sort()
            # Loop over each trace of stream
            for i in range(len(st_high)):
                tr_high = st_high[i]
                tr_low = st_low[i]

                # Make sure you are comparing same components
                assert tr_high.stats.channel[-1] == tr_low.stats.channel[-1]
                component = tr_high.stats.channel[-1]

                # Preprocess data. Afterwards they will be sampled at exactly
                # the same points in time and have the desired units.
                processing.preprocess_traces(
                    tr_high, tr_low,
                    starttime=self.starttime,
                    endtime=self.endtime,
                    data_units=self.data_units,
                    desired_units=self.desired_analysis_units)

                # Calculate each misfit for a range of periods.
                collected_misfits = defaultdict(list)

                stnm_tag =  job.network + '.' + job.station + '.' + component
                if save_debug_plots:
                    plt.close()
                    plt.style.use("ggplot")
                    plt.figure(figsize=(10, len(self.periods)))
                    plt.subplot(len(self.periods) + 1, 1, 1)
                    plt.subplots_adjust(left=0.0, right=1.0, top=1.0,
                                        bottom=0.0, wspace=0.0, hspace=0.0)
                    plt.plot(tr_high.data, color="red", label=self.trace_tags[1])
                    plt.plot(tr_low.data, color="blue", label=self.trace_tags[0])
                    plt.legend(fontsize=8, loc=3)
                        # x-axis tick-marks and labels
                    xtick_len = 20
                    locs = np.arange(self.starttime/tr_low.stats.delta,
                                     self.endtime/tr_low.stats.delta + 1,
                                     xtick_len/tr_low.stats.delta)
                    labels = map(str, locs*tr_low.stats.delta)
                    plt.xticks(locs)
                    plt.yticks([])

                    ax = plt.gca()
                    ax.text(0.02, 0.95, stnm_tag, transform=ax.transAxes,
                            fontdict=dict(fontsize=12, ha='left', va='top'),
                            bbox=dict(boxstyle="round", fc="w", alpha=0.8))
                    ax.text(0.02, 0.65, "raw", transform=ax.transAxes,
                            fontdict=dict(fontsize="small", ha='left', va='top'),
                            bbox=dict(boxstyle="round", fc="w", alpha=0.8))

                # Loop over various periods
                for _i, period in enumerate(self.periods):
                    l_tr = tr_low.copy()
                    h_tr = tr_high.copy()
                    l_tr.filter("lowpass", freq=1.0 / period, corners=3)
                    h_tr.filter("lowpass", freq=1.0 / period, corners=3)

                    this_misfits = {}

                    # Calculate each desired misfit.
                    for name, fct in misfit_functions.items():
                        # Potentially returns multiple measures, e.g. CCs return
                        # coefficient and time shift.
                        mfs = fct(l_tr, h_tr)
                        if isinstance(mfs, dict):
                            mfs = [mfs]
                        for mf in mfs:
                            collected_misfits[mf["name"]].append(mf)
                            this_misfits[mf["name"]] = mf["value"]

                    if save_debug_plots:
                        plt.subplot(len(self.periods) + 1, 1, _i + 2)
                        plt.plot(h_tr.data, color="red")
                        plt.plot(l_tr.data, color="blue")
                        # x-axis label
                        if _i ==1:
                            plt.xticks(locs, labels)
                            plt.tick_params(labelbottom='off',labeltop='on')
                        else:
                            plt.xticks(locs)
                        plt.yticks([])
                        ax = plt.gca()
                        ax.text(0.02, 0.95, "Lowpass period: %.1f" % period,
                                transform=ax.transAxes,
                                fontdict=dict(fontsize="small", ha='left',
                                              va='top'),
                                bbox=dict(boxstyle="round", fc="w", alpha=0.8))
                        keys = sorted(this_misfits.keys())
                        txt = ["%s: %g" % (key, this_misfits[key]) for key in keys]
                        ax.text(0.98, 0.95, "\n".join(txt),
                                ha="right", transform=ax.transAxes,
                                fontdict=dict(fontsize="small", ha='right',
                                              va='top'),
                                bbox=dict(boxstyle="round", fc="w", alpha=0.8))

                if save_debug_plots:
                    filename = os.path.join(
                        debug_folder,
                        "%s_%s_%s.%s" % (job.network, job.station,
                                         component, output_format))
                    plt.savefig(filename)

                # Now assemble a frequency dependent misfit measurement for each
                # final misfit type.
                for key, value in collected_misfits.items():
                    r = {
                        "network": job.network,
                        "station": job.station,
                        "component": component,
                        "latitude": job.latitude,
                        "longitude": job.longitude,
                        "periods": list(self.periods),
                        "misfit_values": [_i["value"] for _i in value],
                        "misfit_name": value[0]["name"],
                        "misfit_pretty_name": value[0]["pretty_name"],
                        "misfit_logarithmic_plot": value[0]["logarithmic_plot"],
                        "minimizing_misfit": value[0]["minimizing_misfit"]
                        }
                    results.append(r)

            if COMM.rank == 0:
                logger.info(
                    "Approximately %i of %i stations have been processed." % (
                        min((jobnum + 1) * MPI.COMM_WORLD.size, total_length),
                        total_length))

        # Gather and store results on disc.
        gathered_results = MPI.COMM_WORLD.gather(results, root=0)
        if COMM.rank == 0:
            results = Results()

            for _i in gathered_results:
                for _j in _i:
                    results.add_result(_j)

            results.dump(os.path.join(output_directory, "results.json"))

        results = COMM.bcast(results, root=0)

        return results


    def _find_waveform_files(self):
        """
        Finds the waveform files for the low and high resolution seismograms.
        This is not needed for asdf format files
        """
        low_res = glob.glob(self.low_res_seismos)
        high_res = glob.glob(self.high_res_seismos)

        # For asdf input files
        if self.wf_format == 'asdf':
            for _, st in enumerate(self.asdf_low.waveforms):
                for tr in st[self.asdf_tags[0]]:
                    net, sta, chan = tr.stats.network, tr.stats.station, tr.stats.channel[-1]
                    sta_tag = net + '_' + sta
                    filename = ''
                    if self.bundle_jobs_by_channels is True:
                        chan = ''
                    self.wf_dataset.add_waveform_to_dataset_low([net, sta, chan],
                                                                filename)

            for _, st in enumerate(self.asdf_high.waveforms):
                for tr in st[self.asdf_tags[1]]:
                    net, sta, chan = tr.stats.network, tr.stats.station, tr.stats.channel[-1]
                    sta_tag = net + '_' + sta
                    filename = ''
                    if self.bundle_jobs_by_channels is True:
                        chan = ''
                    self.wf_dataset.add_waveform_to_dataset_high([net, sta, chan],
                                                                filename)


        # For specfem files
        else:
            for filename in low_res:
                net, sta, chan = get_net_sta_comp(os.path.basename(filename),
                                                  new_format=self.new_specfem_name)
                if self.bundle_jobs_by_channels is True:
                    fname = filename[:-8] + '*'  # replace CHA.semd by *
                    chan = ''
                self.wf_dataset.add_waveform_to_dataset_low(
                    [net, sta, chan], fname)

            for filename in high_res:
                net, sta, chan = get_net_sta_comp(os.path.basename(filename),
                                                  new_format=self.new_specfem_name)
                if self.bundle_jobs_by_channels is True:
                    fname = filename[:-8] + '*'  # replace CHA.semd by *
                    chan = ''
                self.wf_dataset.add_waveform_to_dataset_high(
                    [net, sta, chan], fname)

        c_chan = self.wf_dataset.common_channels
        a_chan = self.wf_dataset.channels_only_in_high_set
        b_chan = self.wf_dataset.channels_only_in_low_set

        logger.info("Found %i Stations files that are part of both "
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
