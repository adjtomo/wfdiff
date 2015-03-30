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

from collections import defaultdict, namedtuple
import copy
import glob
import json
import os

import obspy
import matplotlib.pylab as plt
from mpi4py import MPI
import numpy as np

from . import logger, misfits, processing, visualization, watermark
from .io import read_specfem_stations_file, read_specfem_ascii_waveform_file

plt.style.use("ggplot")

COMM = MPI.COMM_WORLD

Channel = namedtuple("Channel", ["network", "station", "component",
                                 "latitude", "longitude", "filename_high",
                                 "filename_low"])


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
            "network": result["network"],
            "station": result["station"],
            "component": result["component"],
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
        with open(filename, "w") as fh:
            json.dump(measurements, fh, sort_keys=True, indent=4,
                      separators=(",", ": "))

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
        return [
            _i for _i in
            self.__misfit_measurements[misfit]["measurements"].values()
            if _i["component"] == component]

    def plot_misfits(self, output_directory):
        for misfit in self.available_misfits:
            for component in self.get_available_components_for_misfit(misfit):
                visualization.plot_misfit_curves(
                    items=self.filter(misfit, component),
                    logarithmic = self.__misfit_measurements[misfit][
                        "misfit_logarithmic_plot"],
                    component=component,
                    pretty_misfit_name= self.__misfit_measurements[misfit][
                        "misfit_pretty_name"],
                    filename= os.path.join(
                        output_directory,
                        "%s_misfit_curves_component_%s.pdf" % (misfit,
                                                               component)))

    def plot_maps(self, output_directory):
        for misfit in self.available_misfits:
            for component in self.get_available_components_for_misfit(misfit):
                visualization.plot_misfit_map(
                    items=self.filter(misfit, component),
                    component=component,
                    pretty_misfit_name= self.__misfit_measurements[misfit][
                        "misfit_pretty_name"],
                    filename= os.path.join(
                        output_directory,
                        "%s_misfit_map_component_%s.pdf" % (misfit,
                                                            component)))


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
                 t_min, t_max, dt,
                 data_units, desired_analysis_units,
                 starttime=None,
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
        :param data_units: The units of the data, one of ``"displacement"``,
            ``"velocity"``, ``"acceleration"``.
        :type data_units: str
        :param desired_analysis_units: The units in which the analysis
            should be performed. One one of ``"displacement"``,
                ``"velocity"``, ``"acceleration"``.
        :type desired_analysis_units: str
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

    def run(self, misfit_types, output_directory, save_debug_plots=False):
        misfit_functions = {}
        # Check if all the misfit types also have corresponding functions.
        for m_type in misfit_types:
            try:
                fct = getattr(misfits, m_type)
            except AttributeError:
                raise ValueError("Misfit '%s' not known. Known misfits: %s" %
                                 (m_type, ", ".join(["'%s'" % _i for _i in
                                                     misfits.__all__])))
            misfit_functions[m_type] = fct

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

        # Rank zero figures out what to do and distributes it.
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

        # Each rank now calculates things. No load balancing but as data
        # traces very likely all have the same length, this is not needed.
        results = []
        for jobnum, job in enumerate(jobs):
            # Read the waveform traces. Fork depending on SPECFEM ASCII
            # output or not. If not data is read with ObsPy which should be
            # able to read almost anything.
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

                # Calculate each desired misfit.
                for name, fct in misfit_functions.items():
                    # Potentially returns multiple measures, e.g. CCs return
                    # coefficient and time shift.
                    mfs = fct(l_tr, h_tr)
                    if isinstance(mfs, dict):
                        mfs = [mfs]
                    for mf in mfs:
                        collected_misfits[mf["name"]].append(mf)

                # plt.subplot(count, 1, _i + 2)
                # plt.plot(h_tr.data, color="red")
                # plt.plot(l_tr.data, color="blue")
                # plt.title("Period %.1f, misfit: %g" % (period, misfit))

            # plt.show()

            # Now assemble a frequency dependent misfit measurement for each
            # final misfit type.
            for key, value in collected_misfits.items():
                r = {
                    "network": job.network,
                    "station": job.station,
                    "component": job.component,
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
                    "Approximately %i of %i channels have been processed." % (
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
            #results.plot_misfits(output_directory)


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
