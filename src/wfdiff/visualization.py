#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Visualization functions for wfdiff.

:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2014
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

import matplotlib.pylab as plt

from .utils import rightmost_threshold_crossing

plt.style.use("ggplot")


def plot_misfit_curves(items, threshold, threshold_is_upper_limit,
                       logarithmic, component, pretty_misfit_name, filename):
    plt.close()

    crossing_periods = []
    crossing_values = []

    for item in items:
        if logarithmic:
            plt.semilogy(item["periods"], item["misfit_values"])
        else:
            plt.plot(item["periods"], item["misfit_values"])

        # Find the threshold.
        point = rightmost_threshold_crossing(
            item["periods"], item["misfit_values"], threshold,
            threshold_is_upper_limit)
        crossing_periods.append(point[0])
        crossing_values.append(point[1])

    plt.title("%s misfit curves for component %s" % (
        pretty_misfit_name, component))
    plt.xlabel("Lowpass Period [s]")
    plt.ylabel("%s" % pretty_misfit_name)

    x = items[0]["periods"][0] - 0.5, items[0]["periods"][-1] + 0.5

    plt.hlines(threshold, x[0], x[1],
               linestyle="--", color="0.5")
    plt.scatter(crossing_periods, crossing_values, color="0.2", s=10,
                zorder=5)
    plt.xlim(*x)

    plt.savefig(filename)


def plot_histogram(items, threshold, threshold_is_upper_limit,
                   component, pretty_misfit_name, filename):
    plt.close()
    plt.figure(figsize=(12, 4))

    crossing_periods = []

    for item in items:
        # Find the threshold.
        point = rightmost_threshold_crossing(
            item["periods"], item["misfit_values"], threshold,
            threshold_is_upper_limit)
        crossing_periods.append(point[0])

    plt.hist(crossing_periods, bins=20)

    plt.title("Minimum resolvable period for %s and component %s" % (
        pretty_misfit_name, component))
    plt.xlabel("Lowpass Period [s]")
    plt.ylabel("Count")
    plt.tight_layout()

    plt.savefig(filename)


def plot_misfit_map(items, component, pretty_misfit_name, filename):
    longitudes = [_i["longitude"] for _i in items]
    latitudes = [_i["latitudes"] for _i in items]

    plt.close()
    for item in items:
        plt.plot(item["periods"], item["misfit_values"])
        plt.title("%s misfit curves for component %s" % (
            pretty_misfit_name, component))
        plt.xlabel("Lowpass Period [s]")
        plt.ylabel("%s" % pretty_misfit_name)
    plt.savefig(filename)
