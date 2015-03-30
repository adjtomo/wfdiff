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

plt.style.use("ggplot")


def plot_misfit_curves(items, logarithmic, component, pretty_misfit_name,
                       filename):
    plt.close()
    for item in items:
        if logarithmic:
            plt.semilogy(item["periods"], item["misfit_values"])
        else:
            plt.plot(item["periods"], item["misfit_values"])
        plt.title("%s misfit curves for component %s" % (
            pretty_misfit_name, component))
        plt.xlabel("Lowpass Period [s]")
        plt.ylabel("%s" % pretty_misfit_name)
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
