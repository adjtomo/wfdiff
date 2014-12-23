#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Visualization functions for pysatsi.

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
# Importing seaborn will result in much prettier plots.
import seaborn as sns  # NOQA


def plot_tradeoff_curve(tradeoff_curve, show=True, filename=None):
    plt.figure(figsize=(12, 5))
    plt.plot(tradeoff_curve["data_misfit"], tradeoff_curve["model_length"],
             lw=0.5, ls="--", marker="o")
    plt.xlabel("Data misfit")
    plt.ylabel("Model length")
    plt.title("Trade-off curve (Selected damping value: %.1f)" %
              tradeoff_curve["selected_damping_value"])
    # Labels.
    for damping, misfit, m_length in zip(
            tradeoff_curve["damping_candidates"],
            tradeoff_curve["data_misfit"],
            tradeoff_curve["model_length"]):
        plt.annotate(s="%.1f" % damping, xy=(misfit, m_length),
                     textcoords="offset points", xytext=(3, 3))

    # Plot the selected damping value.
    idx = tradeoff_curve["selected_damping_value_index"]
    plt.scatter([tradeoff_curve["data_misfit"][idx]],
                [tradeoff_curve["model_length"][idx]], color="red", s=150,
                marker="o")

    if show:
        plt.show()
    if filename:
        plt.savefig(filename)
