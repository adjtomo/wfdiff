#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
I/O Utilities for veritas.

:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2015
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

import pandas


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
