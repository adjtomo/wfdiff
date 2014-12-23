#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2014
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

import logging


class PysatsiError(Exception):
    """
    Base class for all Pysatsi exceptions. Will probably be used for all
    exceptions to not overcomplicate things as the whole package is pretty
    small.
    """
    pass


class PysatsiWarning(UserWarning):
    """
    Base class for all Pysatsi warnings.
    """
    pass


# Import pysatsi things beneath the exception and warning definitions to avoid
# circular imports.
from .satsi import read_fault_plane_solutions  # NOQA


__version__ = "0.0.0"


# Setup the logger.
logger = logging.getLogger("pysatsi")
logger.setLevel(logging.WARNING)
# Prevent propagating to higher loggers.
logger.propagate = 0
# Console log handler.
ch = logging.StreamHandler()
# Add formatter
FORMAT = "[%(asctime)s] - %(name)s - %(levelname)s: %(message)s"
formatter = logging.Formatter(FORMAT)
ch.setFormatter(formatter)
logger.addHandler(ch)
