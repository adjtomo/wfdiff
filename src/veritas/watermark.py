#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Watermark of the current system to increase reproducibility and provenance.

:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2015
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

from multiprocessing import cpu_count
from pkg_resources import get_distribution
import platform
from socket import gethostname
from time import strftime

# Dependencies.
modules = ["numpy", "scipy", "matplotlib", "obspy", "lxml", "future", "mpi4py",
           "basemap", "pandas"]


def get_watermark():
    watermark = {
        "python_implemenation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "python_compiler": platform.python_compiler(),
        "platform_system": platform.system(),
        "platform_release": platform.release(),
        "platform_version": platform.version(),
        "platform_machine": platform.machine(),
        "platform_processor": platform.processor(),
        "platform_processor_count": cpu_count(),
        "platform_architecture": platform.architecture()[0],
        "platform_hostname": gethostname(),
        "date": strftime('%d/%m/%Y'),
        "time": strftime('%H:%M:%S'),
        "timezone": strftime('%Z')}

    watermark["module_versions"] = {
        module: get_distribution(module).version for module in modules}

    return watermark
