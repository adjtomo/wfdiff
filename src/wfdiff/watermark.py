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
           "cartopy", "pandas"]


def get_watermark():
    # Unicode all the way for py2 & py3 compatibility. The future package
    # overrides str() on Python 2.
    watermark = {
        "python_implemenation": str(platform.python_implementation()),
        "python_version": str(platform.python_version()),
        "python_compiler": str(platform.python_compiler()),
        "platform_system": str(platform.system()),
        "platform_release": str(platform.release()),
        "platform_version": str(platform.version()),
        "platform_machine": str(platform.machine()),
        "platform_processor": str(platform.processor()),
        "platform_processor_count": str(cpu_count()),
        "platform_architecture": str(platform.architecture()[0]),
        "platform_hostname": str(gethostname()),
        "date": str(strftime('%d/%m/%Y')),
        "time": str(strftime('%H:%M:%S')),
        "timezone": str(strftime('%Z'))}

    watermark["module_versions"] = {
        module: str(get_distribution(module).version) for module in modules}

    return watermark
