#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test suite for the pysatsi visualizations.

Run with pytest.

:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2015
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
import inspect
import matplotlib.pyplot as plt
from matplotlib.testing.compare import compare_images as mpl_compare_images
import numpy as np
import os

from pysatsi.visualization import plot_tradeoff_curve


# Most generic way to get the data folder path.
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(
    inspect.getfile(inspect.currentframe()))), "data")

# Baseline images for the plotting test.
IMAGE_DIR = os.path.join(os.path.dirname(DATA_DIR), "baseline_images")


def reset_matplotlib():
    """
    Reset matplotlib to a common default.
    """
    # Force agg backend.
    plt.switch_backend('agg')
    import locale
    locale.setlocale(locale.LC_ALL, str('en_US.UTF-8'))


reset_matplotlib()


def images_are_identical(image_name, temp_dir, dpi=None):
    """
    Partially copied from ObsPy. Used to check images for equality.
    """
    image_name += os.path.extsep + "png"
    expected = os.path.join(IMAGE_DIR, image_name)
    actual = os.path.join(temp_dir, image_name)

    if dpi:
        plt.savefig(actual, dpi=dpi)
    else:
        plt.savefig(actual)
    plt.close()

    assert os.path.exists(expected)
    assert os.path.exists(actual)

    print(actual)

    # Use a reasonably high tolerance to get around difference with different
    # freetype and possibly agg versions. matplotlib uses a tolerance of 13.
    result = mpl_compare_images(expected, actual, 5, in_decorator=True)
    if result is not None:
        print(result)
    assert result is None


def test_plots_2D_tradeoff_curve(tmpdir):
    """
    Tests plotting the 2D tradeoff curve.
    """
    tc = {
        "damping_candidates": np.array([
            0.4, 0.6, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2., 2.2, 2.4,
            2.6, 2.8, 3., 3.5, 4., 5., 6.]),
        "damping_candidates": np.array(
            [0.4, 0.6, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.6, 1.8,
             2., 2.2, 2.4, 2.6, 2.8, 3., 3.5, 4., 5., 6.]),
        "data_misfit": np.array([
                0.43616322, 0.43616382, 0.43616887, 0.43617663, 0.43619133,
                0.43621632, 0.43625504, 0.43631055, 0.43638502, 0.43659432,
                0.43688097, 0.43723418, 0.43764157, 0.43809272, 0.43857934,
                0.43909394, 0.43962895, 0.44100894, 0.44239568, 0.445234,
                0.44837617]),
        "model_length": np.array([
                0.11355583, 0.11228609, 0.10914431, 0.10666015, 0.10355018,
                0.09987661, 0.0957507, 0.09131161, 0.08670465, 0.07749526,
                0.06888049, 0.06120687, 0.05453153, 0.04876978, 0.04379686,
                0.03949697, 0.03577625, 0.02852976, 0.023449, 0.01680472,
                0.01228102]),
        "selected_damping_value":  2.8,
        "selected_damping_value_index": 15}

    plot_tradeoff_curve(tc)
    images_are_identical("tradeoff_curve", str(tmpdir))
