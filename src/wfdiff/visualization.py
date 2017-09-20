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

import matplotlib.cm
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import sys

from obspy.imaging.beachball import beach
from obspy.geodetics import gps2dist_azimuth

from .utils import rightmost_threshold_crossing

try:
    import adjustText     # [unoffical] To prevent overlapping station names on the map 
except ImportError:
    pass

plt.style.use("ggplot")


def plot_misfit_curves(items, threshold, threshold_is_upper_limit,
                       logarithmic, component, pretty_misfit_name, filename):
    plt.close()

    crossing_periods = []
    crossing_values = []
    
    misfit_all = []
    for item in items:
        if logarithmic:
            plt.semilogy(item["periods"], item["misfit_values"])
        else:
            plt.plot(item["periods"], item["misfit_values"], color="blue", 
                     alpha=0.15, lw = 3)

        # Find the threshold.
        point = rightmost_threshold_crossing(
            item["periods"], item["misfit_values"], threshold,
            threshold_is_upper_limit)
        crossing_periods.append(point[0])
        crossing_values.append(point[1])

        misfit_all.append(item['misfit_values'])

    # compute mean and median of misfit for all stations at each filter period
    misfit_all= np.asarray(misfit_all)
    misfit_mean = misfit_all.mean(axis=0)
    misfit_std = misfit_all.std(axis=0)
    misfit_median = np.median(misfit_all, axis=0)
 
    plt.plot(np.asarray(items[0]["periods"]), misfit_mean, color="red", 
             lw = 2, label='mean')
    # Standard deviation doesn't make sense for a non-normal distribution
    #plt.errorbar(np.asarray(items[0]["periods"]), misfit_mean, misfit_std,
    #         lw = 2, zorder=3)
    plt.plot(np.asarray(items[0]["periods"]), misfit_median, color="Chartreuse", 
             lw = 2, label='median', linestyle="--")

    plt.title("%s misfit curves for component %s" % (
        pretty_misfit_name, component))
    plt.xlabel("Lowpass Period [s]")
    plt.ylabel("%s" % pretty_misfit_name)

    x = items[0]["periods"][0] - 0.5, items[0]["periods"][-1] + 0.5

    plt.hlines(threshold, x[0], x[1],
               linestyle="--", color="0.5")
    plt.scatter(crossing_periods, crossing_values, color="orange", s=10,
                zorder=5, alpha=0.3)
    plt.xlim(*x)

    plt.savefig(filename)

def plot_misfit_hist(items, component, pretty_misfit_name, filename):
    # Plot multiple histograms
    # Histograms of misfit distribution for all stations at each filter period
    plt.close()

    misfit_all = np.array([_i["misfit_values"] for _i in items])

    nrows = 3
    ncols = 4
    bins_range = (misfit_all.min(), misfit_all.max())

    fig = plt.figure(figsize=(3*ncols, 3*nrows))

    for i in range(len(items[0]["periods"])):
        ax = fig.add_subplot(nrows, ncols, i+1)
        #ax=plt.gca()
        ax.hist(misfit_all[:,i], range=bins_range, bins=50)
        ax.set_title('t > ' + str(items[0]["periods"][i]) + ' s', fontsize=12)
        ax.set_ylim([0, 80])
        ax.set_ylabel("Nstn")
        ax.set_xlabel("%s" % pretty_misfit_name)

    fig.suptitle("%s distribution for component %s" % (
        pretty_misfit_name, component))
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    fig.savefig(filename)

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


def plot_map(items, threshold, threshold_is_upper_limit,
             component, pretty_misfit_name, filename, event=None):
    # Choose red to green colormap.
    cm = matplotlib.cm.RdYlGn_r

    longitudes = np.array([_i["longitude"] for _i in items])
    latitudes = np.array([_i["latitude"] for _i in items])

    resolvable_periods = []
    station_array = []
    period_range = items[0]["periods"]
 
    for item in items:
        # Find the threshold.
        point = rightmost_threshold_crossing(
            item["periods"], item["misfit_values"], threshold,
            threshold_is_upper_limit)
        resolvable_periods.append(point[0])
        station_array.append(item["network"] + '_' + item["station"])

    resolvable_periods = np.array(resolvable_periods)

    plt.close()
    lat_plot = np.append(latitudes, event.origins[0].latitude)
    lon_plot = np.append(longitudes, event.origins[0].longitude)
    lat_mean = (lat_plot.min() + lat_plot.max())/2
    lon_mean = (lon_plot.min() + lon_plot.max())/2
    m = get_basemap(lon_plot.ptp(), lat_plot.ptp(), lon_mean,
                    lat_mean) 

    x, y = m(longitudes, latitudes)

    data = m.scatter(x, y, c=resolvable_periods, s=30, vmin=period_range[0],
                     vmax=period_range[-1], cmap=cm, alpha=0.9, zorder=10)
    # add station label
    texts = []
    
    for stnm, xi, yi in zip(station_array, x, y):
        texts.append(plt.text(xi, yi, stnm,fontsize=5)) 
    if 'adjustText' in sys.modules:   # require adjust_Text module
        adjustText.adjust_text(texts, force_points=1, force_text=1, expand_points=(1,1), 
                    expand_text=(1,1), arrowprops=dict(arrowstyle="<-", color='black', 
                                                       alpha=0.5, lw=0.5)) 
    ax = plt.gca()

    # plot beachball
    tensor  = event.focal_mechanisms[0].moment_tensor.tensor
    ev_mt = [tensor.m_rr, tensor.m_tt, tensor.m_pp,
             tensor.m_rt, tensor.m_rp, tensor.m_tp]
    ex, ey = m(event.origins[0].longitude, event.origins[0].latitude)
    b = beach(ev_mt, xy=(ex, ey), width=int(5000*event.magnitudes[0].mag), 
              linewidth=0.5, facecolor='deepskyblue')
    ax.add_collection(b)

    cbar = m.colorbar(data, location="right", pad="15%")
    cbar.set_label("Minimum Resolvable Period [s]")

    plt.title("%s minimum resolvable period for component %s" % (
              pretty_misfit_name, component), fontsize="small")
    plt.savefig(filename)


def plot_misfit_map(items, component, pretty_misfit_name, filename, event=None):

    # Choose red to green colormap.
    cm = matplotlib.cm.RdYlGn_r

    longitudes = np.array([_i["longitude"] for _i in items])
    latitudes = np.array([_i["latitude"] for _i in items])
    
    misfit_all = np.array([_i["misfit_values"] for _i in items])

    plt.close()
    lat_plot = np.append(latitudes, event.origins[0].latitude)
    lon_plot = np.append(longitudes, event.origins[0].longitude)
    lat_mean = (lat_plot.min() + lat_plot.max())/2
    lon_mean = (lon_plot.min() + lon_plot.max())/2
    m = get_basemap(lon_plot.ptp(), lat_plot.ptp(), lon_mean,
                    lat_mean) 
    x, y = m(longitudes, latitudes)

    misfit_all= np.asarray(misfit_all)

    nrows = 3
    ncols = 4
    bins_range = (misfit_all.min(), misfit_all.max())

    fig = plt.figure(figsize=(3*ncols, 3*nrows))

    # Get beachball info
    tensor  = event.focal_mechanisms[0].moment_tensor.tensor
    ev_mt = [tensor.m_rr, tensor.m_tt, tensor.m_pp,
             tensor.m_rt, tensor.m_rp, tensor.m_tp]
    ex, ey = m(event.origins[0].longitude, event.origins[0].latitude)

    for i in range(len(items[0]["periods"])):
        ax = fig.add_subplot(nrows, ncols, i+1)
        m = get_basemap(lon_plot.ptp(), lat_plot.ptp(), lon_mean,
                        lat_mean, stepsize=4, resolution='c')
        data = m.scatter(x, y, c=misfit_all[:,i], s=30, vmin=misfit_all.min(),
                         vmax=misfit_all.max(), cmap=cm, alpha=0.9, zorder=10)
        # Add beachball
        b = beach(ev_mt, xy=(ex, ey), width=int(8000*event.magnitudes[0].mag), 
              linewidth=.5, facecolor='deepskyblue')
        ax.add_collection(b)
        # Add colorbar for the last subplot
        #if i == (len(items[0]["periods"]) -1):
        #    divider = make_axes_locatable(ax)
        #    cax = divider.append_axes("right", size="5%", pad=0.05)
        #    cbar = m.colorbar(data, location="right", pad="15%", shrink=0.75)
        #    cbar.set_label(pretty_misfit_name)

        ax.set_title('t > ' + str(items[0]["periods"][i]) + ' s', fontsize=12)
        ax.add_collection(b)

    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.9, 0.1, 0.01, 0.85])
    fig.colorbar(data, cax=cbar_ax)
    fig.tight_layout(rect=[0.08, 0.03, 0.9, 0.9])

    fig.suptitle("%s distribution for component %s" % (
            pretty_misfit_name, component))
    filename = 'output_NGLL_test_cook_basin3_vsmin_1000_ismooth1/maps_' + component +'.pdf'
    
    fig.savefig(filename)


def get_basemap(longitudinal_extent, latitudinal_extent, center_longitude,
                center_latitude, stepsize = None, resolution = None, ax=None):
    """
    Helper function attempting to automatically choose a good map projection.

    Likely not perfect.
    """
    if ax is None:
        ax = plt.gca()

    max_extent = max(longitudinal_extent, latitudinal_extent)

    # Use a global plot for very large domains.
    if max_extent >= 180.0:
        m = Basemap(projection='moll', lon_0=0, resolution="c", ax=ax)
        stepsize = 45.0
    # Orthographic projection for 75.0 <= extent < 180.0
    elif max_extent >= 75.0:
        m = Basemap(projection="ortho", lon_0=center_longitude,
                    lat_0=center_latitude, resolution="c", ax=ax)
        stepsize = 10.0
    # Lambert azimuthal equal area projection. Equal area projections
    # are useful for interpreting features and this particular one also
    # does not distort features a lot on regional scales.
    else:
        # Calculate map region 
        lat_min = center_latitude - (latitudinal_extent / 2.0)
        lat_max = center_latitude + (latitudinal_extent / 2.0)
        lon_min = center_longitude - (longitudinal_extent / 2.0)
        lon_max = center_longitude + (longitudinal_extent / 2.0)
      
        # Try to pick suitable tick-marks increment and resolution
        # on the basis of size of map region
        if stepsize is None or resolution is None:
            if longitudinal_extent > 50.0:
                stepsize = 10.0
                resolution = "i"
            elif 20.0 < longitudinal_extent <= 50.0:
                stepsize = 5.0
                resolution = "i"
            elif 5.0 < longitudinal_extent <= 20.0:
                stepsize = 2.0
                resolution = "i"
            elif 2.0 < longitudinal_extent < 5.0:
                stepsize = 1.0
                resolution = "h"
            else:
                stepsize = 0.5
                resolution = "h"
        
        # Change map dimensions from degree to meters
        width, _, _  = gps2dist_azimuth(center_latitude, lon_min, center_latitude, lon_max) 
        height, _, _ = gps2dist_azimuth(lat_min, center_longitude, lat_max, center_longitude) 
        # add little extra margin around the map
        map_margin = 100000
        width += map_margin
        height += map_margin
       
        m = Basemap(projection='laea', resolution=resolution, width=width,
                    height=height, lat_ts=center_latitude, lat_0=center_latitude,
                    lon_0=center_longitude, ax=ax)

    _plot_features(m, stepsize)
    return m


def _plot_features(map_object, stepsize):
    """
    Helper function aiding in consistent plot styling.
    """
    import matplotlib.pyplot as plt

    map_object.drawmapboundary(fill_color='#bbbbbb')
    map_object.fillcontinents(color='white', lake_color='#cccccc', zorder=0)
    plt.gcf().patch.set_alpha(0.0)

    # Style for parallels and meridians.
    LINESTYLE = {
        "linewidth": 0.5,
        "dashes": [],
        "color": "#999999"}

    # Parallels.
    if map_object.projection in ["moll", "laea"]:
        label = True
    else:
        label = False
    parallels = np.arange(-90.0, 90.0, stepsize)
    map_object.drawparallels(parallels, labels=[label, False, False, False],
                             zorder=200, **LINESTYLE)
    # Meridians.
    if map_object.projection in ["laea"]:
        label = True
    else:
        label = False
    meridians = np.arange(0.0, 360.0, stepsize)
    map_object.drawmeridians(
        meridians, labels=[False, False, False, label], zorder=200,
        **LINESTYLE)

    map_object.drawcoastlines(color="#444444", linewidth=0.7)
    map_object.drawcountries(linewidth=0.4, color="#777777")
    if stepsize < 10.0:
        map_object.drawstates(linewidth=0.2, color="#999999")
        map_object.drawrivers(linewidth=0.1, color="#888888",
                              linestyle="solid")
