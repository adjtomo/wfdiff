#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Visualization functions for wfdiff.

:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2014
    Julien Thurin (jthurin@alaska.edu), 2022 - Cartopy implementation
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

import matplotlib.cm
import matplotlib.pylab as plt
import numpy as np
import sys

from obspy.imaging.beachball import beach
from obspy.geodetics import gps2dist_azimuth
from cartopy import crs as ccrs

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
    if event is not None:
        lat_plot = np.append(latitudes, event.origins[0].latitude)
        lon_plot = np.append(longitudes, event.origins[0].longitude)
    else:
        lat_plot = latitudes
        lon_plot = longitudes

    lat_mean = (lat_plot.min() + lat_plot.max())/2
    lon_mean = (lon_plot.min() + lon_plot.max())/2

    m = get_basemap(lon_plot.ptp(), lat_plot.ptp(), lon_mean,
                    lat_mean) 

    x, y, _ = m.projection.transform_points(ccrs.PlateCarree(), np.asanyarray(longitudes), np.asanyarray(latitudes)).T

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
    if event is not None:
        tensor  = event.focal_mechanisms[0].moment_tensor.tensor
        ev_mt = [tensor.m_rr, tensor.m_tt, tensor.m_pp,
                 tensor.m_rt, tensor.m_rp, tensor.m_tp]
        # ex, ey = m(event.origins[0].longitude, event.origins[0].latitude)
        ex, ey, _ = m.projection.transform_points(ccrs.PlateCarree(), np.asarray(event.origins[0].longitude),np.asarray(event.origins[0].latitude)).T
        b = beach(ev_mt, xy=(ex, ey), width=int(5000*event.magnitudes[0].mag), 
                  linewidth=0.5, facecolor='deepskyblue')
        ax.add_collection(b)

    cbar = plt.colorbar(data, location="right", pad=0.05)
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

    if event is not None:
        lat_plot = np.append(latitudes, event.origins[0].latitude)
        lon_plot = np.append(longitudes, event.origins[0].longitude)
    else:
        lat_plot = latitudes
        lon_plot = longitudes

    lat_mean = (lat_plot.min() + lat_plot.max())/2
    lon_mean = (lon_plot.min() + lon_plot.max())/2

    m = get_basemap(lon_plot.ptp(), lat_plot.ptp(), lon_mean,
                    lat_mean) 
    x, y, _ = m.projection.transform_points(ccrs.PlateCarree(), np.asanyarray(longitudes), np.asanyarray(latitudes)).T
    plt.close()

    misfit_all= np.asarray(misfit_all)

    nrows = 3
    ncols = 4
    bins_range = (misfit_all.min(), misfit_all.max())

    fig,axes = plt.subplots(ncols=ncols,nrows=nrows,figsize=(3*ncols, 3*nrows),
                      subplot_kw={'projection': m.projection},gridspec_kw = {'wspace':0.2, 'hspace':0.2})
    axes = axes.flatten()
    # Get beachball info
    if event is not None:
        tensor  = event.focal_mechanisms[0].moment_tensor.tensor
        ev_mt = [tensor.m_rr, tensor.m_tt, tensor.m_pp,
                 tensor.m_rt, tensor.m_rp, tensor.m_tp]
        # ex, ey = m(event.origins[0].longitude, event.origins[0].latitude)
        ex, ey, _ = m.projection.transform_points(ccrs.PlateCarree(), np.asarray(event.origins[0].longitude),np.asarray(event.origins[0].latitude)).T

    for i in range(len(items[0]["periods"])):
        ax = axes[i]
        m = get_basemap(lon_plot.ptp(), lat_plot.ptp(), lon_mean,
                        lat_mean, stepsize=4, resolution='110m', ax=ax)
        data = m.scatter(x, y, c=misfit_all[:,i], s=30, vmin=misfit_all.min(),
                         vmax=misfit_all.max(), cmap=cm, alpha=0.9, zorder=10)
        # Add beachball
        if event is not None:
            b = beach(ev_mt, xy=(ex, ey), width=int(8000*event.magnitudes[0].mag), 
                      linewidth=.5, facecolor='deepskyblue')
            ax.add_collection(b)

        ax.set_title(f't > {items[0]["periods"][i]:3.1f} s', fontsize=12)
        ax.set_aspect('equal')

    # Remove unused axes starting from index i
    for j in range(i+1, len(axes)):
        axes[j].remove()


    fig.subplots_adjust(right=0.85, top=0.90, bottom=0.05, left=0.05)
    cbar_ax = fig.add_axes([0.9, 0.1, 0.01, 0.85])
    fig.colorbar(data, cax=cbar_ax)

    fig.suptitle("%s distribution for component %s" % (
            pretty_misfit_name, component))
    
    fig.savefig(filename)

def get_basemap(longitudinal_extent, latitudinal_extent, center_longitude,
                center_latitude, stepsize = None, resolution = None, ax=None):
    """
    Helper function to automatically choose a good map projection with cartopy

    """
    # if stepsize is None:
    #     stepsize = 10
    # if resolution is None:
    #     resolution = '110m'
        
    max_extent = max(longitudinal_extent, latitudinal_extent)

    # Set up the projection
    if max_extent > 180:
        projection = ccrs.Mollweide(central_longitude=center_longitude)
        if ax is None:
            ax = plt.axes(projection=projection)
        else:
            ax.projection = projection
        stepsize = 10
        resolution = '110m'
        ax.set_global()

    elif max_extent > 75:
        projection = ccrs.Orthographic(central_longitude=center_longitude,
                                       central_latitude=center_latitude)
        if ax is None:
            ax = plt.axes(projection=projection)
        else:
            ax.projection = projection
        stepsize = 10
        resolution = '110m'

    else:
        projection = ccrs.LambertAzimuthalEqualArea(central_longitude=center_longitude,
                                                    central_latitude=center_latitude)

        # Calculate map region boundaries
        lat_min = center_latitude - latitudinal_extent / 2
        lat_max = center_latitude + latitudinal_extent / 2
        lon_min = center_longitude - longitudinal_extent / 2
        lon_max = center_longitude + longitudinal_extent / 2
        # Define best resolution and stepsize for the map given region size
        
        if longitudinal_extent > 60:
            if resolution == None:
                resolution = '110m'
            if stepsize == None:
                stepsize = 10
        elif longitudinal_extent > 30:
            if resolution == None:
               resolution = '50m'
            if stepsize == None:
                stepsize = 5
        elif longitudinal_extent > 10:
            if resolution == None:
                resolution = '50m'
            if stepsize == None:
                stepsize = 5
        elif longitudinal_extent > 5:
            if resolution == None:
                resolution = '50m'
            if stepsize == None:
                stepsize = 2
        elif longitudinal_extent > 2:
            if resolution == None:
                resolution = '50m'
            if stepsize == None:
                stepsize = 1
        elif longitudinal_extent > 1:
            if resolution == None:
                resolution = '50m'
            if stepsize == None:
                stepsize = 0.5
        
        # Change map dimensions from degree to meters
        lon_min, lat_min = projection.transform_point(lon_min, lat_min,
                                                    ccrs.PlateCarree())
        lon_max, lat_max = projection.transform_point(lon_max, lat_max,
                                                    ccrs.PlateCarree())
        # add little extra margin around the map
        map_margin = 100000
        lat_min -= map_margin
        lat_max += map_margin
        lon_min -= map_margin
        lon_max += map_margin

        # Set up the map
        if ax is None:
            ax = plt.axes(projection=projection)
        else:
            ax.projection = projection
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=projection)
    
    _plot_feature(ax, stepsize, resolution)
    
    return ax

def _plot_feature(ax, stepsize, resolution,**kwargs):
    """
    Helper function to plot cartopy features

    """
    import matplotlib.pyplot as plt
    import cartopy.feature as cfeature
    # Add grey ocean features
    ax.add_feature(cfeature.OCEAN.with_scale(resolution), facecolor='0.9', edgecolor='0.9')
    # Add grey continents features  
    ax.add_feature(cfeature.LAND.with_scale(resolution),facecolor='0.8', edgecolor='0.8')
    # Add light blue lakes features
    ax.add_feature(cfeature.LAKES.with_scale(resolution), facecolor='0.9', edgecolor='0.9')

    ax.add_feature(cfeature.COASTLINE.with_scale(resolution), linewidth=0.3)
    ax.add_feature(cfeature.BORDERS.with_scale(resolution), linewidth=0.3)
    ax.add_feature(cfeature.RIVERS.with_scale(resolution), linewidth=0.3, edgecolor='0.9')
    if ax.projection.srs.split()[0].split('=')[1] == 'laea':
        gls = ax.gridlines(xlocs=np.arange(-180, 180, stepsize),
                        ylocs=np.arange(-90, 90, stepsize), draw_labels=True, rotate_labels=False, x_inline=False, y_inline=False, **kwargs)
        gls.top_labels = False
        gls.right_labels = False
        gls.xlabel_style = {'size': 8, 'color': 'black'}
        gls.ylabel_style = {'size': 8, 'color': 'black'}
    else:
        gls = ax.gridlines(xlocs=np.arange(-180, 180, stepsize),
                ylocs=np.arange(-90, 90, stepsize), draw_labels=False, **kwargs)
    