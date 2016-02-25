"""
awot.graph.common
=================

Common graphing routines.

"""

from __future__ import print_function
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.dates import (
    DateFormatter, SecondLocator, MinuteLocator,
    HourLocator, DayLocator)
from matplotlib import ticker as mtic
import matplotlib.cm as cm
from datetime import datetime


########################
#   Global Variables   #
########################
DATE_STRING_FORMAT = ("Date string format: YYYY-MM-DDTHH:MM:SS, "
                      "(e.g. '2014-08-20T12:30:00')")
EARTH_RADIUS = 6371.

###################
#   Map modules   #
###################


def create_basemap(corners=None, proj=None, resolution='l',
                   area_thresh=None, lon_0=None, lat_0=None,
                   meridians=True, parallels=True,
                   lon_spacing=2., lat_spacing=2.,
                   coastlines=True, countries=True, states=False,
                   counties=False,
                   rivers=False, etopo=False, ax=None):
    """
    Create a basemap instance for plotting.

    Parameters
    ----------
    corners : float tuple
        Array of map corners
        [lowerLeftLon, lowerLeftLat, upperRightLon, upperRightLat].
        Default will find min max of arrays.
    proj : str
        Map projection to use.
    resolution : string
        Map resolution, 'l' or low res is default.
    area_thresh : float
        Minimum area threshold for basemap plot instance.
    lon_0 : float
        Longitudinal center of desire map [degrees].
    lat_0 : float
        Latitudinal center of desire map [degrees].
    meridians : bool
        Flag to turn on meridian (lonigitude) line plotting.
    parallels : bool
        Flag to turn on parallels (latitude) line plotting.
    lon_spacing : float
        Spacing for meridians, i.e. longitude [degrees].
    lat_spacing : float
        Spacing for parallels, i.e. latitude [degrees].
    coastlines : bool
        Flag to turn on basemap coastline plotting.
    countries : bool
        Flag to turn on basemap country plotting.
    states : bool
        Flag to turn on basemap state plotting.
    counties : bool
        Flag to turn on basemap county plotting.
    rivers : bool
        Flag to turn on basemap river plotting.
    etopo : bool
        Flag to turn on basemap etopo (topography) plotting.
    ax : Matplotlib axis instance
            Axis to plot. None will use the current axis.
    """
    # parse parameters
    ax = _parse_ax(ax)

    # Check inputs
    if corners is None:
        corners = [-180., -90., 180., 90.]

    if proj is None:
        proj = 'cea'

    # Create a basemap instance
    bm = Basemap(projection=proj, resolution=resolution,
                 area_thresh=area_thresh,
                 llcrnrlon=corners[0], urcrnrlon=corners[2],
                 llcrnrlat=corners[1], urcrnrlat=corners[3],
                 lon_0=lon_0, lat_0=lat_0, ax=ax)

    # Check the customizations for the basemap
    if meridians:
        bm.drawmeridians(np.arange(corners[0], corners[
                         2], lon_spacing), labels=[1, 0, 0, 1])
    if parallels:
        bm.drawparallels(np.arange(corners[1], corners[
                         3], lat_spacing), labels=[1, 0, 0, 1])
    if countries:
        bm.drawcountries()
    if coastlines:
        bm.drawcoastlines()
    if states:
        bm.drawstates()
    if counties:
        bm.drawcounties()
    if rivers:
        bm.drawrivers()
    return bm

####################
#   Plot modules   #
####################


def plot_polar_contour(values, azimuths, zeniths, nlevs=30,
                       vmin=-24., vmax=60., cmap='jet',
                       ax=None):
    """
    Plot a polar contour plot, with 0 degrees at the North.

    NOTE: Need to call .common.create_polar_fig_ax to
    create properly projected axis.

    Parameters
    ----------
    values: float
        A list (or other iterable - eg. a NumPy array) of
        the values to plot on the contour plot (the `z` values).
    azimuths: float
        A list of azimuths (in degrees).
    zeniths: float
        A list of zeniths (that is, radii).
    vmin : float
        Luminance minimum value, None for default value.
    vmax : float
        Luminance maximum value, None for default value.
    cmap : str
        Matplotlib colormap name.

    The shapes of these lists are important, and are designed for a particular
    use case (but should be more generally useful).
    The values list should be `len(azimuths) * len(zeniths)` long with
    data for the first azimuth for all the zeniths.
    Then the second azimuth for all the zeniths etc.

    This is designed to work nicely with data that is produced
    using a loop as follows:

    values = []
    for azimuth in azimuths:
      for zenith in zeniths:
        # Do something and get a result
        values.append(result)

    After that code the azimuths, zeniths and values lists will be ready
    to be passed into this function.

    Obtained from http://blog.rtwilson.com/
        producing-polar-contour-plots-with-matplotlib/
    """
    # parse parameters
    ax = _parse_ax(ax)
    # Make sure that the data is numpy array
    zeniths = np.array(zeniths)
    values = np.array(values)

    # Resize it to work
    values = values.reshape(len(azimuths), len(zeniths))

    # Create 2D variables to plot contour against and convert degrees to
    # radians
    r, theta = np.meshgrid(zeniths, np.radians(azimuths))

    # Default Matplotlib is math polar coords, want meteorological polar coords
    ax.set_theta_zero_location("N")  # This makes "North" set to 0 degrees
    ax.set_theta_direction(-1)  # This makes the angles increase clockwise

    p = ax.pcolormesh(theta, r, values, cmap=cmap, vmin=vmin, vmax=vmax)
    return p


def plot_fill_surface(xarr, surface, color=None, ymin=None, ax=None):
    """
    Add filled surface to plot (e.g. time-height).

    Parameters
    ----------
    xarr : array
        Array of x values.
    surface : array
        Array of surface height in meters. Same size as xarr.
    color : str
        Color string to use to fill below surface. If none
        defaults to dark gray.
    min : float
        Minimum number to use for fell. If None defaults to 0.
    ax : Matplotlib axis instance
        Axis to plot. None will use the current axis.
    """
    # parse parameters
    ax = _parse_ax(ax)

    if ymin is None:
        ymin = 0.
    if color is None:
        color = '0.85'

    p = ax.fill_between(xarr, ymin, surface, facecolor=color)
    return


def plot_xy(var1, var2, color=None, lw=None, ls=None, marker=None, msize=None,
            ax=None):
    """
    Returns an X-Y plot of variable1 vs. variable2.

    Parameters
    ----------
    var1 : array
        Variable to plot on x-axis.
    var2 : array
        Variable to plot on y-axis.
    color : str
        Line color.
    lw : float
        Linewidth to display
    ls : str
        Linestyle to use, can be abbreviation or name.
    marker : str
        Marker to display.
    msize : float
        Marker size.
    ax : Matplotlib axis instance
        Axis to plot. None will use the current axis.
    """
    # Parse parameters
    ax = _parse_ax(ax)

    # Set parameters if none
    if color is None:
        color = 'k'
    if lw is None:
        lw = 2

    ax.plot(var1, var2, color=color, ls=ls, lw=lw,
            marker=marker, markersize=msize, markeredgecolor=color)
    return


def image_2d(xvar, yvar, data_var,
             plot_log10_var=False, vmin=None, vmax=None, clevs=25,
             cmap=None, norm=None,
             x_major_ticks=None, x_minor_ticks=None,
             y_major_ticks=None, y_minor_ticks=None,
             x_min=None, x_max=None, y_min=None, y_max=None,
             title=None, titleFontSize=None,
             xlab=None, xlabFontSize=None, xpad=None,
             ylab=None, ylabFontSize=None, ypad=None,
             color_bar=True, cb_orient=None,
             cb_fontsize=None, cb_ticklabel_size=None,
             cb_pad=None, cb_tick_int=None, cb_label=None,
             ax=None, fig=None):
    """
    Returns a time series plot, with time on X-axis and variable on Y-axis.

    Parameters
    ----------
    xvar : float
        Array for x-axis (2D same as data_var).
    ax_var : float
        Array for y-axis (2D same as data_var).
    data_var : float
        Data variable.
    vmin : float
        Minimum contour value to display.
    vmax : float
        Maximum contour value to display.
    clevs : int
        Number of levels to use in colorbar tick calculation.
    cmap : str
        Matplotlib color map to use.
    norm : Matplotlib.colors.Normalize instance
        Matplotlib normaliztion instance used to scale luminance data.
    title : str
        Plot title.
    titleFontSize : int
        Font size to use for Title label.
    xlab : str
        X-axis label.
    ylab : str
        Y-axis label.
    xpad : int
        Padding for X-axis label.
    ypad : int
        Padding for Y-axis label.
    xlabFontSize : int
        Font size to use for X-axis label.
    ylabFontSize : int
        Font size to use for Y-axis label.
    x_major_ticks : float
        Values for x-axis major tickmark spacing.
    x_minor_ticks : float
        Values for x-axis minor tickmark spacing.
    y_major_ticks : float
        Values for y-axis major tickmark spacing.
    y_minor_ticks : float
        Values for y-axis minor tickmark spacing.
    x_min : float
        Minimum value for x-axis.
    x_max : float
        Maximum value for x-axis.
    x_min : float
        Minimum value for y-axis.
    x_max : float
        Maximum value for y-axis.
    color_bar : bool
        True to add colorbar, False does not.
    cb_pad : str
        Pad to move colorbar, in the form "5%",
        pos is to right for righthand location.
    cb_orient : str
        Colorbar orientation, either 'vertical' or 'horizontal'.
    cb_fontsize : int
        Font size of the colorbar label.
    cb_ticklabel_size : int
        Font size of colorbar tick labels.
    cb_tick_int : int
        Interval to use for colorbar tick labels, higher number "thins" labels.
    cb_label : str
        String to use as colorbar label.
    ax : Matplotlib axis instance
        Axis to plot. None will use the current axis.
    fig : Matplotlib figure instance
        Figure on which to add the plot. None will use the current figure.
    """
    # parse parameters
    ax = _parse_ax(ax)
    # If no cmap is specified, grab current
    if cmap is None:
        cmap = cm.get_cmap()

    # Set up the axes
    _set_axes(ax, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max,
              x_major_ticks=x_major_ticks, x_minor_ticks=x_minor_ticks,
              y_major_ticks=y_major_ticks, y_minor_ticks=y_minor_ticks,
              title=title, titleFontSize=titleFontSize,
              xlab=xlab, ylab=ylab, xpad=xpad, ypad=ypad,
              xlabFontSize=xlabFontSize, ylabFontSize=ylabFontSize)

    # Create the plot
    p = ax.pcolormesh(xvar, yvar, data_var,
                      vmin=vmin, vmax=vmax, norm=norm, cmap=cmap)

    # Add Colorbar
    if color_bar:
        cb = add_colorbar(ax, p, orientation=cb_orient, pad=cb_pad,
                          label=cb_label, fontsize=cb_fontsize,
                          ticklabel_size=cb_ticklabel_size,
                          clevs=clevs, tick_interval=cb_tick_int)
    return


def plot_date_ts(time, Var, color='k', marker='o', msize=1.5, lw=2,
                 date_format='%H:%M', tz=None, xdate=True,
                 date_minor_string='minute',
                 other_major_ticks=None, other_minor_ticks=None,
                 other_min=None, other_max=None,
                 title=None, titleFontSize=None,
                 xlab=None, xlabFontSize=None, xpad=None,
                 ylab=None, ylabFontSize=None, ypad=None,
                 ax=None):
    """
    Returns a time series plot, with time on X-axis and variable on Y-axis.

    Parameters
    ----------
    time : float
        Time array to plot on x-axis.
    Var : float
        Variable to plot as time series.
    color : str
        Color of marker.
    marker : str
        Marker to display.
    msize : float
        Marker size.
    date_format : str
        Format of the time string for x-axis labels.
    tz : str
        Time zone info to use when creating axis labels (see datetime).
    xdate : boolean
        True to use X-axis as date axis, false implies Y-axis is date axis.
    title : str
        Plot title.
    titleFontSize : int
        Font size to use for Title label.
    xlab : str
        X-axis label.
    ylab : str
        Y-axis label.
    xpad : int
        Padding for X-axis label.
    ypad : int
        Padding for Y-axis label.
    xlabFontSize : int
        Font size to use for X-axis label.
    ylabFontSize : int
        Font size to use for Y-axis label.
    date_minor_string : str
        Sting to set minor ticks of date axis,
        'second','minute','hour','day' supported.
    other_major_ticks : float
        Values for major tickmark spacing, non-date axis.
    other_minor_ticks : float
        Values for minor tickmark spacing, non-date axis.
    other_min : float
        Minimum value for non-date axis.
    other_max : float
        Maximum value for non-date axis.
    ax : Matplotlib axis instance
        Axis to plot. None will use the current axis.
    """
    # parse parameters
    ax = _parse_ax(ax)

    if xdate:
        ydate = False
    else:
        ydate = True

    # Set up the axes
    _set_ts_axes(ax, date_format=date_format, tz=tz, xdate=xdate,
                 date_minor_string=date_minor_string,
                 other_major_ticks=other_major_ticks,
                 other_minor_ticks=other_minor_ticks,
                 other_min=other_min, other_max=other_max,
                 title=title, titleFontSize=titleFontSize,
                 xlab=xlab, xlabFontSize=xlabFontSize, xpad=xpad,
                 ylab=ylab, ylabFontSize=ylabFontSize, ypad=ypad)
    # Create the plot
    ax.plot_date(time, Var, tz=tz, xdate=xdate, ydate=ydate,
                 mfc=color, mec=color, marker=marker,
                 markersize=msize, lw=lw)
    return


def image_2d_date(time, ax_var, data_var,
                  plot_log10_var=False,
                  vmin=None, vmax=None, clevs=25,
                  cmap=None, norm=None,
                  date_format='%H:%M', tz=None, xdate=True,
                  date_minor_string='minute',
                  other_major_ticks=None, other_minor_ticks=None,
                  other_min=None, other_max=None,
                  title=None, titleFontSize=None,
                  xlab=None, xlabFontSize=None, xpad=None,
                  ylab=None, ylabFontSize=None, ypad=None,
                  color_bar=True, cb_orient=None,
                  cb_fontsize=None, cb_ticklabel_size=None,
                  cb_pad=None, cb_tick_int=None,
                  cb_label=None,
                  ax=None, fig=None):
    """
    Returns a time series plot, with time on X-axis and variable on Y-axis.

    Parameters
    ----------
    time : float
        time array (2D same as data_var) to plot on x-axis.
    ax_var : float
        Variable (2D same as data_var) to use as other axis variable for plot.
    data_var : float
        Variable to plot as time series.
    vmin : float
        Minimum contour value to display.
    vmax : float
        Maximum contour value to display.
    clevs : int
        Number of levels to use in colorbar tick calculation.
    cmap : str
        Matplotlib color map to use.
    norm : Matplotlib.colors.Normalize instance
        Matplotlib normaliztion instance used to scale luminance data.
    date_format : str
        Format of the time string for x-axis labels.
    tz : str
        Time zone info to use when creating axis labels (see datetime).
    xdate : bool
        True to use X-axis as date axis, false implies Y-axis is date axis.
    title : str
        Plot title.
    titleFontSize : int
        Font size to use for Title label.
    xlab : str
        X-axis label.
    ylab : str
        Y-axis label.
    xpad : int
        Padding for X-axis label.
    ypad : int
        Padding for Y-axis label.
    xlabFontSize : int
        Font size to use for X-axis label.
    ylabFontSize : int
        Font size to use for Y-axis label.
    date_minor_string : str
        Sting to set minor ticks of date axis,
        'second','minute','hour','day' supported.
    other_major_ticks : float
        Values for major tickmark spacing, non-date axis.
    other_minor_ticks : float
        Values for minor tickmark spacing, non-date axis.
    other_min : float
        Minimum value for non-date axis.
    other_max : float
        Maximum value for non-date axis.
    color_bar : bool
        True to add colorbar, False does not.
    cb_pad : str
        Pad to move colorbar, in the form "5%",
        pos is to right for righthand location.
    cb_orient : str
        Colorbar orientation, either 'vertical' or 'horizontal'.
    cb_fontsize : int
        Font size of the colorbar label.
    cb_ticklabel_size : int
        Font size of colorbar tick labels.
    cb_tick_int : int
        Interval to use for colorbar tick labels, higher number "thins" labels.
    cb_label : str
        String to use as colorbar label.
    ax : Matplotlib axis instance
        Axis to plot. None will use the current axis.
    fig : Matplotlib figure instance
        Figure on which to add the plot. None will use the current figure.
    """
    # parse parameters
    ax = _parse_ax(ax)
    # If no cmap is specified, grab current
    if cmap is None:
        cmap = cm.get_cmap()

    # Set the axes variables depending on which is time axis
    if xdate:
        ydate = False
        xvar = time
        yvar = ax_var
    else:
        ydate = True
        xvar = ax_var
        yvar = time

    # Set up the axes
    _set_ts_axes(ax, date_format=date_format, tz=tz, xdate=xdate,
                 date_minor_string=date_minor_string,
                 other_major_ticks=other_major_ticks,
                 other_minor_ticks=other_minor_ticks,
                 other_min=other_min, other_max=other_max,
                 title=title, titleFontSize=titleFontSize,
                 xlab=xlab, xlabFontSize=xlabFontSize, xpad=xpad,
                 ylab=ylab, ylabFontSize=ylabFontSize, ypad=ypad)

    # Create the plot
    p = ax.pcolormesh(xvar, yvar, data_var,
                      vmin=vmin, vmax=vmax, norm=norm, cmap=cmap)

    # Add Colorbar
    if color_bar:
        cb = add_colorbar(ax, p, orientation=cb_orient, pad=cb_pad,
                          label=cb_label, fontsize=cb_fontsize,
                          ticklabel_size=cb_ticklabel_size,
                          clevs=clevs, tick_interval=cb_tick_int)
    return

#######################
#   General Methods   #
#######################


def find_nearest_indices(array, values):
    """
    Find the nearest value indices in an array to input value(s).

    Parameters
    ----------
    array : float array
        Input array to search.
    values  : float (array)
        Value(s) for which to search.
    """
    # Set the values to a 1D structure
    values = np.atleast_1d(values)
    # Find the nearest neighbor indices
    indices = np.abs(np.subtract.outer(array, values)).argmin(0)
    return indices if len(indices) > 1 else indices[0]


def _set_axes(ax, x_min=None, x_max=None,
              y_min=None, y_max=None,
              x_major_ticks=None, x_minor_ticks=None,
              y_major_ticks=None, y_minor_ticks=None,
              title=None, titleFontSize=None,
              xlab=None, xlabFontSize=None, xpad=None,
              ylab=None, ylabFontSize=None, ypad=None):
    """
    Returns a time series plot, with time on X-axis and variable on Y-axis.

    Parameters
    ----------
    ax : Matplotlib axis instance
        Axis to use.
    x_min : float
        Minimum value for X-axis.
    x_max : float
        Maximum value for X-axis.
    y_min : float
        Minimum value for Y-axis.
    y_max : float
        Maximum value for Y-axis.
    x_major_ticks : float
        Values for x-axis major tickmark spacing.
    x_minor_ticks : float
        Values for x-axis tickmark spacing.
    y_major_ticks : float
        Values for y-axis major tickmark spacing.
    y_minor_ticks : float
        Values for y-axis tickmark spacing.
    title : str
        Plot title.
    titleFontSize : int
        Font size to use for Title label.
    xlab : str
        X-axis label.
    ylab : str
        Y-axis label.
    xpad : int
        Padding for X-axis label.
    ypad : int
        Padding for Y-axis label.
    xlabFontSize : int
        Font size to use for X-axis label.
    ylabFontSize : int
        Font size to use for Y-axis label.
    """
    # Potentially set the X- and Y-axis limits
    if x_min is not None:
        ax.set_xlim(left=x_min)
    if x_max is not None:
        ax.set_xlim(right=x_max)

    if y_min is not None:
        ax.set_ylim(bottom=y_min)
    if y_max is not None:
        ax.set_ylim(top=y_max)

    # Set the major and minor y-axis ticks/tickmarks
    if x_major_ticks is not None:
            ax.xaxis.set_major_locator(mtic.MultipleLocator(x_major_ticks))
    if x_minor_ticks is not None:
            ax.xaxis.set_minor_locator(mtic.MultipleLocator(x_minor_ticks))
    if y_major_ticks is not None:
            ax.yaxis.set_major_locator(mtic.MultipleLocator(y_major_ticks))
    if y_minor_ticks is not None:
            ax.yaxis.set_minor_locator(mtic.MultipleLocator(y_minor_ticks))

    # Turn the tick marks outward
    ax.tick_params(which='both', direction='out')

    # Set the axis labels if provided
    if ylab is not None:
        if ylabFontSize is None:
            ylabFontSize = 16
        if ypad is None:
            ypad = 7
        ax.set_ylabel(ylab, labelpad=ypad, fontsize=ylabFontSize)
    if xlab is not None:
        if xlabFontSize is None:
            xlabFontSize = 16
        if xpad is None:
            xpad = 7
        ax.set_xlabel(xlab, labelpad=xpad, fontsize=xlabFontSize)

    # Set the title
    if title is not None:
        if titleFontSize is None:
            titleFontSize = 16
        ax.set_title(title, fontsize=titleFontSize)


def _set_ts_axes(ax, date_format='%H:%M', tz=None, xdate=True,
                 date_minor_string='minute',
                 other_major_ticks=None, other_minor_ticks=None,
                 other_min=None, other_max=None,
                 title=None, titleFontSize=None,
                 xlab=None, xlabFontSize=None, xpad=None,
                 ylab=None, ylabFontSize=None, ypad=None):
    """
    Returns a time series plot, with time on X-axis and variable on Y-axis.

    Parameters
    ----------
    ax : Matplotlib axis instance
        Axis to use.
    date_format : str
        Format of the time string for x-axis labels.
    tz : str
        Time zone info to use when creating axis labels (see datetime).
    xdate : bool
        True to use X-axis as date axis, false implies Y-axis is date axis.
    title : str
        Plot title.
    titleFontSize : int
        Font size to use for Title label.
    xlab : str
        X-axis label.
    ylab : str
        Y-axis label.
    xpad : int
        Padding for X-axis label.
    ypad : int
        Padding for Y-axis label.
    xlabFontSize : int
        Font size to use for X-axis label.
    ylabFontSize : int
        Font size to use for Y-axis label.
    date_minor_string : str
        Sting to set minor ticks of date axis,
        'second','minute','hour','day' supported.
    other_major_ticks : float
        Values for major tickmark spacing, non-date axis.
    other_minor_ticks : float
        Values for minor tickmark spacing, non-date axis.
    other_min : float
        Minimum value for non-date axis.
    other_max : float
        Maximum value for non-date axis.
    """
    # Set the date format
    date_Fmt = DateFormatter(date_format, tz=tz)

    # Check which axis to set for date axis and set tick parameters
    if xdate:
        # Set the x-axis date format and ticks
        ax.xaxis.set_major_formatter(date_Fmt)
        if date_minor_string == 'second':
            ax.xaxis.set_minor_locator(SecondLocator())
        elif date_minor_string == 'minute':
            ax.xaxis.set_minor_locator(MinuteLocator())
        elif date_minor_string == 'hour':
            ax.xaxis.set_minor_locator(HourLocator())
        elif date_minor_string == 'day':
            ax.xaxis.set_minor_locator(DayLocator())

        # Set the major and minor y-axis ticks/tickmarks
        try:
            other_major_ticks
            ax.yaxis.set_major_locator(mtic.MultipleLocator(other_major_ticks))
        except:
            pass
        try:
            other_minor_ticks
            ax.yaxis.set_minor_locator(mtic.MultipleLocator(other_minor_ticks))
        except:
            pass

        if other_min is not None:
            ax.set_ylim(bottom=other_min)
        if other_max is not None:
            ax.set_ylim(top=other_max)

    else:
        # Set the y-axis date format and ticks
        ax.yaxis.set_major_formatter(date_Fmt)
        if date_minor_string == 'second':
            ax.yaxis.set_minor_locator(SecondLocator())
        elif date_minor_string == 'minute':
            ax.yaxis.set_minor_locator(MinuteLocator())
        elif date_minor_string == 'hour':
            ax.yaxis.set_minor_locator(HourLocator())
        elif date_minor_string == 'day':
            ax.yaxis.set_minor_locator(DayLocator())

        # Set the major and minor x-axis ticks/tickmarks
        try:
            other_major_ticks
            ax.xaxis.set_major_locator(mtic.MultipleLocator(other_major_ticks))
        except:
            pass
        try:
            other_minor_ticks
            ax.xaxis.set_minor_locator(mtic.MultipleLocator(other_minor_ticks))
        except:
            pass

        if other_min is not None:
            ax.set_xlim(left=other_min)
        if other_max is not None:
            ax.set_xlim(right=other_max)

    # Turn the tick marks outward
    ax.tick_params(which='both', direction='out')

    # Set the axis labels if provided
    if ylab is not None:
        if ylabFontSize is None:
            ylabFontSize = 16
        if ypad is None:
            ypad = 7
        ax.set_ylabel(ylab, labelpad=ypad, fontsize=ylabFontSize)
    if xlab is not None:
        if xlabFontSize is None:
            xlabFontSize = 16
        if xpad is None:
            xpad = 7
        ax.set_xlabel(xlab, labelpad=xpad, fontsize=xlabFontSize)

    # Set the title
    if title is not None:
        if titleFontSize is None:
            titleFontSize = 16
        ax.set_title(title, fontsize=titleFontSize)
    return


def add_colorbar(ax, plot_instance, orientation=None, pad=None,
                 label=None, fontsize=None, ticklabel_size=None,
                 clevs=None, tick_interval=None):
    '''
    Add Colorbar to a plot instance.

    Parameters
    ----------
    ax : Matplotlib axis instance
        Axis to plot.
    plot : Matplotlib plot instance
        Plot instance used to establish colorbar.
    orientation : str
        Colorbar orientation, either 'vertical' or 'horizontal'.
    pad : str
        Pad to move colorbar, in the form "5%",
        positive is to right for righthand location.
    label : str
        String to use as colorbar label.
    fontsize : int
        Font size of the colorbar label.
    ticklabel_size : int
        Font size of colorbar tick labels.
    clevs : int
        Number of colorbar levels to use in tick calculation.
    tick_interval : int
        Interval to use for colorbar tick labels, higher number "thins" labels.
    '''
    if orientation is None:
        orientation = 'vertical'
    if pad is None:
        pad = .05
    if tick_interval is None:
        tick_interval = 2
    cb = plt.colorbar(plot_instance, orientation=orientation,
                      pad=pad, ax=ax)
    if label is not None:
        cb.set_label(label, fontsize=fontsize)
    # Set the tick label size
    cb.ax.tick_params(labelsize=ticklabel_size)
    # Set the number of ticks in the colorbar based upon number of contours
    if (clevs is not None) & (tick_interval is not None):
        tick_locator = mtic.MaxNLocator(nbins=int(clevs / tick_interval))
        cb.locator = tick_locator
        cb.update_ticks()
    return cb


def create_polar_fig_ax(nrows=1, ncols=1, figsize=(5, 5)):
    '''
    Returns the figure and axes instance of a polar plot.

    Parameters
    ----------
    nrows : int
        Number of rows.
    ncols : int
        Number of columns.
    figsize : tuple
        (xsize, ysize) in inches of figure.
    '''
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols,
                           subplot_kw=dict(projection='polar'),
                           fig_kw=dict(figsize=figsize))
    return fig, ax

##################
#  Save methods  #
##################


def save_figure(name='awot_plot', type="png", dpi=300):
    '''Save the current plot.

    Parameters
    ------------
    name : str
        Figure name.
    type : str
        Figure format, default to .png file type.
    dpi : int
        Resolution in dots per inch.
    '''
    plt.savefig(name+'.'+type, format=figType, dpi=dpi)

#################
#  Get methods  #
#################


def get_masked_data(data, mask_procedure, mask_tuple):
    '''Get a masked variable from dictionary'''
    if mask_procedure.lower() == 'less':
        vmin = mask_tuple[1]
        data = np.ma.masked_less(data, vmin)
    elif mask_procedure.lower() == 'less_equal':
        vmin = mask_tuple[1]
        data = np.ma.masked_less_equal(data, vmin)
    elif mask_procedure.lower() == 'greater':
        vmax = mask_tuple[1]
        data = np.ma.masked_greater(data, vmax)
    elif mask_procedure.lower() == 'greater_equal':
        vmax = mask_tuple[1]
        data = np.ma.masked_greater_equal(data, vmax)
    elif mask_procedure.lower() == 'equal':
        vmeq = mask_tuple[1]
        data = np.ma.masked_equal(data, vmeq)
    elif mask_procedure.lower() == 'inside':
        vmin = mask_tuple[1]
        vmax = mask_tuple[2]
        data = np.ma.masked_inside(data, vmin, vmax)
    elif mask_procedure.lower() == 'outside':
        vmin = mask_tuple[1]
        vmax = mask_tuple[2]
        data = np.ma.masked_outside(data, vmin, vmax)
    else:
        print("Check the mask_procedure operation string!")
    return data


def _get_start_datetime(time, start_time):
    '''Get a start time as datetime instance for subsetting.'''
    # Check to see if time is subsetted
    if start_time is None:
        dt_start = time['data'].min()
    else:
        startStr = [start_time[0:4], start_time[5:7], start_time[8:10],
                    start_time[11:13], start_time[14:16],
                    start_time[17:19], '0']
        startInt = [int(s) for s in startStr]
        try:
            dt_start = datetime(
                startInt[0], startInt[1], startInt[2], startInt[3],
                startInt[4], startInt[5], startInt[6])
        except:
            import warnings
            warnings.warn(common.DATE_STRING_FORMAT)
            return

    # Check to see if date time specified is beyond start
    if dt_start < time['data'].min():
        import warnings
        warnings.warn("WARNING: Specified START time occurs before the first "
                      "time instance. Using start of time array instead.")
        dt_start = time['data'].min()
    return dt_start


def _get_end_datetime(time, end_time):
    '''Get a start time as datetime instance for subsetting.'''
    # Check to see if the time is subsetted
    if end_time is None:
        dt_end = time['data'].max()
    else:
        endStr = [end_time[0:4], end_time[5:7], end_time[8:10],
                  end_time[11:13], end_time[14:16], end_time[17:19], '0']
        endInt = [int(s) for s in endStr]
        try:
            dt_end = datetime(endInt[0], endInt[1], endInt[2], endInt[3],
                              endInt[4], endInt[5], endInt[6])
        except:
            import warnings
            warnings.warn(common.DATE_STRING_FORMAT)
            return

    # Check to see if date time specified is beyond start
    if dt_end > time['data'].max():
        import warnings
        warnings.warn("WARNING: Specified END time occurs after the last "
                      "time instance. Using end of time array instead.")
        dt_end = time['data'].max()
    return dt_end


def _get_variable_dict(dict, field):
    ''' Return variable dictionary of specfied field. '''
    Var = dict[field]
    return Var


def _get_variable_dict_data(dict, field):
    ''' Return variable dictionary and data of specfied field. '''
    Var, data = dict[field], dict[field]['data'][:]
    return Var, data

########################
#   Parsing Methods    #
########################


def _parse_ax_fig(ax, fig):
    """ Parse and return ax and fig parameters. """
    if ax is None:
        ax = plt.gca()
    if fig is None:
        fig = plt.gcf()
    return ax, fig


def _parse_ax(ax):
    """ Parse and return ax parameters. """
    if ax is None:
        ax = plt.gca()
    return ax


def _parse_fig(fig):
    """ Parse and return fig parameters. """
    if fig is None:
        fig = plt.gcf()
    return fig

#####################
#   Check methods   #
#####################


def _check_basemap(instance, strong=True):
    """
    Check for a basemap instance.

    Parameters
    ----------
    instance: Class instance
        An instance to check for a basemap instance attached.
    strong: bool
        If True and the check finds no basemap instance,
        a ValueError is raised.
        If False and the check finds no basemap instance,
        a Warning is issued.
    """
    if instance.basemap is None:
        if strong:
            raise ValueError('Please supply basemap instance')
            return None
        else:
            print("WARNING: A basemap instance may be required for some plots")
            return False


def _check_field(instance, field):
    """
    Check to see if a field has a value, return if it does not
    Parameters
    ----------
    instance: dict
        A data dictionary to check.
    field : string
        The key name of the field to check
    """
    if instance[field] is None:
        print("This field has no value!!")
        return
