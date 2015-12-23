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

###################
#   Map modules   #
###################


def create_basemap(corners=None, proj=None, resolution='l',
                   area_thresh=None, lon_0=None, lat_0=None,
                   meridians=True, parallels=True, dLon=2., dLat=2.,
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
    meridians : bool
        Flag to turn on meridian (lonigitude) line plotting.
    parallels : bool
        Flag to turn on parallels (latitude) line plotting.
    dLon : float
        Spacing for meridians.
    dLat : float
        Spacing for parallels.
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
                         2], dLon), labels=[1, 0, 0, 1])
    if parallels:
        bm.drawparallels(np.arange(corners[1], corners[
                         3], dLat), labels=[1, 0, 0, 1])
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

def plot_date_ts(Time, Var, color='k', marker='o', msize=1.5, lw=2,
                 dForm='%H:%M', tz=None, xdate=True,
                 date_MinTicker='minute',
                 other_MajTicks=None, other_MinTicks=None,
                 other_min=None, other_max=None,
                 title=None,
                 xlab=None, xlabFontSize=16, xpad=7,
                 ylab=None, ylabFontSize=16, ypad=7,
                 ax=None):
    """
    Returns a time series plot, with time on X-axis and variable on Y-axis.

    Parameters
    ----------
    Time : float
        Time array to plot on x-axis.
    Var : float
        Variable to plot as time series.
    color : str
        Color of marker.
    marker : str
        Marker to display.
    msize : float
        Marker size.

    dForm : str
        Format of the time string for x-axis labels.
    tz : str
        Time zone info to use when creating axis labels (see datetime).
    xdate : boolean
        True to use X-axis as date axis, false implies Y-axis is date axis.
    title : str
        Plot title.
    xlab : str
        X-axis label.
    ylab : str
        Y-axis label.
    xpad : int
        Padding for X-axis label.
    ypad : int
        Padding for Y-axis label.
    date_MinTicker : str
        Sting to set minor ticks of date axis,
        'second','minute','hour','day' supported.
    other_MajTicks : float
        Values for major tickmark spacing, non-date axis.
    other_MinTicks : float
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
    _set_ts_axes(dForm=dForm, tz=tz, xdate=xdate,
                 date_MinTicker=date_MinTicker,
                 other_MajTicks=other_MajTicks, other_MinTicks=other_MinTicks,
                 other_min=other_min, other_max=other_max,
                 title=title,
                 xlab=xlab, xlabFontSize=xlabFontSize, xpad=xpad,
                 ylab=ylab, ylabFontSize=ylabFontSize, ypad=ypad,
                 ax=ax)
    # Create the plot
    ax.plot_date(Time, Var, tz=tz, xdate=xdate, ydate=ydate,
                 mfc=color, mec=color, marker=marker,
                 markersize=msize, lw=lw)
    return


def image_2d_date(Time, AxVar, PlotVar,
                    ptype='pcolormesh', plot_log10_var=False,
                    vmin=None, vmax=None, clevs=25,
                    cmap=None,
                    color_bar=True, cb_orient='vertical',
                    cb_pad=.05, cb_tick_int=2,
                    cb_label=None,
                    dForm='%H:%M', tz=None, xdate=True,
                    date_MinTicker='minute',
                    other_MajTicks=None, other_MinTicks=None,
                    other_min=None, other_max=None,
                    title=None,
                    xlab=None, xlabFontSize=16, xpad=7,
                    ylab=None, ylabFontSize=16, ypad=7,
                    ax=None, fig=None):
    """
    Returns a time series plot, with time on X-axis and variable on Y-axis.

    Parameters
    ----------
    Time : float
        Time array (2D same as PlotVar) to plot on x-axis.
    AxVar : float
        Variable (2D same as PlotVar) to use as other axis variable for plot.
    PlotVar : float
        Variable to plot as time series.
    ptype : str
        Type of plot to make, takes 'plot', 'contour', or 'pcolormsh'.
    vmin : float
        Minimum contour value to display.
    vmax : float
        Maximum contour value to display.
    cmap : str
        Matplotlib color map to use.
    color_bar : bool
        True to add colorbar, False does not.
    cb_pad : str
        Pad to move colorbar, in the form "5%",
        pos is to right for righthand location.
    cb_loc : str
        Location of colorbar, default is 'right', also available:
        'bottom', 'top', 'left'.
    cb_tick_int : int
        Interval to use for colorbar tick labels, higher number "thins" labels.
    cb_label : str
        String to use as colorbar label.
    dForm : str
        Format of the time string for x-axis labels.
    tz : str
        Time zone info to use when creating axis labels (see datetime).
    xdate : bool
        True to use X-axis as date axis, false implies Y-axis is date axis.
    title : str
        Plot title.
    xlab : str
        X-axis label.
    ylab : str
        Y-axis label.
    xpad : int
        Padding for X-axis label.
    ypad : int
        Padding for Y-axis label.
    date_MinTicker : str
        Sting to set minor ticks of date axis,
        'second','minute','hour','day' supported.
    other_MajTicks : float
        Values for major tickmark spacing, non-date axis.
    other_MinTicks : float
        Values for minor tickmark spacing, non-date axis.
    other_min : float
        Minimum value for non-date axis.
    other_max : float
        Maximum value for non-date axis.
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
        XVar = Time
        YVar = AxVar
    else:
        ydate = True
        XVar = AxVar
        YVar = Time

    # Set up the axes
    _set_ts_axes(dForm=dForm, tz=tz, xdate=xdate,
                 date_MinTicker=date_MinTicker,
                 other_MajTicks=other_MajTicks, other_MinTicks=other_MinTicks,
                 other_min=other_min, other_max=other_max,
                 title=title,
                 xlab=xlab, xlabFontSize=xlabFontSize, xpad=xpad,
                 ylab=ylab, ylabFontSize=ylabFontSize, ypad=ypad,
                 ax=ax)

    # Create the plot
    if ptype == 'pcolormesh':
        p = ax.pcolormesh(XVar, YVar, PlotVar,
                          vmin=vmin, vmax=vmax, cmap=cmap)
    elif ptype == 'contour':
        p = ax.contourf(XVar, YVar, PlotVar,
                        vmin=vmin, vmax=vmax, cmap=cmap)

    # Add Colorbar
    if color_bar:
        cb = fig.colorbar(p, orientation=cb_orient,
                          pad=cb_pad, ax=ax)  # ,ticks=clevels)
        if cb_label is not None:
            cb.set_label(cb_label)
        # Set the number of ticks in the colorbar based upon number of contours
        tick_locator = mtic.MaxNLocator(nbins=int(clevs / cb_tick_int))
        cb.locator = tick_locator
        cb.update_ticks()
    return


def _set_ts_axes(dForm='%H:%M', tz=None, xdate=True,
                 date_MinTicker='minute',
                 other_MajTicks=None, other_MinTicks=None,
                 other_min=None, other_max=None,
                 title=None,
                 xlab=None, xlabFontSize=16, xpad=7,
                 ylab=None, ylabFontSize=16, ypad=7,
                 ax=None):
    """
    Returns a time series plot, with time on X-axis and variable on Y-axis.

    Parameters
    ----------
    dForm : str
        Format of the time string for x-axis labels.
    tz : str
        Time zone info to use when creating axis labels (see datetime).
    xdate : bool
        True to use X-axis as date axis, false implies Y-axis is date axis.
    title : str
        Plot title.
    xlab : str
        X-axis label.
    ylab : str
        Y-axis label.
    xpad : int
        Padding for X-axis label.
    ypad : int
        Padding for Y-axis label.
    date_MinTicker : str
        Sting to set minor ticks of date axis,
        'second','minute','hour','day' supported.
    other_MajTicks : float
        Values for major tickmark spacing, non-date axis.
    other_MinTicks : float
        Values for minor tickmark spacing, non-date axis.
    other_min : float
        Minimum value for non-date axis.
    other_max : float
        Maximum value for non-date axis.
    ax : Matplotlib axis instance
        Axis to plot. None will use the current axis.
    """
    # Set the date format
    date_Fmt = DateFormatter(dForm, tz=tz)

    # Check which axis to set for date axis and set tick parameters
    if xdate:
        # Set the x-axis date format and ticks
        ax.xaxis.set_major_formatter(date_Fmt)
        if date_MinTicker == 'second':
            ax.xaxis.set_minor_locator(SecondLocator())
        elif date_MinTicker == 'minute':
            ax.xaxis.set_minor_locator(MinuteLocator())
        elif date_MinTicker == 'hour':
            ax.xaxis.set_minor_locator(HourLocator())
        elif date_MinTicker == 'day':
            ax.xaxis.set_minor_locator(DayLocator())

        # Set the major and minor y-axis ticks/tickmarks
        try:
            other_MajTicks
            ax.yaxis.set_major_locator(mtic.MultipleLocator(other_MajTicks))
        except:
            pass
        try:
            other_MinTicks
            ax.yaxis.set_minor_locator(mtic.MultipleLocator(other_MinTicks))
        except:
            pass

        if other_min is not None:
            ax.set_ylim(bottom=other_min)
        if other_max is not None:
            ax.set_ylim(top=other_max)

    else:
        # Set the y-axis date format and ticks
        ax.yaxis.set_major_formatter(date_Fmt)
        if date_MinTicker == 'second':
            ax.yaxis.set_minor_locator(SecondLocator())
        elif date_MinTicker == 'minute':
            ax.yaxis.set_minor_locator(MinuteLocator())
        elif date_MinTicker == 'hour':
            ax.yaxis.set_minor_locator(HourLocator())
        elif date_MinTicker == 'day':
            ax.yaxis.set_minor_locator(DayLocator())

        # Set the major and minor x-axis ticks/tickmarks
        try:
            other_MajTicks
            ax.xaxis.set_major_locator(mtic.MultipleLocator(other_MajTicks))
        except:
            pass
        try:
            other_MinTicks
            ax.xaxis.set_minor_locator(mtic.MultipleLocator(other_MinTicks))
        except:
            pass

        if other_min is not None:
            ax.set_xlim(left=other_min)
        if other_max is not None:
            ax.set_xlim(right=other_max)

    # Turn the tick marks outward
    ax.tick_params(which='both', direction='out')

    # Set the Y label
    if ylab is not None:
        ax.set_ylabel(ylab, labelpad=ypad, fontsize=ylabFontSize)
    if xlab is not None:
        ax.set_xlabel(xlab, labelpad=xpad, fontsize=xlabFontSize)

    # Set the title
    if title is not None:
        ax.set_title(title)
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
            print(
                "Format of date string should be e.g. '2014-08-20 12:30:00'")
            return
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
            print(
                "Check the format of date string (e.g. '2014-08-20 12:30:00')")
            return
    return dt_end

def _get_variable_dict(dict, field):
    '''Get the variable from the fields dictionary.'''
    Var = dict[field]
    return Var

def _get_variable_dict_data(dict, field):
    '''Get the variable from the fields dictionary.'''
    Var, data = dict[field], dict[field]['data'][:]
    return Var, data

def _get_earth_radius():
    return 6371.

########################
#   Parsing Methods    #
########################


def _parse_ax_fig(ax, fig):
    """Parse and return ax and fig parameters."""
    if ax is None:
        ax = plt.gca()
    if fig is None:
        fig = plt.gcf()
    return ax, fig


def _parse_ax(ax):
    """Parse and return ax parameters."""
    if ax is None:
        ax = plt.gca()
    return ax


def _parse_fig(fig):
    """Parse and return fig parameters."""
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
