"""
awot.graph.common
==================

Common graphing routines.

Created by Nick Guy.
"""
# HISTORY::
# 20 Aug 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
# 
#-------------------------------------------------------------------
# Load the needed packages

import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from  matplotlib.dates import DateFormatter, SecondLocator, MinuteLocator, HourLocator, DayLocator
from matplotlib import ticker as mtic

#######################
# Map creation module #
#######################

def create_basemap_instance(corners=None, proj=None, resolution='l', area_thresh=1000.,
                   meridians=True, parallels=True, dLon=2., dLat=2.,
                   coastlines=True, countries=True, states=False, counties=False, 
                   rivers=False, etopo=False, ax=None):
    """
    Create a basemap instance for plotting.
    
    PARAMETERS::
    ----------
    corners : float tuple
        Array of map corners [lowerLeftLon, lowerLeftLat, upperRightLon, upperRightLat].
        Default will find min max of arrays
    proj : string
        Map projection to use
    resolution : string
        Map resolution, 'l' or low res is default
    area_thresh : float
        Minimum area threshold for basemap plot instance
    meridians : boolean
        Flag to turn on meridian (lonigitude) line plotting
    parallels : boolean
        Flag to turn on parallels (latitude) line plotting
    dLon : float
        Spacing for meridians
    dLat : float
        Spacing for parallels
    coastlines : boolean
        Flag to turn on basemap coastline plotting
    countries : boolean
        Flag to turn on basemap country plotting
    states : boolean
        Flag to turn on basemap state plotting
    counties : boolean
        Flag to turn on basemap county plotting
    rivers : boolean
        Flag to turn on basemap river plotting
    etopo : boolean
        Flag to turn on basemap etopo (topography) plotting
    ax : Axes instance
        Optional axes instance to add basemap to 
    """
                
    # parse parameters
    ax = _parse_ax(ax)
        
    # Check inputs
    if corners is None:
        corners = [-180., -90., 180., 90.]
    
    if proj is None:
        proj = 'cea'
        
    # Create a basemap instance
    bm = Basemap(projection=proj, resolution=resolution, area_thresh = area_thresh,
            llcrnrlon=corners[0], urcrnrlon=corners[2],
            llcrnrlat=corners[1], urcrnrlat=corners[3],
            ax=ax)
                   
    # Check the customizations for the basemap
    if meridians:
        bm.drawmeridians(np.arange(corners[0], corners[2], dLon), labels=[1,1,0,1])           
    if parallels:
        bm.drawparallels(np.arange(corners[1], corners[3], dLat), labels=[1,1,0,1])             
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

#########################
# Plot creation modules #
#########################

def plot_polar_contour(values, azimuths, zeniths, nlevs=30,
                       vmin=-24., vmax=60., cmap='jet',
                       ax=None):
    """Plot a polar contour plot, with 0 degrees at the North.
    
    NOTE: Need to call .common.create_polar_fig_ax to create properly projected axis
 
    Parameters::
    ----------
    values: float
        A list (or other iterable - eg. a NumPy array) of the values to plot on the
        contour plot (the `z` values)
    azimuths: float
        A list of azimuths (in degrees)
    zeniths: float
        A list of zeniths (that is, radii)
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
 
    This is designed to work nicely with data that is produced using a loop as follows:
 
    values = []
    for azimuth in azimuths:
      for zenith in zeniths:
        # Do something and get a result
        values.append(result)
 
    After that code the azimuths, zeniths and values lists will be ready 
    to be passed into this function.
    
    Obtained from http://blog.rtwilson.com/producing-polar-contour-plots-with-matplotlib/
 
    """
    # parse parameters
    ax = _parse_ax(ax)
    
    # Make sure that the data is numpy array
    zeniths = np.array(zeniths) 
    values = np.array(values) 
    
    # Resize it to work
    values = values.reshape(len(azimuths), len(zeniths)) 

    # Create 2D variables to plot contour against and convert degrees to radians
    r, theta = np.meshgrid(zeniths, np.radians(azimuths))
    
    # Default Matplotlib is math polar coords, want meteorological polar coords
    ax.set_theta_zero_location("N") # This makes "North" set to 0 degrees
    ax.set_theta_direction(-1) # This makes the angles increase clockwise
    
    p = ax.pcolormesh(theta,r,values,cmap=cmap,vmin=vmin,vmax=vmax)
    
    return p
    
##########

def plot_date_ts(Time, Var, colF='ko', msize=1.5, lw=2,
                dForm='%H:%M',tz=None, xdate=True, 
                date_MinTicker='minute',
                other_MajTicks=None, other_MinTicks=None,
                title=None,
                xlab=' ', xlabFontSize=16, xpad=7,
                ylab=' ', ylabFontSize=16, ypad=7,
                ax=None):
    """Returns a time series plot, with time on X-axis and variable on Y-axis.
    Parameters::
    ----------
    Time : float
        Time array to plot on x-axis
    Var : float
        Variable to plot as time series
    colF : str
        Color and marker shortcut (see python documentation)
    msize : float
        Marker size
    dForm : str
        Format of the time string for x-axis labels
    tz : str
        Time zone info to use when creating axis labels (see datetime)
    xdate : boolean
        True to use X-axis as date axis, false implies Y-axis is date axis
    title : str
        Plot title
    xlab : str
        X-axis label
    ylab : str
        Y-axis label
    xpad : int
        Padding for X-axis label
    ypad : int
        Padding for Y-axis label
    date_MinTicker : str
        Sting to set minor ticks of date axis,
        'second','minute','hour','day' supported
    other_MajTicks : float
        Values for major tickmark spacing, non-date axis
    other_MinTicks : float
        Values for minor tickmark spacing, non-date axis
    ax : Axes instance
        Optional axes instance to plot the graph
    """         
    # parse parameters
    ax = _parse_ax(ax)
    
    # Set the date format
    date_Fmt = DateFormatter(dForm,tz=tz)
    
    # Check which axis to set for date axis and set tick parameters
    if xdate:
        ydate = False
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
    
    else:
        ydate = True
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

    # Create the plot
    ax.plot_date(Time, Var, fmt=colF, tz=tz, xdate=xdate, ydate=ydate,
                 markersize=msize, lw=lw)
    

    ax.tick_params(which='both',direction='out') # Turn the tick marks outward
    
    # Set the Y label
    ax.set_ylabel(ylab,labelpad=ypad,fontsize=ylabFontSize)
    if xlab is None:
        pass
    else:
        ax.set_xlabel(xlab,labelpad=xpad,fontsize=xlabFontSize)
#    except:
#        pass
        
    # Set the title
    if title is None:
        pass
    else:
        ax.set_title(title)
    
    return
    
#####################
# General Methods #
#####################

def find_nearest_indices(array, values):
    """Find the nearest value indices in an array to input value(s)
    Parameters::
    ------------
    array : float array
        Input array to search
    values  : float (array)
        Value(s) to search for
    """
    # Set the values to a 1D structure
    values = np.atleast_1d(values)
    # Find the nearest neighbor indices
    indices = np.abs(np.subtract.outer(array, values)).argmin(0)
    return indices if len(indices) > 1 else indices[0]
    
#######

def create_polar_fig_ax(nrows=1, ncols=1, figsize=(5, 5)):
    ''' Returns the figure and axes instance of a polar plot
    Parameters::
    ----------
    nrows : int
        Number of rows
    ncols : int
        Number of columns
    figsize : tuple
        (xsize, ysize) in inches of figure 
    '''
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, 
                           subplot_kw=dict(projection='polar'), 
                           fig_kw=dict(figsize=figsize))
    
    return fig, ax
    
#######################
# Get masking methods #
#######################
       
def get_masked_data(data, mask_procedure, mask_tuple):
    '''Get a masked variable from dictionary'''
    if mask_procedure == 'less':
        vmin = mask_tuple[1]
        data = np.ma.masked_less(data, vmin)
    elif mask_procedure == 'less_equal':
        vmin = mask_tuple[1]
        data = np.ma.masked_less_equal(data, vmin)
    if mask_procedure == 'greater':
        vmin = mask_tuple[1]
        data = np.ma.masked_greater(data, vmax)
    elif mask_procedure == 'greater_equal':
        vmin = mask_tuple[1]
        data = np.ma.masked_greater_equal(data, vmax)
    elif mask_procedure == 'equal':
        veq = mask_tuple[1]
        data = np.ma.masked_equal(data, vmeq)
    elif mask_procedure == 'inside':
        vmin = mask_tuple[1]
        vmax = mask_tuple[2]
        data = np.ma.masked_inside(data, vmin, vmax)
    elif mask_procedure == 'oustide':
        vmin = mask_tuple[1]
        vmax = mask_tuple[2]
        data = np.ma.masked_outside(data, vmin, vmax)
    else:
        print "Check the mask_procedure operation string!"
            
    return data
    
#####################
# Parseing Methods #
#####################   

def _parse_ax_fig(ax, fig):
        """ Parse and return ax and fig parameters. """
        if ax is None:
            ax = plt.gca()
        if fig is None:
            fig = plt.gcf()
        return ax, fig        

######

def _parse_ax(ax):
        """ Parse and return ax parameters. """
        if ax is None:
            ax = plt.gca()
        return ax        

######

def _parse_fig(fig):
        """ Parse and return fig parameters. """
        if fig is None:
            fig = plt.gcf()
        return fig