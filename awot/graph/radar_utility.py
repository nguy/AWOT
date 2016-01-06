"""
awot.graph.radar_utility
=========================

A group of scripts to create various radar utility plots.
"""

# Load the needed packages
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import date2num, num2date
from scipy.stats import mstats

from .common import (_parse_ax, _set_axes,
                     _get_start_datetime, _get_end_datetime,
#                     _get_variable_dict, _get_variable_dict_data,
                     add_colorbar)

from ..util.helper import add_dict_to_awot_fields

class RadarUtilityPlot(object):
    """Create vertical plot of radar data."""

    def __init__(self, radar, basemap=None,
                lon_name=None, lat_name=None, height_name=None,
                time_name=None, surface_name=None):
        """
        Initialize the class to create plots

        Parameters
        ----------
        radar : dict
            AWOT radar instance.
        basemap : basemap instance
        lon_name : str
            Key in radar instance for longitude variable.
            None uses AWOT default.
        lat_name : str
            Key in radar instance for latitude variable.
            None uses AWOT default.
        height_name : str
            Key in radar instance for height variable.
            None uses AWOT default.
        time_name : str
            Key in radar instance for time variable.
            None uses AWOT default.
        surface_name : str
            Key in radar instance for surface height variable.
            None uses AWOT default.
        """
        self.radar = radar
        self.fields = self.radar['fields']

        if lon_name is None:
            self.longitude = self.radar['longitude']
        else:
            self.longitude = self.radar[lon_name]
        if lat_name is None:
            self.latitude = self.radar['latitude']
        else:
            self.latitude = self.radar[lat_name]
        if height_name is None:
            self.height = self.radar['height']
        else:
            self.height = self.radar[height_name]

        # Attempt to pull in time if found
        if time_name is None:
            try:
                self.time = self.radar['time']
            except:
                self.time = None
        else:
            self.time = self.radar[time_name]

        # Attempt to pull in surface height if found
        if surface_name is None:
            try:
                self.surface = self.radar['surface']
            except:
                self.surface = None
        else:
            self.surface = self.radar[surface_name]

        # See if time or height arrays are 2D, if not create
        fieldshape = self.radar['fields'][self.radar['fields'].keys()[0]]['data'].shape
        self.timefield = self.time.copy()
        self.heightfield = self.height.copy()

        # This step is slow converting back-and-forth between datenum
        # instances, but Numpy does not support meshgrid with datenum dtype
        if len(self.timefield['data'].shape) == 1:
            timenum = date2num(self.time['data'][:], self.time['units'])
            self.timefield['data'] = num2date(np.resize(timenum, fieldshape),
                                              self.time['units'])

        if len(self.heightfield['data'].shape) == 1:
            self.heightfield['data'] = np.resize(self.height['data'][:], fieldshape)

##################
#  plot methods  #
##################

    def plot_bivariate_frequency(self, xfield, yfield,
                             xbinsminmax=None, nbinsx=50,
                             ybinsminmax=None, nbinsy=50,
                             start_time=None, end_time=None,
                             vmin=None, vmax=None, cmap=None,
                             mask_below=None,
                             plot_percent=False, plot_colorbar=True,
                             x_min=None, x_max=None,
                             y_min=None, y_max=None,
                             xlab=None, xlabFontSize=None, xpad=None,
                             ylab=None, ylabFontSize=None, ypad=None,
                             title=None, titleFontSize=None,
                             cb_fontsize=None, cb_ticklabel_size=None,
                             cb_orient=None, cb_pad=None,
                             cb_levs=None, cb_tick_int=None,
                             ax=None, fig=None):
        """
        Create a bivariate frequency distribution plot of two variables.

        Parameters
        ----------
        xfield : str
            Name of the field to use for the x array in calculation.
        yfield : str
            Name of the field to use for the y array in calculation.
        xbinsminmax : 2-tuple
            A tuple with the minimum and maximax values to
            use with xarr. None will use min/max of xarr.
        nbinsx : int
            The number of bins to use with xarr, default is 50.
        ybinsminmax : array
            A tuple with the minimum and maximax values to
            use with yarr. None will use min/max of yarr.
        nbinsy : int
            The number of bins to use with yarr, default is 50.
        start_time : str
            UTC time to use as start time for subsetting in datetime format.
            (e.g. 2014-08-20 12:30:00)
        end_time : str
            UTC time to use as an end time for subsetting in datetime format.
            (e.g. 2014-08-20 16:30:00)
        vmin : float
            Minimum value to display.
        vmax : float
            Maximum value to display.
        cmap : str
            Matplotlib colormap string.
        mask_below : float
            If provided, values less than mask_below will be masked.
        plot_percent : boolean
            True to display percentage. Default is to display fraction.
        plot_colorbar : boolean
            True to diaplay colorbar. False does not display colorbar.
        x_min : float
            Minimum value for X-axis.
        x_max : float
            Maximum value for X-axis.
        y_min : float
            Minimum value for Y-axis.
        y_max : float
            Maximum value for Y-axis.
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
        title : str
            Plot title.
        titleFontSize : int
            Font size to use for Title label.
        cb_fontsize : int
            Font size of the colorbar label.
        cb_ticklabel_size : int
            Font size of colorbar tick labels.
        cb_pad : str
            Pad to move colorbar, in the form "5%",
            pos is to right for righthand location.
        cb_orient : str
            Colorbar orientation, either 'vertical' or 'horizontal'.
        cb_levs : int
            Number of colorbar levels to use in tick calculation.
        cb_tick_int : int
            Interval to use for colorbar tick labels,
            higher number "thins" labels.
        ax : Matplotlib axis instance
            Axis to plot. None will use the current axis.
        fig : Matplotlib figure instance
            Figure on which to add the plot. None will use the current figure.
        """
        # parse parameters
        ax = _parse_ax(ax)
        # Snag the data from requested fields
        xarr, yarr = self._get_bivariate_data(
                      xfield, yfield, start_time, end_time)

        if xbinsminmax is None:
            xbinsminmax = (np.ma.min(xarr), np.ma.max(xarr))
        if ybinsminmax is None:
            ybinsminmax = (np.ma.min(yarr), np.ma.max(yarr))
        binsx = np.linspace(xbinsminmax[0], xbinsminmax[1],
                            nbinsx, endpoint=True)
        binsy = np.linspace(ybinsminmax[0], ybinsminmax[1],
                            nbinsy, endpoint=True)

        CFAD, xedges, yedges = np.histogram2d(
                 xarr.ravel(), yarr.ravel(), bins=(binsx, binsy), normed=True)
        X, Y = np.meshgrid(xedges, yedges)
        if mask_below is not None:
            CFAD = np.ma.masked_where(CFAD < mask_below, CFAD)

        cb_label = "Frequency"
        if plot_percent:
            CFAD = CFAD * 100.
            cb_label = cb_label + " (%)"

        # Set the axes
        _set_axes(x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max,
                  title=title, titleFontSize=titleFontSize,
                  xlab=xlab, ylab=ylab, xpad=xpad, ypad=ypad,
                  xlabFontSize=xlabFontSize, ylabFontSize=ylabFontSize,
                  ax=ax)
        # Plot the data
        p = ax.pcolormesh(X, Y, CFAD.T, vmin=vmin, vmax=vmax, cmap=cmap)

        if plot_colorbar:
            cb = add_colorbar(ax, p, orientation=cb_orient, pad=cb_pad,
                 label=cb_label, fontsize=cb_fontsize,
                 ticklabel_size=cb_ticklabel_size,
                 clevs=cb_levs, tick_interval=cb_tick_int)

        # Create a dictionary to return
        bivar = {'frequency' : CFAD.T,
                 'frequency_label' : cb_label,
                 'xaxis' : X,
                 'yaxis' : Y}
        return bivar

    def plot_cfad(self, field, height_axis=1,
                  xbinsminmax=None, nbinsx=50,
                  start_time=None, end_time=None,
                  vmin=None, vmax=None, cmap=None,
                  mask_below=None, plot_percent=False,
                  plot_colorbar=True,
                  x_min=None, x_max=None,
                  y_min=None, y_max=None,
                  xlab=None, xlabFontSize=None, xpad=None,
                  ylab=None, ylabFontSize=None, ypad=None,
                  title=None, titleFontSize=None,
                  cb_fontsize=None, cb_ticklabel_size=None,
                  cb_orient=None, cb_pad=None,
                  cb_levs=None, cb_tick_int=None,
                  ax=None, fig=None):
        """
        Create a frequency by altitude distribution plot of two variables.
        This is the traditional method of calculating a frequency distribution
        at each height of input array by iterating through the height array
        and input data array.

        Parameters
        ----------
        field : str
            Name of the field to use in CFAD calculation.
        height_axis : int
            The dimension to use as the height axis.
        xbinsminmax : 2-tuple
            A tuple with the minimum and maximax values to
            use with xarr. None will use min/max of xarr.
        nbinsx : int
            The number of bins to use with xarr, default is 50.
        start_time : str
            UTC time to use as start time for subsetting in datetime format.
            (e.g. 2014-08-20 12:30:00)
        end_time : str
            UTC time to use as an end time for subsetting in datetime format.
            (e.g. 2014-08-20 16:30:00)
        vmin : float
            Minimum value to display.
        vmax : float
            Maximum value to display.
        cmap : str
            Matplotlib colormap string.
        mask_below : float
            If provided, values less than mask_below will be masked.
        plot_percent : boolean
            True to display percentage. Default is to display fraction.
        plot_colorbar : boolean
            True to diaplay colorbar. False does not display colorbar.
        x_min : float
            Minimum value for X-axis.
        x_max : float
            Maximum value for X-axis.
        y_min : float
            Minimum value for Y-axis.
        y_max : float
            Maximum value for Y-axis.
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
        title : str
            Plot title.
        titleFontSize : int
            Font size to use for Title label.
        cb_fontsize : int
            Font size of the colorbar label.
        cb_ticklabel_size : int
            Font size of colorbar tick labels.
        cb_pad : str
            Pad to move colorbar, in the form "5%",
            pos is to right for righthand location.
        cb_orient : str
            Colorbar orientation, either 'vertical' or 'horizontal'.
        cb_levs : int
            Number of colorbar levels to use in tick calculation.
        cb_tick_int : int
            Interval to use for colorbar tick labels,
            higher number "thins" labels.
        ax : Matplotlib axis instance
            Axis to plot. None will use the current axis.
        fig : Matplotlib figure instance
            Figure on which to add the plot. None will use the current figure.
        """
        # parse parameters
        ax = _parse_ax(ax)
        # Snag the data from requested field
        xdict, tsub, xarr = self._get_fields_variable_dict_data_time_subset(
            field, start_time, end_time)

        if xbinsminmax is None:
            xbinsminmax = (np.ma.min(xarr), np.ma.max(xarr))
        binsx = np.linspace(xbinsminmax[0], xbinsminmax[1], nbinsx, endpoint=True)

        cb_label = "Frequency"
        percent = False

        if plot_percent:
            cb_label = cb_label + " (%)"
            percent = True

        # Create CFAD array to fill
        nh = len(self.height['data'][:])
        CFAD = np.empty((nh, len(binsx)-1))
        for nn in range(nh):
            if height_axis == 0:
                array = xarr[nn, ...]
            elif height_axis == 1:
                array = xarr[:, nn, ...]
            elif height_axis == 2:
                array = xarr[..., nn]
            CFAD[nn, :], bin_edges = np.histogram(
                   array, bins=binsx, density=percent)

        # Mask any invalid or negative numbers
        CFAD = np.ma.masked_invalid(CFAD)
        CFAD = np.ma.masked_less(CFAD, 0.)

        X, Y = np.meshgrid(bin_edges, self.height['data'][:])
        if mask_below is not None:
            CFAD = np.ma.masked_where(CFAD < mask_below, CFAD)

        # Set the axes
        _set_axes(x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max,
                  title=title, titleFontSize=titleFontSize,
                  xlab=xlab, ylab=ylab, xpad=xpad, ypad=ypad,
                  xlabFontSize=xlabFontSize, ylabFontSize=ylabFontSize,
                  ax=ax)
        # Plot the data
        p = ax.pcolormesh(X, Y, CFAD, vmin=vmin, vmax=vmax, cmap=cmap)

        if plot_colorbar:
            cb = add_colorbar(ax, p, orientation=cb_orient, pad=cb_pad,
                 label=cb_label, fontsize=cb_fontsize,
                 ticklabel_size=cb_ticklabel_size,
                 clevs=cb_levs, tick_interval=cb_tick_int)

        # Create a dictionary to return
        cfad_dict = {'frequency' : CFAD,
                     'frequency_label' : cb_label,
                     'xaxis' : X,
                     'yaxis' : Y}
        return cfad_dict

#     def plot_quantiles(self, field, quantiles=None,
#                   height_axis=1,
#                   xbinsminmax=None, nbinsx=50,
#                   start_time=None, end_time=None,
#                   mask_below=None, plot_percent=False,
#                   plot_colorbar=True,
#                   contour_levels_color='k',
#                   x_min=None, x_max=None,
#                   y_min=None, y_max=None,
#                   xlab=None, xlabFontSize=None, xpad=None,
#                   ylab=None, ylabFontSize=None, ypad=None,
#                   title=None, titleFontSize=None,
#                   cb_fontsize=None, cb_ticklabel_size=None,
#                   setup_axes=True,
#                   ax=None, fig=None):
#         """
#         Create a frequency by altitude distribution plot of two variables.
#         This is the traditional method of calculating a frequency distribution
#         at each height of input array by iterating through the height array
#         and input data array.
#
#         Parameters
#         ----------
#         field : str
#             Name of the field to use in CFAD calculation.
#         height_axis : int
#             The dimension to perform quantile calculation over (non-height).
#         xbinsminmax : 2-tuple
#             A tuple with the minimum and maximax values to
#             use with xarr. None will use min/max of xarr.
#         nbinsx : int
#             The number of bins to use with xarr, default is 50.
#         start_time : str
#             UTC time to use as start time for subsetting in datetime format.
#             (e.g. 2014-08-20 12:30:00)
#         end_time : str
#             UTC time to use as an end time for subsetting in datetime format.
#             (e.g. 2014-08-20 16:30:00)
#         mask_below : float
#             If provided, values less than mask_below will be masked.
#         plot_percent : boolean
#             True to display percentage. Default is to display fraction.
#         plot_colorbar : boolean
#             True to diaplay colorbar. False does not display colorbar.
#         x_min : float
#             Minimum value for X-axis.
#         x_max : float
#             Maximum value for X-axis.
#         y_min : float
#             Minimum value for Y-axis.
#         y_max : float
#             Maximum value for Y-axis.
#         title : str
#             Plot title.
#         titleFontSize : int
#             Font size to use for Title label.
#         xlab : str
#             X-axis label.
#         ylab : str
#             Y-axis label.
#         xpad : int
#             Padding for X-axis label.
#         ypad : int
#             Padding for Y-axis label.
#         xlabFontSize : int
#             Font size to use for X-axis label.
#         ylabFontSize : int
#             Font size to use for Y-axis label.
#         cb_fontsize : int
#             Font size of the colorbar label.
#         setup_axes : bool
#             If True, the axes will be set up using either keyword information
#             or attempted to determine automatically.
#             If False, existing axes will be used, either passed by keyword ax
#             or current working axis.
#         ax : Matplotlib axis instance
#             Axis to plot. None will use the current axis.
#         fig : Matplotlib figure instance
#             Figure on which to add the plot. None will use the current figure.
#         """
#         # parse parameters
#         ax = _parse_ax(ax)
#         # Snag the data from requested field
#         xdict, tsub, xarr = self._get_fields_variable_dict_data_time_subset(
#             field, start_time, end_time)
#
#         if quantiles is None:
#             quantiles = [5, 10, 25, 50, 75, 90]
#
#         if xbinsminmax is None:
#             xbinsminmax = (np.ma.min(xarr), np.ma.max(xarr))
#         binsx = np.linspace(xbinsminmax[0], xbinsminmax[1], nbinsx, endpoint=True)
#
#         cb_label = "Frequency"
#         percent = False
#
#         if plot_percent:
#             cb_label = cb_label + " (%)"
#
#         # Create CFAD array to fill
#         nh = len(self.height['data'][:])
# #         CFAD = np.empty((nh, len(binsx)-1))
#         qArr = np.empty((len(quantiles), nh))
#         xarr = np.rollaxis(xarr, height_axis)
#         for nn in range(nh):
#             array = xarr[nn, ...]
#             qArr[:, nn] = mstats.mquantiles(
#                     array, prob=quantiles)
#             print(array.min(), array.max())
#             print(qArr.min(), qArr.max())
#
# #        qArr = mstats.mquantiles(xarr, prob=quantiles, axis=axis)
#
#         if mask_below is not None:
#             CFAD = np.ma.masked_where(CFAD < mask_below, CFAD)
#
#         # Set the axes
#         if setup_axes:
#             _set_axes(x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max,
#                       title=title, titleFontSize=titleFontSize,
#                       xlab=xlab, ylab=ylab, xpad=xpad, ypad=ypad,
#                       xlabFontSize=xlabFontSize, ylabFontSize=ylabFontSize,
#                       ax=ax)
#         # Plot the data
#         for num in range(len(quantiles)):
#             p = ax.plot(qArr[num, :], self.heightfield['data'][0, :])#, color='k')
#             ytextloc = self.heightfield['data'][0, :].max() * 1.05
#             ax.text(qArr[num, -1], ytextloc, str(quantiles[num]))
# #            print(qArr[num, :].min(), qArr[num, :].max())
# #            print(qArr[num, -1], ytextloc)
#         return

###################
#   Get methods   #
###################

    def _get_fields_variable_dict_data_time_subset(self, field, start_time, end_time):
        '''
        Get the variable from the fields dictionary.
        Subset the time when in time series format.
        '''
        Var, data = self.fields[field], self.fields[field]['data'][:]

        # Check to see if time is subsetted
        dt_start = _get_start_datetime(self.time, start_time)
        dt_end = _get_end_datetime(self.time, end_time)
        tsub = self.time['data'][(self.time['data'] >= dt_start) &
                                 (self.time['data'] <= dt_end)]
        datasub = data[(self.time['data'] >= dt_start) &
                       (self.time['data'] <= dt_end)]
        return Var, tsub, datasub

    def _get_variable_dict_data_time_subset(self, field, start_time, end_time):
        '''
        Get the variable from the fields dictionary.
        Subset the time when in time series format.
        '''
        Var, data = self.radar[field], self.radar[field]['data'][:]

        # Check to see if time is subsetted
        dt_start = _get_start_datetime(self.time, start_time)
        dt_end = _get_end_datetime(self.time, end_time)
        tsub = self.time['data'][(self.time['data'] >= dt_start) &
                                 (self.time['data'] <= dt_end)]
        datasub = data[(self.time['data'] >= dt_start) &
                       (self.time['data'] <= dt_end)]
        return Var, tsub, datasub

    def _get_bivariate_data(self, field1, field2,
                                  start_time, end_time):
        '''
        Return the the data for field1 and field 2.
        Subset the time when in time series format.
        '''
        if field1 != 'height':
            try:
                data1 = self.fields[field1]['data']
            except:
                data1 = self.radar[field1]['data']
        else:
            data1 = self.heightfield['data']
        if field2 != 'height':
            try:
                data2 = self.fields[field2]['data']
            except:
                data2 = self.radar[field2]['data']
        else:
            data2 = self.heightfield['data']

        # Check to see if time is subsetted
        dt_start = _get_start_datetime(self.time, start_time)
        dt_end = _get_end_datetime(self.time, end_time)

        data1sub = data1[(self.timefield['data'] >= dt_start) &
                       (self.timefield['data'] <= dt_end)]

        data2sub = data2[(self.timefield['data'] >= dt_start) &
                       (self.timefield['data'] <= dt_end)]
        return data1sub, data2sub

    def _get_2d_height_time_subset(self, start_time, end_time):
        '''Get subsetted data if requested.'''
        # Check to see if time is subsetted
        dt_start = _get_start_datetime(self.time, start_time)
        dt_end = _get_end_datetime(self.time, end_time)

        hsub = self.height['data'][(self.time['data'] >= dt_start) &
                                    (self.time['data'] <= dt_end), :]
        hsub = np.ma.masked_invalid(hsub)
        return hsub

    def _get_lat_index(self, value):
        '''Calculate the exact index position within latitude array.'''
        # Find the spacing
        dp = self.latitude['data'][1] - self.latitude['data'][0]

        # Calculate the relative position
        pos = (value - self.latitude['data'][0]) / dp
        return pos

    def _get_lon_index(self, value):
        '''Calculate the exact index position within latitude array'''
        # Find the spacing
        dp = self.longitude['data'][1] - self.longitude['data'][0]

        # Calculate the relative position
        pos = (value - self.longitude['data'][0]) / dp
        return pos
