"""
awot.graph.microphysical_vertical
=========================

A group of scripts to create vertical radar plots.
"""
# FUNCTIONS::
# polar_sweep - Plot polar coordinate data on polar coordinate axis
# polar_sweep_grid - Plotting transformed data to Cartesian output
# sweep_to_Cart - Polar coordinates transformed to Cartesian
# sweep_aircraft_relative - Polar coord data transformed
#                           to aircraft-relative frame
# sweep_track_relative - Polar coord data transformed to track-relative frame
# sweep_earth_relative - Polar coord data transformed to earth-relative frame

# Load the needed packages
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter, date2num
from matplotlib import ticker
from matplotlib.colors import from_levels_and_colors
import scipy.ndimage as scim

from . import common
from .coord_transform import (radar_coords_to_cart_track_relative,
                              radar_coords_to_cart_earth_relative,
                              radar_coords_to_cart_aircraft_relative)


class MicrophysicalVerticalPlot(object):
    """
    Create a vertical plot of microphysical data.

    This class was created for microphysical retrievals from radar systems,
    for example the LATMOS French Falcon.
    """

    def __init__(self, microphysdata, basemap=None,
                 lon_name=None, lat_name=None, height_name=None,
                 time_name=None):
        """
        Intitialize the class to create plots

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
        """
        # Now initialize the RadarHorizontalPlot Class
        self.microphys_data = microphysdata
        self.basemap = basemap
#        common._check_basemap(self)
        self.fields = microphys_data['fields']

        if lon_name is None:
            self.longitude = self.microphys_data['longitude']
        else:
            self.longitude = self.microphys_data[lon_name]
        if lat_name is None:
            self.latitude = self.microphys_data['latitude']
        else:
            self.latitude = self.microphys_data[lat_name]
        if height_name is None:
            self.height = self.microphys_data['height']
        else:
            self.height = self.microphys_data[height_name]
        if time_name is None:
            self.time = self.microphys_data['time']
        else:
            self.time = self.microphys_data[time_name]

#############################
#   Vertical plot methods   #
#############################

    def time_height_image(self, field,
                          mask_procedure=None, mask_tuple=None,
                          plot_log10_var=False,
                          cminmax=(0., 60.), clevs=25,
                          vmin=None, vmax=None,
                          cmap='gist_ncar', discrete_cmap_levels=None,
                          date_format='%H:%M', tz=None, xdate=True,
                          date_minor_string='minute',
                          height_MajTicks=None, height_MinTicks=None,
                          height_min=None, height_max=None,
                          start_time=None, end_time=None,
                          title=None,
                          xlab=' ', xlabFontSize=16, xpad=7,
                          ylab=' ', ylabFontSize=16, ypad=7,
                          color_bar=True, cb_orient='vertical',
                          cb_pad=.05, cb_tick_int=2,
                          cb_label=None,
                          cb_fontsize=None, cb_ticklabel_size=None,
                          ax=None, fig=None):
        """
        Wrapper function to produce a contoured time series plot
        of variable indicated.

        Parameters
        ----------
        field : array
            Variable to plot as time series.
        mask_procedure : str
            String indicating how to apply mask via numpy, possibilities are:
            'less', 'less_equal', 'greater', 'greater_equal',
            'equal', 'inside', 'outside'.
        mask_tuple : (str, float[, float])
            Tuple containing the field name and value(s) below which to mask
            field prior to plotting.
        plot_log10_var : boolean
                True plots the log base 10 of Data field.
        cminmax : tuple
            (min, max) values for controur levels.
        clevs : int
            Number of contour levels.
        vmin : float
            Minimum contour value to display.
        vmax : float
            Maximum contour value to display.
        cmap : str
            Matplotlib color map to use.
        discrete_cmap_levels : array
            A list of levels to be used for display. If chosen discrete
            color will be used in the colorbar instead of a linear luminance
            mapping.
        date_format : str
            Format of the time string for x-axis labels.
        tz : str
            Time zone info to use when creating axis labels (see datetime).
        xdate : bool
            True to use X-axis as date axis, False implies Y-axis is date axis.
        date_minor_string : str
            Sting to set minor ticks of date axis,
            'second','minute','hour','day' supported.
        height_MajTicks : float
            Values for major tickmark spacing for height axis.
        height_MinTicks : float
            Values for minor tickmark spacing for height axis.
        height_min : float
            Minimum value for height axis.
        other_max : float
            Maximum value for height axis.
        start_time : str
            UTC time to use as start time for subsetting in datetime format.
            (e.g. 2014-08-20 12:30:00)
        end_time : str
            UTC time to use as an end time for subsetting in datetime format.
            (e.g. 2014-08-20 16:30:00)
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
        color_bar : boolean
            True adds colorbar, False does not.
        cb_pad : str
            Pad to move colorbar, in the form "5%",
            pos is to right for righthand location.
        cb_orient : str
            Colorbar orientation, either 'vertical' or 'horizontal'.
        cb_tick_int : int
            Interval to use for colorbar tick labels,
            higher number "thins" labels.
        cb_label : str
            Label for colorbar (e.g. units 'dBZ').
        cb_fontsize : int
            Font size of the colorbar label.
        cb_ticklabel_size : int
            Font size of colorbar tick labels.
        ax : Matplotlib axes instance
            Optional axes instance to plot the graph.
        fig : Matplotlib figure instance
            Figure which to add the plot.
            None will use the current figure.
        """
        # parse parameters
        ax, fig = common._parse_ax_fig(ax, fig)

        # Return masked or unmasked variable
        # Subsetted if desired
        Var, Data = self._get_variable_dict_data_time_subset(
            field, start_time, end_time)
        if mask_procedure is not None:
            Data = common.get_masked_data(Data, mask_procedure, mask_tuple)

        if plot_log10_var:
            Data = np.log10(Data)
            if cb_label is not None:
                cb_label = r'log$_{10}$[' + cb_label + ']'

        # Print out min and max values to screen
        print("Minimum value of %s = %g" % (field, np.ma.min(Data)))
        print("Maximum value of %s = %g" % (field, np.ma.max(Data)))

        # Get vmin/vmax if not given
        if vmin is None:
            vmin = np.ma.min(Data)
        if vmax is None:
            vmax = np.ma.max(Data)

        # Create contour level array
        clevels = np.linspace(cminmax[0], cminmax[1], clevs)

        tSub2D, Ht2D = np.meshgrid(date2num(tsub), self.height['data'][:])

        # Get the colormap and calculate data spaced by number of levels
        norm = None
        if discrete_cmap_levels is not None:
            cm = plt.get_cmap(cmap)
            try:
                levpos = np.rint(np.squeeze(
                    [np.linspace(0, 255,
                                 len(discrete_cmap_levels))])).astype(int)
                # Convert levels to colormap values
                cmap, norm = from_levels_and_colors(
                    discrete_cmap_levels, cm(levpos), extend='max')
            except:
                print("Keyword error: 'discrete_cmap_levels' must "
                      "be a list of float or integer")

        # Plot the time series
        ts = common.image_2d_date(tSub2D, Ht2D, Data,
                                  vmin=vmin, vmax=vmax, clevs=clevs,
                                  date_format=date_format, tz=tz, xdate=xdate,
                                  date_minor_string=date_minor_string,
                                  other_major_ticks=height_MajTicks,
                                  other_minor_ticks=height_MinTicks,
                                  other_min=height_min, other_max=height_max,
                                  title=title,
                                  xlab=xlab, xlabFontSize=xlabFontSize,
                                  xpad=xpad,
                                  ylab=ylab, ylabFontSize=ylabFontSize,
                                  ypad=ypad,
                                  color_bar=color_bar, cb_orient=cb_orient,
                                  cb_pad=cb_pad, cb_tick_int=cb_tick_int,
                                  cb_label=cb_label,
                                  cb_fontsize=cb_fontsize,
                                  cb_ticklabel_size=cb_ticklabel_size,
                                  ax=ax, fig=fig)
        return

    def track_height_image(self, field, track_key=None,
                           mask_procedure=None, mask_tuple=None,
                           plot_log10_var=False,
                           cminmax=(0., 60.), clevs=25,
                           vmin=None, vmax=None,
                           cmap='gist_ncar', discrete_cmap_levels=None,
                           x_min=None, x_max=None,
                           height_MajTicks=None, height_MinTicks=None,
                           height_min=None, height_max=None,
                           title=None,
                           xlab=' ', xlabFontSize=16, xpad=7,
                           ylab=' ', ylabFontSize=16, ypad=7,
                           color_bar=True, cb_orient='vertical',
                           cb_pad=.05, cb_tick_int=2,
                           cb_label=None,
                           cb_fontsize=None, cb_ticklabel_size=None,
                           ax=None, fig=None):
        """
        Wrapper function to produce a contoured plot along a track
        of variable indicated.

        Parameters
        ----------
        field : array
            Variable to plot as time series.
        track_key : str
            Key name of track distance variable.
        mask_procedure : str
            String indicating how to apply mask via numpy, possibilities are:
            'less', 'less_equal', 'greater', 'greater_equal',
            'equal', 'inside', 'outside'.
        mask_tuple : (str, float[, float])
            Tuple containing the field name and value(s) below which to mask
            field prior to plotting.
        plot_log10_var : boolean
                True plots the log base 10 of Data field.
        cminmax : tuple
            (min, max) values for controur levels.
        clevs : int
            Number of contour levels.
        vmin : float
            Minimum contour value to display.
        vmax : float
            Maximum contour value to display.
        cmap : str
            Matplotlib color map to use.
        discrete_cmap_levels : array
            A list of levels to be used for display. If chosen discrete
            color will be used in the colorbar instead of a linear luminance
            mapping.
        height_MajTicks : float
            Values for major tickmark spacing for height axis.
        height_MinTicks : float
            Values for minor tickmark spacing for height axis.
        height_min : float
            Minimum value for height axis.
        height_max : float
            Maximum value for height axis.
        x_min : float
            Minimum value for track distance axis.
        x_max : float
            Maximum value for track distance axis.
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
        color_bar : boolean
            True adds colorbar, False does not.
        cb_pad : str
            Pad to move colorbar, in the form "5%",
            pos is to right for righthand location.
        cb_orient : str
            Colorbar orientation, either 'vertical' or 'horizontal'.
        cb_tick_int : int
            Interval to use for colorbar tick labels,
            higher number "thins" labels.
        cb_label : str
            Label for colorbar (e.g. units 'dBZ').
        cb_fontsize : int
            Font size of the colorbar label.
        cb_ticklabel_size : int
            Font size of colorbar tick labels.
        ax : Matplotlib axes instance
            Optional axes instance to plot the graph.
        fig : Matplotlib figure instance
            Figure which to add the plot.
            None will use the current figure.
        """
        # parse parameters
        ax, fig = common._parse_ax_fig(ax, fig)

        # Return masked or unmasked variable
        # Subsetted if desired
        Var, tsub, Data = self._get_variable_dict_data(field)
        if mask_procedure is not None:
            Data = common.get_masked_data(Data, mask_procedure, mask_tuple)

        if track_key is None:
            try:
                track = self.flight_data['track_distance_air']
            except:
                ValueError('Did not find suitable track distance variable')
        else:
            track = self.flight_data[track_key]
        trackd = track['data'][:]

        if plot_log10_var:
            Data = np.log10(Data)
            if cb_label is not None:
                cb_label = r'log$_{10}$[' + cb_label + ']'

        # Print out min and max values to screen
        print("Minimum value of %s = %g" % (field, np.ma.min(Data)))
        print("Maximum value of %s = %g" % (field, np.ma.max(Data)))

        # Get vmin/vmax if not given
        if vmin is None:
            vmin = np.ma.min(Data)
        if vmax is None:
            vmax = np.ma.max(Data)

        # Create contour level array
        clevels = np.linspace(cminmax[0], cminmax[1], clevs)

        track2D, Ht2D = np.meshgrid(trackd, self.height['data'][:])

        # Get the colormap and calculate data spaced by number of levels
        norm = None
        if discrete_cmap_levels is not None:
            cm = plt.get_cmap(cmap)
            try:
                levpos = np.rint(np.squeeze(
                    [np.linspace(0, 255,
                                 len(discrete_cmap_levels))])).astype(int)
                # Convert levels to colormap values
                cmap, norm = from_levels_and_colors(
                    discrete_cmap_levels, cm(levpos), extend='max')
            except:
                print("Keyword error: 'discrete_cmap_levels' must "
                      "be a list of float or integer")

        # Plot the time series
        ts = common.image_2d(track2D, Ht2D, Data,
                             vmin=vmin, vmax=vmax, clevs=clevs,
                             date_format=date_format, tz=tz, xdate=xdate,
                             date_minor_string=date_minor_string,
                             other_major_ticks=height_MajTicks,
                             other_minor_ticks=height_MinTicks,
                             other_min=height_min, other_max=height_max,
                             title=title,
                             xlab=xlab, xlabFontSize=xlabFontSize,
                             xpad=xpad,
                             ylab=ylab, ylabFontSize=ylabFontSize,
                             ypad=ypad,
                             color_bar=color_bar, cb_orient=cb_orient,
                             cb_pad=cb_pad, cb_tick_int=cb_tick_int,
                             cb_label=cb_label,
                             cb_fontsize=cb_fontsize,
                             cb_ticklabel_size=cb_ticklabel_size,
                             ax=ax, fig=fig)
        return

#################
#  Get methods  #
#################

    def _get_variable_dict_data(self, field):
        '''
        Get the variable from the fields dictionary.
        Subset the time when in time series format.
        '''
        Var, data = self.fields[field], self.fields[field]['data'][:]
        return Var, data

    def _get_variable_dict_data_time_subset(self, field, start_time, end_time):
        '''
        Get the variable from the fields dictionary.
        Subset the time when in time series format.
        '''
        Var, data = self.fields[field], self.fields[field]['data'][:]

        # Check to see if time is subsetted
        dt_start = common._get_start_datetime(self.time, start_time)
        dt_end = common._get_end_datetime(self.time, end_time)

        # Create temporary 2D arrays for subsetting
        tsub = self.time['data'][(self.time['data'] >= dt_start) &
                                 (self.time['data'] <= dt_end)]
        datasub = data[(self.time['data'] >= dt_start) &
                       (self.time['data'] <= dt_end)]
        datasub = np.ma.masked_invalid(datasub)
        return Var, tsub, datasub

    def _get_variable_subset(self, data, start_time, end_time):
        '''
        Get the variable from the fields dictionary.
        Subset the time when in time series format.
        '''
        # Check to see if time is subsetted
        dt_start = common._get_start_datetime(self.time, start_time)
        dt_end = common._get_end_datetime(self.time, end_time)

        # Create temporary 2D arrays for subsetting
        tsub = self.time['data'][(self.time['data'] >= dt_start) &
                                 (self.time['data'] <= dt_end)]
        datasub = data[(self.time['data'] >= dt_start) &
                       (self.time['data'] <= dt_end)]
        datasub = np.ma.masked_invalid(datasub)
        return datasub

    def _get_2d_height_time_subset(self, start_time, end_time):
        '''Get subsetted data if requested.'''
        # Check to see if time is subsetted
        dt_start = common._get_start_datetime(self.time, start_time)
        dt_end = common._get_end_datetime(self.time, end_time)

        # Create temporary 2D arrays for subsetting
        hsub = self.height['data'][[(self.time['data'] >= dt_start) &
                                   (self.time['data'] <= dt_end)], :]
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
        '''Calculate the exact index position within latitude array.'''
        # Find the spacing
        dp = self.longitude['data'][1] - self.longitude['data'][0]

        # Calculate the relative position
        pos = (value - self.longitude['data'][0]) / dp
        return pos
