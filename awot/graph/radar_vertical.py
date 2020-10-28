"""
awot.graph.radar_vertical
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


class RadarVerticalPlot(object):
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
        self.basemap = basemap
#        common._check_basemap(self)
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

#############################
#  Vertical plot methods    #
#############################

    def plot_cross_section(self, field, start_pt, end_pt, xs_length=500,
                           mask_procedure=None, mask_tuple=None,
                           plot_km=False, title=" ", title_size=20,
                           cminmax=(0., 60.), clevs=25, vmin=15., vmax=60.,
                           cmap='gist_ncar', discrete_cmap_levels=None,
                           color_bar=True, clabel='dBZ', cb_pad=None,
                           cb_orient=None, cb_fontsize=None,
                           cb_ticklabel_size=None, cb_tick_int=None,
                           ax=None, fig=None):
        '''
        Plot a cross-section between two points.

        Parameters
        ----------
        field : str
            3-D variable (e.g. Reflectivity [dBZ]) to use in plot.
        start_pt, end_pt : tuple
            (lat, lon) Tuple of start, end points for cross-section.
        xs_length : int
            Number of to use for the cross section.
        mask_procedure : str
            String indicating how to apply mask via numpy, possibilities are:
            'less', 'less_equal', 'greater', 'greater_equal',
            'equal', 'inside', 'outside'.
        mask_tuple : (str, float[, float])
            Tuple containing the field name and value(s) below which to mask
            field prior to plotting, for example to mask all data where.
        plot_km : boolean
            True to convert meters to kilometers for cappi_height. False
            retains meters information.
        title : str
            Plot title.
        title_size : int
            Font size of title to display.
        cminmax : tuple
            (min,max) values for controur levels.
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
        color_bar : bool
            True to add colorbar, False does not.
        clabel : str
            Label for colorbar (e.g. units 'dBZ').
        cb_fontsize : int
            Font size of the colorbar label.
        cb_ticklabel_size : int
            Font size of colorbar tick labels.
        cb_pad : str
            Pad to move colorbar, in the form "5%", pos is to right for
            righthand location.
        cb_orient : str
            Colorbar orientation, either 'vertical' or 'horizontal'.
        cb_tick_int : int
            Interval to use for colorbar tick labels,
            higher number "thins" labels.
        ax : Matplotlib axis instance
            Axis to plot. None will use the current axis.
        fig : Matplotlib figure instance
            Figure on which to add the plot. None will use the current figure.
        '''
        # parse parameters
        ax, fig = common._parse_ax_fig(ax, fig)

        # Return masked or unmasked variable
        Var, Data = common._get_variable_dict_data(self.fields, field)
        if mask_procedure is not None:
            Data = common.get_masked_data(Data, mask_procedure, mask_tuple)

        # Create contour level array
        clevels = np.linspace(cminmax[0], cminmax[1], clevs)

        # Create lon and lat arrays for display
        xslon = np.linspace(start_pt[0], end_pt[0], xs_length)
        xslat = np.linspace(start_pt[1], end_pt[1], xs_length)

        # Create an array to hold the interpolated cross-section
        xs_data = np.empty([xs_length, len(self.height['data'][:])])

        # Create arrays for cross-section lon-lat points
        startloclon = self._get_lon_index(start_pt[0])
        startloclat = self._get_lat_index(start_pt[1])
        endloclon = self._get_lon_index(end_pt[0])
        endloclat = self._get_lat_index(end_pt[1])

        xsY = np.linspace(startloclat, endloclat, xs_length)
        xsX = np.linspace(startloclon, endloclon, xs_length)

        # Loop through each level to create cross-section and stack them
        for nlev in range(len(self.height['data'][:])):
            # Extract the values along the line, using cubic interpolation
            xs_data[:, nlev] = scim.map_coordinates(Data[nlev, :, :],
                                                    np.vstack((xsY, xsX)),
                                                    prefilter=False)
            # , mode='nearest')

        # Calculate the distance array along the cross-section
        Xdist = np.absolute(
            (np.pi * common.EARTH_RADIUS / 180.) * (xslon - xslon[0]))
        Ydist = np.absolute(
            (np.pi * common.EARTH_RADIUS / 180.) * (xslat - xslat[0]))
        xsDist = np.sqrt(Xdist**2 + Ydist**2)

        # Define the angle of the cross-secton
        Dlon = (start_pt[1] - end_pt[1])
        Dlat = (start_pt[0] - end_pt[0])
        Ang = np.arctan2(Dlat, Dlon)
        if Ang < 0:
            AngNref = 2 * np.pi + Ang
        else:
            AngNref = Ang

        # Convert Height, distance arrays to 2D
        if plot_km:
            Ht2D, Dist2D = np.meshgrid(self.height['data'][:]/1000., xsDist)
            ylabel = 'Altitude (km)'
        else:
            Ht2D, Dist2D = np.meshgrid(self.height['data'][:], xsDist)
            ylabel = 'Altitude (m)'

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
                print("--> Keyword error: 'discrete_cmap_levels' must "
                      "be a list of float or integer")

        p = ax.pcolormesh(Dist2D, Ht2D,
                          np.ma.masked_less_equal(xs_data, -800.),
                          vmin=vmin, vmax=vmax, norm=norm, cmap=cmap)

        ax.set_xlabel('Distance along track (km)')
        ax.set_ylabel(ylabel)

        # Add title
        ax.set_title(title, fontsize=title_size)

        # Add Colorbar
        if color_bar:
            cbStr = Var['long_name'] + ' ' + self.height['units']
            cb = common.add_colorbar(ax, p, orientation=cb_orient, pad=cb_pad,
                                     label=cbStr, fontsize=cb_fontsize,
                                     ticklabel_size=cb_ticklabel_size,
                                     clevs=clevs, tick_interval=cb_tick_int)

        # Add title
        ax.set_title(title, fontsize=title_size)

    def time_height_image(self, field, plot_km=False,
                          mask_procedure=None, mask_tuple=None,
                          plot_log10_var=False,
                          cminmax=(0., 60.), clevs=25, vmin=15., vmax=60.,
                          cmap='gist_ncar', discrete_cmap_levels=None,
                          date_format='%H:%M', tz=None, xdate=True,
                          date_minor_string='minute',
                          height_MajTicks=None, height_MinTicks=None,
                          height_min=None, height_max=None,
                          fill_surface=False, fill_min=None, fill_color=None,
                          start_time=None, end_time=None, color_bar=True,
                          cb_orient='vertical',
                          cb_pad=.05, cb_tick_int=2, cb_label=None,
                          cb_fontsize=None, cb_ticklabel_size=None,
                          title=None, xlab=' ', xlabFontSize=16, xpad=7,
                          ylab=' ', ylabFontSize=16, ypad=7,
                          ax=None, fig=None):
        """
        Wrapper function to produce a time series vs. height plot
        of variable indicated.

        Parameters
        ----------
        field : float
            Variable to plot as time series.
        plot_km : boolean
            True to convert meters to kilometers for cappi_height. False
            retains meters information.
        mask_procedure : str
            String indicating how to apply mask via numpy, possibilities are:
            'less', 'less_equal', 'greater', 'greater_equal',
            'equal', 'inside', 'outside'.
        mask_tuple : (str, float[, float])
            Tuple containing the field name and value(s) below which to mask
            field prior to plotting, for example to mask all data where.
        plot_log10_var : bool
                True plots the log base 10 of Data field.
        cminmax : tuple
            (min,max) values for controur levels.
        clevs : int
            Number of contour levels.
        vmin : float
            Minimum value to display.
        vmax : float
            Maximum value to display.
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
            True to use X-axis as date axis, false implies Y-axis is date axis.
        date_minor_string : str
            Sting to set minor ticks of date axis,
            'second','minute','hour','day' supported.
        height_MajTicks : float
            Values for major tickmark spacing on height axis.
        height_MinTicks : float
            Values for minor tickmark spacing on height axis.
        height_min : float
            Minimum value for height axis.
        height_max : float
            Maximum value for height axis.
        fill_surface : boolean
            True to fill in surface, False to leave alone.
        fill_min : float
            Minimum surface elvation to shade. Only applied
            if fill_surface is True.
        fill_color : float
            Color to use if fill_surface is True.
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
        color_bar : bool
            True to add colorbar, False does not.
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
        ax : Matplotlib axis instance
            Axis to plot. None will use the current axis.
        fig : Matplotlib figure instance
            Figure on which to add the plot. None will use the current figure.
        """
        # parse parameters
        ax, fig = common._parse_ax_fig(ax, fig)

        # Return masked or unmasked variable
        # Subsetted if desired
        Var, tsub, Data = self._get_variable_dict_data_time_subset(
            field, start_time, end_time)
        if mask_procedure is not None:
            Data = common.get_masked_data(Data, mask_procedure, mask_tuple)
        if len(self.height['data'].shape) == 2:
            Height = self._get_2d_height_time_subset(
                          start_time, end_time)

        if plot_log10_var:
            Data = np.log10(Data)
            if cb_label is not None:
                cb_label = r'log$_{10}$[' + cb_label + ']'

        # Create contour level array
        clevels = np.linspace(cminmax[0], cminmax[1], clevs)

        if len(self.height['data'].shape) == 2:
            tSub2D, junk = np.meshgrid(date2num(tsub),
                                       self.height['data'][0, :])
            Ht2D = Height.T
            del junk
        else:
            tSub2D, Ht2D = np.meshgrid(date2num(tsub), self.height['data'][:])

        if plot_km:
            Ht2D = Ht2D / 1000.

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
                print("--> Keyword error: 'discrete_cmap_levels' must "
                      "be a list of float or integer")

        # Plot the time series
        ts = common.image_2d_date(tSub2D, Ht2D, Data.T,
                                  vmin=vmin, vmax=vmax, clevs=clevs,
                                  cmap=cmap, norm=norm,
                                  date_format=date_format, tz=tz, xdate=xdate,
                                  date_minor_string=date_minor_string,
                                  other_major_ticks=height_MajTicks,
                                  other_minor_ticks=height_MinTicks,
                                  other_min=height_min,
                                  other_max=height_max,
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
        if fill_surface:
            if self.surface is not None:
                sfc = self._get_variable_subset(self.surface['data'][:],
                                                start_time, end_time)
                ft = common.plot_fill_surface(tsub, sfc, ymin=fill_min,
                                              color=fill_color, ax=ax)
            else:
                print("--> No surface height information, cannot fill...")
        return

    def track_height_image(self, field, track_key=None,
                           plot_height_km=False, plot_track_km=False,
                           mask_procedure=None, mask_tuple=None,
                           plot_log10_var=False,
                           cminmax=(0., 60.), clevs=25, vmin=None, vmax=None,
                           cmap='gist_ncar', discrete_cmap_levels=None,
                           track_min=None, track_max=None,
                           height_min=None, height_max=None,
                           start_time=None, end_time=None,
                           zero_track_start=False,
                           track_MajTicks=None, track_MinTicks=None,
                           height_MajTicks=None, height_MinTicks=None,
                           fill_surface=False, fill_min=None, fill_color=None,
                           color_bar=True, cb_orient='vertical',
                           cb_pad=.05, cb_tick_int=2,
                           cb_label=None,
                           cb_fontsize=None, cb_ticklabel_size=None,
                           title=None, xlab=' ', xlabFontSize=16, xpad=7,
                           ylab=' ', ylabFontSize=16, ypad=7,
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
        height_min : float
            Minimum value for height axis.
        height_max : float
            Maximum value for height axis.
        start_time : str
            UTC time to use as start time for subsetting in datetime format.
            (e.g. 2014-08-20 12:30:00)
        end_time : str
            UTC time to use as an end time for subsetting in datetime format.
            (e.g. 2014-08-20 16:30:00)
        zero_track_start : bool
            If True the first element in the track will be forced to zero.
            Of interest if using start_time keyword.
            Defaults to False.
        track_MajTicks : float
            Values for major tickmark spacing on track axis.
        track_MinTicks : float
            Values for minor tickmark spacing on track axis.
        height_MajTicks : float
            Values for major tickmark spacing on height axis.
        height_MinTicks : float
            Values for minor tickmark spacing on height axis.
        track_min : float
            Minimum value for track distance axis.
        track_max : float
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
        Var, tsub, Data = self._get_variable_dict_data_time_subset(
            field, start_time, end_time)
        if mask_procedure is not None:
            Data = common.get_masked_data(Data, mask_procedure, mask_tuple)

        if track_key is None:
            try:
                track = self.radar['track_distance_air']
            except:
                ValueError('Did not find suitable track distance variable')
        else:
            track = self.radar[track_key]
        trackd = self._get_variable_subset(track['data'][:], start_time,
                                           end_time)

        if plot_track_km:
            if track['units'] != 'km':
                trackd = trackd / 1000.

        if zero_track_start:
            trackd = trackd - trackd[0]

        if plot_log10_var:
            Data = np.log10(Data)
            if cb_label is not None:
                cb_label = r'log$_{10}$[' + cb_label + ']'

        # Print out min and max values to screen
        print(f"--> Minimum/Maximum value of {field} = {np.ma.min(Data)}/{np.ma.max(Data)}")

        # Get vmin/vmax if not given
        if vmin is None:
            vmin = np.ma.min(Data)
        if vmax is None:
            vmax = np.ma.max(Data)

        # Create contour level array
        clevels = np.linspace(cminmax[0], cminmax[1], clevs)

        if len(self.height['data'].shape) == 2:
            track2D, _ = np.meshgrid(trackd, self.height['data'][0, :])
            Ht2D = self.height['data'][:].T
        else:
            track2D, Ht2D = np.meshgrid(trackd, self.height['data'][:])

        if plot_height_km:
            Ht2D = Ht2D / 1000.

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
                print("--> Keyword error: 'discrete_cmap_levels' must "
                      "be a list of float or integer")

        # Plot the time series
        ts = common.image_2d(track2D, Ht2D, Data.T,
                             vmin=vmin, vmax=vmax, clevs=clevs,
                             cmap=cmap, norm=norm,
                             x_min=track_min, x_max=track_max,
                             y_min=height_min, y_max=height_max,
                             x_minor_ticks=track_MinTicks,
                             x_major_ticks=track_MajTicks,
                             y_minor_ticks=height_MinTicks,
                             y_major_ticks=height_MajTicks,
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
        if fill_surface:
            if self.surface is not None:
                sfc = self._get_variable_subset(self.surface['data'][:],
                                                start_time, end_time)
                ft = common.plot_fill_surface(trackd, sfc, ymin=fill_min,
                                              color=fill_color, ax=ax)
            else:
                print("--> No surface height information, cannot fill...")
        return

###################
#   Get methods   #
###################

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
        tsub = self.time['data'][(self.time['data'] >= dt_start) &
                                 (self.time['data'] <= dt_end)]
        datasub = data[(self.time['data'] >= dt_start) &
                       (self.time['data'] <= dt_end)]
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
