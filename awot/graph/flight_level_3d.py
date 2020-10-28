"""
awot.graph.flight_level_3d
==========================

A group of scripts to create 3D plots of flight level data.

TODO:
Fix the time_stamp module to work with savefig.  Okay on screen output,
but for some reason it bombs when trying to save.

Improve time series variable plotting?
"""


import numpy as np
import os
import matplotlib.pyplot as plt
#from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Line3D
from matplotlib.colors import Normalize
from matplotlib.dates import DateFormatter, date2num
from matplotlib.dates import (
    SecondLocator, MinuteLocator, HourLocator, DayLocator)
from matplotlib import ticker
from datetime import datetime
import scipy.ndimage as scim

from . import common
from .. import util


class FlightLevel3D(object):
    """Class for flight level plots."""

    def __init__(self, flightdata, basemap=None,
                 lon_name=None, lat_name=None,
                 alt_name=None, time_name=None,
                 uwind_name=None, vwind_name=None,
                 ):
        '''
        Intitialize the class to create plots.

        Parameters
        ----------
        flightdata : dict
            AWOT flight data dictionary instance.
        basemap : basemap instance
            Basemap instance to use for plotting.
        lon_name : str
            Variable name to use for longitude array.
            None uses AWOT default mapping.
        lat_name : str
            Variable name to use for latitude array.
            None uses AWOT default mapping.
        alt_name : str
            Variable name to use for altitude array.
            None uses AWOT default mapping.
        time_name : str
            Variable name to use for time array.
            None uses AWOT default mapping.
        uwind_name : str
            Variable name to use for zonal wind array.
            None uses AWOT default mapping.
        vwind_name : str
            Variable name to use for meridional wind array.
            None uses AWOT default mapping.
        '''
        if lon_name is None:
            lonkey = 'longitude'
        else:
            lonkey = lon_name
        if lat_name is None:
            latkey = 'latitude'
        else:
            latkey = lat_name
        if alt_name is None:
            altkey = 'altitude'
        else:
            altkey = alt_name
        if time_name is None:
            timekey = 'time'
        else:
            timekey = time_name
        if uwind_name is None:
            ukey = 'Uwind'
        else:
            ukey = uwind_name
        if vwind_name is None:
            vkey = 'Vwind'
        else:
            vkey = vwind_name
        self.longitude = flightdata[lonkey]
        self.latitude = flightdata[latkey]
        self.altitude = flightdata[altkey]
        self.flight_number = flightdata['flight_number']
        self.project = flightdata['project']
        self.platform = flightdata['platform']
        self.Uwind = flightdata[ukey]
        self.Vwind = flightdata[vkey]
        self.time = flightdata[timekey]
        self.flight_data = flightdata
        self.basemap = basemap
        common._check_basemap(self, strong=False)

        if self.basemap is not None:
            # Calculate x,y map position coordinates
            self.x, self.y = self.basemap(
                self.longitude['data'][:], self.latitude['data'][:])

#####################
#  Basemap methods  #
#####################


    def create_basemap(self, corners=None, proj=None, resolution='l',
                       area_thresh=None, lon_0=None, lat_0=None,
                       xpad=None, ypad=None,
                       meridians=True, parallels=True,
                       lon_spacing=2., lat_spacing=2.,
                       coastlines=True, countries=True, states=False,
                       counties=False, rivers=False, etopo=False):
        ''' Create a basemap instance. '''
        bm, ax = common.create_basemap_3d(corners=corners, proj=proj,
                                          resolution=resolution,
                                          area_thresh=area_thresh,
                                          lon_0=lon_0, lat_0=lat_0,
                                          meridians=meridians,
                                          parallels=parallels,
                                          lon_spacing=lon_spacing,
                                          lat_spacing=lat_spacing,
                                          xpad=xpad, ypad=ypad,
                                          coastlines=coastlines,
                                          countries=countries, states=states,
                                          counties=counties, rivers=rivers,
                                          etopo=etopo)
        self.basemap = bm
        self.x, self.y = self.basemap(
                self.longitude['data'][:], self.latitude['data'][:])
        return

#################
#  Track plots  #
#################

    def plot_trackmap_3d(self, track_cmap=None,
                         track_color='k', lw=1.5, ls=None, alpha=1.,
                         start_time=None, end_time=None,
                         min_altitude=None, max_altitude=None,
                         addlegend=False, legLoc='lower right', legLab=None,
                         addtitle=False, title=None, ax=None, fig=None,
#                         save_kmz=False, kmz_filepath=None, kmz_filename=None,
                         **kwargs):
        """
        Create a plot of the aircraft flight track, with optional enhancements.
        Note that the legend only works with single color track at this time.

        Parameters
        ----------
        track_cmap : Matplotlib colormap
            If color_by_altitude True, it will use this colormap.
        track_color : str
            If color_by_altitude False, this will be the track color.
        lw : float or int
            Line width to use for track.
        ls : str
            Line style string. Default to solid line.
        alpha : float
            Alpha value for transparency (0 = transparent, 1. = solid).
        start_time : str
            UTC time to use as start time for subsetting in datetime format.
            (e.g. 2014-08-20 12:30:00)
        end_time : str
            UTC time to use as an end time for subsetting in datetime format.
            (e.g. 2014-08-20 16:30:00)
        min_altitude : float
            Minimum altitude to use in mapping color.
        max_altitude : float
            Maximum altitude to use in mapping color.
        addlegend : bool
            Defaults to No legend drawn.
        legLoc : str
            See matplotlib legend object documentation.
        addtitle : bool
            Defaults to no title for plot.
        title : str
            See matplotlib axis object documentation.
        ax : Matplotlib axis instance
            Axis to plot on.
            None will use the current axis.
        fig : Matplotlib figure instance
            Figure whih to add the colorbar.
            None will use the current figure.

        **kwargs will pass specific arguments to subprograms,
          see Basemap and Matplotlib.
        """
        # parse parameters
        ax, fig = common._parse_ax_fig(ax, fig)

        # Check inputs
        if legLab is None:
            legLab = self.flight_number

        # Get start and end times (this deals with subsets)
        dt_start = self._get_datetime(start_time, get_start=True)
        dt_end = self._get_datetime(end_time, get_end=True)

        # Subset the data (will use min/max if None given)
        xsub, ysub = self._get_x_y_time_subset(
            start_time, end_time, return_time=False)
        lonsub, latsub = self._get_lon_lat_time_subset(start_time, end_time)
        timeSub, altsub = self._get_time_var_time_subset(
            'altitude', start_time, end_time)

        # Clean up the masked data for plotting
        all_mask = latsub.mask + lonsub.mask + altsub.mask
        lonmask = np.ma.array(lonsub, mask=all_mask).compressed()
        latmask = np.ma.array(latsub, mask=all_mask).compressed()
        altmask = np.ma.array(altsub, mask=all_mask).compressed()

        xmask, ymask = self.basemap(lonmask, latmask)

        # Set a couple of keyword defaults for KMZ option
        show_legend, legend_label = False, ' '

        # Plot the track either coloring by altitude or as single color
        if track_cmap is None:
            track_cmap = plt.get_cmap()
        if ls is None:
            ls = '-'
        p = ax.plot(xmask, ymask, altmask, color=track_color,
                    linewidth=lw, linestyle=ls, alpha=alpha)

        # Label the axes
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_zlabel('Altitude (m)')
        if min_altitude is None:
            min_altitude = np.min(altmask)
        if max_altitude is None:
            max_altitude = np.max(altmask)
        ax.set_xlim3d(np.min(xmask), np.max(xmask))
        ax.set_ylim3d(np.min(ymask), np.max(ymask))
        ax.set_zlim3d(min_altitude, max_altitude)

        # Save the KMZ file if requested
#         if save_kmz:
#             lonrange = (np.min(lonmask), np.max(lonmask))
#             latrange = (np.min(latmask), np.max(latmask))
#             times = [dt_start, dt_end]
#             if kmz_filepath is None:
#                 kmz_filepath = os.getcwd()
#             if kmz_filename is None:
#                 kmz_filename = ('awot_' + self.platform + '_' +
#                                 self.flight_number + '_altitude' + '.kmz')
#             util.write_kmz.write_kmz(fig, ax, p, lonrange, latrange, times,
#                                      file_path=kmz_filepath,
#                                      file_name=kmz_filename,
#                                      show_legend=show_legend,
#                                      legend_label=legend_label)
        if addlegend:
            self.basemap.ax.legend(loc=legLoc, fontsize=11, frameon=True,
                                   handlelength=2, handletextpad=1,
                                   labelspacing=.2)

        if addtitle:
            if title is None:
                title = self.project + ' ' + self.platform
            self.basemap.ax.set_title(title)

    def plot_variable3d(self, field, track_cmap=None,
                        track_lw=1.5, alpha=1.,
                        start_time=None, end_time=None,
                        min_altitude=None, max_altitude=None,
                        min_value=None, max_value=None,
                        add_cb=True, cborientation='vertical',
                        cbshrink=0.5, cbaspect=20,
                        cblabel=None, addlegend=False, legLoc='lower right',
                        legLab=None, addtitle=False, title=None,
                        ax=None, fig=None,
#                        save_kmz=False, kmz_filepath=None, kmz_filename=None,
                        **kwargs):
        """
        Create a plot of the aircraft flight track, with optional enhancements.
        Note that the legend only works with single color track at this time.

        Parameters
        ----------
        field : float
            Variable passed to plot, defaults to temperature.
        track_cmap : Matplotlib colormap instance
            Colormap to use for coding values.
        track_lw : float or int
            Line width to use for track.
        alpha : float
            Alpha value for transparency (0 = transparent, 1. = solid).
        start_time : str
            UTC time to use as start time for subsetting in datetime format.
            (e.g. 2014-08-20 12:30:00)
        end_time : str
            UTC time to use as an end time for subsetting in datetime format.
            (e.g. 2014-08-20 16:30:00)
        min_altitude : float
            Minimum altitude to use in mapping color.
        max_altitude : float
            Maximum altitude to use in mapping color.
        min_value : float
            Minimum value to use in mapping color.
        max_value : float
            Maximum value to use in mapping color.
        add_cb : bool
            True plots a colorbar, Fasle does not.
        cborientation : str
            Orientation of colorbar, options: 'vertical' or 'horizontal'.
        cbshrink : float
            Fractional size to shrink colorbar.
        cbaspect : float
            Ratio between height and width of colorbar.
        addlegend : bool
            Defaults to No legend drawn.
        legLoc : str
            See matplotlib legend object documentation.
        addtitle : bool
            Defaults to no title for plot.
        title : str
            See matplotlib axis object documentation.
        ax : Matplotlib axis instance
            Axis to plot. None will use the current axis.
        fig : Matplotlib figure instance
            Figure which to add the colorbar.
            None will use the current figure.

        **kwargs will pass specific arguments to subprograms,
          see Basemap and Matplotlib.
        """

        # Check to see if field exists
        common._check_field(self.flight_data, field)

        # parse parameters
        ax, fig = common._parse_ax_fig(ax, fig)

        # Get start and end times (this deals with subsets)
        dt_start = self._get_datetime(start_time, get_start=True)
        dt_end = self._get_datetime(end_time, get_end=True)

        # Subset the data (will use min/max if None given)
        xsub, ysub = self._get_x_y_time_subset(
            start_time, end_time, return_time=False)
        lonsub, latsub = self._get_lon_lat_time_subset(start_time, end_time)
        altsub = self._get_var_time_subset(
            'altitude', start_time, end_time)
        varsub = self._get_var_time_subset(
            field, start_time, end_time)

        # Clean up the masked data for plotting
        all_mask = latsub.mask + lonsub.mask + varsub.mask
        lonmask = np.ma.array(lonsub, mask=all_mask).compressed()
        latmask = np.ma.array(latsub, mask=all_mask).compressed()
        altmask = np.ma.array(altsub, mask=all_mask).compressed()
        varmask = np.ma.array(varsub, mask=all_mask).compressed()

        xmask, ymask = self.basemap(lonmask, latmask)

        # Check inputs
        if min_value is None:
            min_value = varmask.min()

        if max_value is None:
            max_value = varmask.max()

        if legLab is None:
            legLab = self.flight_number

        # Set a couple of keyword defaults for KMZ option
        show_legend, legend_label = False, ' '

        # Plot the track
        if track_cmap is None:
            track_cmap = plt.get_cmap()
        p = self._colorline(xmask, ymask, altmask, varmask, cmap=track_cmap,
                            vmin=min_value, vmax=max_value,
                            linewidth=track_lw, alpha=alpha)
        ax.add_collection3d(p, zdir='z')

        # Label the axes
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_zlabel('Altitude (m)')
        if min_altitude is None:
            min_altitude = np.min(altmask)
        if max_altitude is None:
            max_altitude = np.max(altmask)
        ax.set_zlim(min_altitude, max_altitude)

        if add_cb:
            cb = fig.colorbar(p, orientation=cborientation,
                              shrink=cbshrink, aspect=cbaspect)
            cb.set_label(cblabel)
            show_legend = True
            legend_label = cblabel

        # Save the KMZ file if requested
#         if save_kmz:
#             lonrange = (np.min(lonmask), np.max(lonmask))
#             latrange = (np.min(latmask), np.max(latmask))
#             times = [dt_start, dt_end]
#             if kmz_filepath is None:
#                 kmz_filepath = os.getcwd()
#             if kmz_filename is None:
#                 kmz_filename = ('awot_' + self.platform + '_' +
#                                 self.flight_number + '_' + cblabel +
#                                 '_' + '.kmz')
#             util.write_kmz.write_kmz(fig, ax, p, lonrange, latrange, times,
#                                      file_path=kmz_filepath,
#                                      file_name=kmz_filename,
#                                      show_legend=show_legend,
#                                      legend_label=legend_label)

        if addlegend:
            self.basemap.ax.legend(
                loc=legLoc, fontsize=11, frameon=True,
                handlelength=2, handletextpad=1, labelspacing=.2)

        if addtitle:
            if title is None:
                title = self.project + ' ' + self.platform
            ax.set_title(title)

    def plot_radar_cross_section(
            self, radar, field, mask_procedure=None, mask_tuple=None,
            start_time=None, end_time=None, title=" ", title_size=20,
            cminmax=(0., 60.), clevs=25, vmin=15., vmax=60.,
            cmap='gist_ncar', clabel='dBZ',
            color_bar=True, cb_orient='vertical', cb_pad=.05, cb_tick_int=2,
            cb_fontsize=None, cb_ticklabel_size=None,
            x_axis_array='distance', date_format='%H:%M', tz=None,
            date_minor_string='minute', plot_km=False,
            ax=None, fig=None):
        '''
        Plot a cross-section along a flight path segment.

        Parameters
        ----------
        radar : AWOT radar object
            RadarHorizontalPlot object.
        field : str
            3-D variable inf RadarHorizontalPlot object
            (e.g. Reflectivity ['dBZ']) to use in plot.
        start_pt, end_pt : tuple
            (lat, lon) Tuple of start, end points for cross-section.
        xs_length : int
            Number of
        start_time : string
            UTC time to use as start time for subsetting in datetime format.
            (e.g. 2014-08-20 12:30:00)
        end_time : string
            UTC time to use as an end time for subsetting in datetime format.
            (e.g. 2014-08-20 16:30:00)
        mask_procedure : str
            String indicating how to apply mask via numpy, possibilities are:
            'less', 'less_equal', 'greater', 'greater_equal', 'equal',
            'inside', 'outside'.
        mask_tuple : (str, float[, float])
            Tuple containing the field name and value(s) below which to mask
            field prior to plotting.
        cminmax : tuple
            (min,max) values for controur levels.
        clevs : integer
            Number of contour levels.
        vmin : float
            Minimum contour value to display.
        vmax : float
            Maximum contour value to display.
        cmap : string
            Matplotlib color map to use.
        clabel : string
            Label for colorbar (e.g. units 'dBZ').
        title : string
            Plot title.
        title_size : int
            Font size of title to display.
        color_bar : boolean
            True to add colorbar, False does not.
        cb_fontsize : int
            Font size of the colorbar label.
        cb_ticklabel_size : int
            Font size of colorbar tick labels.
        cb_pad : str
            Pad to move colorbar, in the form "5%", pos is to right
            for righthand location.
        cb_orient : str
            Colorbar orientation, either 'vertical' or 'horizontal'.
        cb_tick_int : int
            Interval to use for colorbar tick labels,
            higher number "thins" labels.
        x_axis_array: str
            X-axis array to plot against either 'distance' [default] or 'time'.
        date_format : str
            Format of the time string for x-axis labels.
        tz : str
            Time zone info to use when creating axis labels (see datetime).
        date_minor_string : str
            Sting to set minor ticks of date axis,
            'second','minute','hour','day' supported.
        plot_km : boolean
            True to convert meters to kilometers for cappi_height. False
            retains meters information.
        ax : Matplotlib axis instance
            Axis to plot on. None will use the current axis.
        fig : Matplotlib figure instance
            Figure which to add the plot. None will use the current figure.
        '''
        # parse parameters
        ax, fig = common._parse_ax_fig(ax, fig)

        # Get start and end times (this deals with subsets)
        dt_start = self._get_datetime(start_time, get_start=True)
        dt_end = self._get_datetime(end_time, get_end=True)

        # Subset the data (will use min/max if None given)
        xsub, ysub = self._get_x_y_time_subset(
            start_time, end_time, return_time=False)
        lonsub, latsub = self._get_lon_lat_time_subset(start_time, end_time)
        timeSub = self._get_time_subset(start_time, end_time)

        # Return masked or unmasked variable
        Var, Data = self._get_radar_variable_dict_data(radar, field)
        if mask_procedure is not None:
            Data = common.get_masked_data(Data, mask_procedure, mask_tuple)

        # Create contour level array
        clevels = np.linspace(cminmax[0], cminmax[1], clevs)

        # Create an array to hold the interpolated cross-section
        xs_data = np.empty([len(lonsub), len(radar['height']['data'][:])])

        # Create arrays for cross-section lon-lat index points
        xsY = np.empty(len(lonsub))
        xsX = np.empty(len(lonsub))
        Xdist = np.empty(len(lonsub))
        Ydist = np.empty(len(lonsub))
        xsDist = np.empty(len(lonsub))

        for ii in range(len(lonsub)):
            xsX[ii] = self._get_lon_index(lonsub[ii], radar)
            xsY[ii] = self._get_lat_index(latsub[ii], radar)

            # Calculate the distance array along the cross-section
            # Need to keep a running tally moving through track array
            if ii == 0:
                Xdist[ii] = 0.
                Ydist[ii] = 0.
                xsDist[ii] = 0.
            else:
                Xdist[ii] = np.absolute(
                    (np.pi * common.EARTH_RADIUS / 180.) *
                    (lonsub[ii] - lonsub[ii - 1]))
                Ydist[ii] = np.absolute(
                    (np.pi * common.EARTH_RADIUS / 180.) *
                    (latsub[ii] - latsub[ii - 1]))
                xsDist[ii] = (np.sqrt(Xdist[ii]**2 + Ydist[ii]**2)
                              ) + xsDist[ii - 1]

        # Loop through each level to create cross-section and stack them
        for nlev in range(len(radar['height']['data'][:])):
            # Extract the values along the line, using cubic interpolation
            xs_data[:, nlev] = scim.map_coordinates(
                Data[nlev, :, :], np.vstack((xsY, xsX)), prefilter=False)
            # , mode='nearest')

        # Calculate the distance array along the cross-section
        if x_axis_array == 'distance':
            Xax = xsDist
        elif x_axis_array == 'time':
            Xax = date2num(timeSub)

        # Convert Height, distance arrays to 2D
        if plot_km:
            Ht2D, Xax2D = np.meshgrid(radar['height']['data'][:]/1000., Xax)
            ylabel = 'Altitude (km)'
        else:
            Ht2D, Xax2D = np.meshgrid(radar['height']['data'][:], Xax)
            ylabel = 'Altitude (m)'

        p = ax.pcolormesh(Xax2D, Ht2D, np.ma.masked_less_equal(xs_data, -800.),
                          vmin=vmin, vmax=vmax, cmap=cmap)

        ax.set_ylabel(ylabel)
        if x_axis_array == 'distance':
            ax.set_xlabel('Distance along track (km)')
        elif x_axis_array == 'time':
            ax.xaxis_date()
            # Set the date format
            date_Fmt = DateFormatter(date_format, tz=tz)
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

        # Add title
        ax.set_title(title, fontsize=title_size)

        # Add Colorbar
        if color_bar:
            cbStr = Var['long_name'] + ' (' + Var['units'] + ')'
            cb = common.add_colorbar(ax, p, orientation=cb_orient, pad=cb_pad,
                                     label=cbStr, fontsize=cb_fontsize,
                                     ticklabel_size=cb_ticklabel_size,
                                     clevs=clevs, tick_interval=cb_tick_int)

        # Add title
        ax.set_title(title, fontsize=title_size)

############################
#  Basemap add-on methods  #
############################

    def draw_boundary(self, **kwargs):
        """
        Draw a boundary around the map.

        Parameters
        ----------
        See basemap documentation.
        """
        self.basemap.drawmapboundary(**kwargs)

    def draw_scale(self, location='lower left', lon=None, lat=None,
                   lon0=None, lat0=None, length=None, **kwargs):
        """
        Draw a reference scale on the map.

        Parameters
        ----------
        location : str
           Where to put the scale, default is 'lower left'.
           Options are 'lower left', 'upper left', 'lower right',
           'upper right', 'upper middle', 'lower middle'.

        All other parameters, See basemap documentation.

        Note that the location parameter can be overridden by setting
        the lon and lat variables explicitly.
        """
        # Check on where to put the label, defaults to the else statement
        # that calculates 'lower left'
        if (location.lower() == 'upper left') or \
                (location.lower() == 'upper_left'):
            xScaleOffset = (
                self.basemap.urcrnrlon - self.basemap.llcrnrlon) * 0.10
            yScaleOffset = (
                self.basemap.urcrnrlat - self.basemap.llcrnrlat) * 0.90
        elif (location.lower() == 'upper middle') or \
                (location.lower() == 'upper_middle'):
            xScaleOffset = (
                self.basemap.urcrnrlon - self.basemap.llcrnrlon) * 0.50
            yScaleOffset = (
                self.basemap.urcrnrlat - self.basemap.llcrnrlat) * 0.90
        elif (location.lower() == 'upper right') or \
                (location.lower() == 'upper_right'):
            xScaleOffset = (
                self.basemap.urcrnrlon - self.basemap.llcrnrlon) * 0.90
            yScaleOffset = (
                self.basemap.urcrnrlat - self.basemap.llcrnrlat) * 0.90
        elif (location.lower() == 'lower right') or \
                (location.lower() == 'lower_right'):
            xScaleOffset = (
                self.basemap.urcrnrlon - self.basemap.llcrnrlon) * 0.90
            yScaleOffset = (
                self.basemap.urcrnrlat - self.basemap.llcrnrlat) * 0.10
        elif (location.lower() == 'lower middle') or \
                (location.lower() == 'lower_middle'):
            xScaleOffset = (
                self.basemap.urcrnrlon - self.basemap.llcrnrlon) * 0.50
            yScaleOffset = (
                self.basemap.urcrnrlat - self.basemap.llcrnrlat) * 0.10
        else:
            xScaleOffset = (
                self.basemap.urcrnrlon - self.basemap.llcrnrlon) * 0.10
            yScaleOffset = (
                self.basemap.urcrnrlat - self.basemap.llcrnrlat) * 0.10

        if lon is None:
            lon = self.basemap.llcrnrlon + xScaleOffset
        if lat is None:
            lat = self.basemap.llcrnrlat + yScaleOffset
        if lon0 is None:
            lon0 = (self.basemap.llcrnrlon + self.basemap.urcrnrlon) / 2.
        if lat0 is None:
            lat0 = (self.basemap.llcrnrlat + self.basemap.urcrnrlat) / 2.
        if length is None:
            length = 100

        self.basemap.drawmapscale(lon=lon, lat=lat, lon0=lon0, lat0=lat0,
                                  length=length, barstyle='fancy',
                                  labelstyle='simple', **kwargs)

    def draw_barbs(self, barbspacing=10, barbcolor='k', flagcolor='k',
                   blength=8, lw=0.5, start_time=None, end_time=None,
                   **kwargs):
        """
        Draw a reference scale on the map.

        Parameters
        ----------
        barbspacing : int
           Plot every nth via the barbspacing provided.
        See basemap documentation.
        """
        # Get start and end times (this deals with subsets)
        dt_start = self._get_datetime(start_time, get_start=True)
        dt_end = self._get_datetime(end_time, get_end=True)

        # Subset the data
        X = self.x[(self.time['data'][:] >= dt_start) &
                   (self.time['data'][:] <= dt_end)]
        Y = self.y[(self.time['data'][:] >= dt_start) &
                   (self.time['data'][:] <= dt_end)]
        Uwnd = self.Uwind['data'][(self.time['data'][:] >= dt_start) &
                                  (self.time['data'][:] <= dt_end)]
        Vwnd = self.Vwind['data'][(self.time['data'][:] >= dt_start) &
                                  (self.time['data'][:] <= dt_end)]

        # Only plot every nth barb from the barbspacing parameter
        self.basemap.barbs(X[::barbspacing], Y[::barbspacing],
                           Uwnd[::barbspacing], Vwnd[::barbspacing],
                           barbcolor=barbcolor, flagcolor=flagcolor,
                           linewidth=lw, **kwargs)

    def plot_point(self, lon, lat, symbol='ro', label_text=None,
                   label_offset=(None, None), text_size=None, **kwargs):
        """
        Plot a point on the current map.

        Additional arguments are passed to basemap.plot.

        Parameters
        ----------
        lon : float
            Longitude of point to plot.
        lat : float
            Latitude of point to plot.
        symbol : str
            Matplotlib compatible string which specified the symbol of the
            point.
        label_text : str
            Optional text to label symbol with. If None no label will be added.
        label_offset : [float, float]
            Offset in lon, lat degrees for the bottom left corner of the label
            text relative to the point. A value of None will use 0.01 default.
        """
        lon_offset, lat_offset = label_offset
        if lon_offset is None:
            lon_offset = 0.01
        if lat_offset is None:
            lat_offset = 0.01
        if text_size is None:
            text_size = 12

        # Plot the symbol
        self.basemap.plot(lon, lat, symbol, latlon=True, **kwargs)
        # Attach the text
        if label_text is not None:
            # basemap does not have a text method so we must determine
            # the x and y points and plot them on the basemap's axis.
            x_text, y_text = self.basemap(lon + lon_offset,
                                          lat + lat_offset)
            self.basemap.ax.text(x_text, y_text, label_text, size=text_size)

    def plot_line_geo(self, line_lons, line_lats, line_style='r-', **kwargs):
        """
        Plot a line segments on the current map given values in lat and lon.

        Additional arguments are passed to basemap.plot.

        Parameters
        ----------
        line_lons : array
            Longitude of line segment to plot.
        line_lats : array
            Latitude of line segment to plot.
        line_style : str
            Matplotlib compatible string which specifies the line style.
        """
        X, Y = self.basemap(line_lons, line_lats)
        self.basemap.plot(X, Y, line_style, **kwargs)

    def time_stamps(self, labelspacing=1800, symbol='k*',
                    size=12, color='k', label_offset=(None, None),
                    start_time=None, end_time=None, ax=None, **kwargs):
        """
        Add time stamps along track according to labelspacing.

        Parameters
        ----------
        labelspacing : int
            Plot every nth via the labelpacing provided.
            The default assumes seconds and labels every 30 minutes.
        symbol : str
            Color-symbol combo per matplotlib.pyplot.plot
        size : int
            Font size used.
        color : str
            Matplotlib color.
        label_offset : [float, float]
            Offset in lon, lat degrees for the bottom left corner of the label
            text relative to the point. A value of None will use 0.01 default.
        start_time : str
            UTC time to use as start time for subsetting in datetime format.
            (e.g. 2014-08-20 12:30:00)
        end_time : str
            UTC time to use as an end time for subsetting in datetime format.
            (e.g. 2014-08-20 16:30:00)
        ax : Matplotlib axis instance
            Axis on which to plot.
        """
        # parse parameters
        ax = common._parse_ax(ax)

        # Only plot every nth barb from the barbspacing parameter
        lon_offset, lat_offset = label_offset
        if lon_offset is None:
            lon_offset = 0.01
        if lat_offset is None:
            lat_offset = 0.01

        # Get start and end times (this deals with subsets)
        dt_start = self._get_datetime(start_time, get_start=True)
        dt_end = self._get_datetime(end_time, get_end=True)

        time = self.time['data'][(self.time['data'][:] >= dt_start) &
                                 (self.time['data'][:] <= dt_end)]
        lon = self.longitude['data'][(self.time['data'][:] >= dt_start) &
                                     (self.time['data'][:] <= dt_end)]
        lat = self.latitude['data'][(self.time['data'][:] >= dt_start) &
                                    (self.time['data'][:] <= dt_end)]

        # Find the number of points to plot
        labeltimes = time[::labelspacing]
        lons = lon[::labelspacing]
        lats = lat[::labelspacing]
        xpos_mrkr, ypos_mrkr = self.basemap(lons, lats)
        xpos_text, ypos_text = self.basemap(lons + lon_offset,
                                            lats + lat_offset)

        for nn in range(len(labeltimes)):
            # Check to make sure there are no missing lon/lat values and
            # only plot the valid points
            # This causes fatal error when trying to save figure later
            if np.logical_and(np.isfinite(xpos_mrkr[nn]),
                              np.isfinite(ypos_mrkr[nn])):

                # Plot the symbol
                ax.plot(xpos_mrkr[nn], ypos_mrkr[nn], symbol, **kwargs)
                label_text = labeltimes[nn].strftime("%H:%M")
                # Attach the text
                ax.text(xpos_text[nn], ypos_text[nn], label_text,
                        fontsize=size, color=color)

    def label_time_point(self, label_time, label_text=None,
                         add_symbol=False, symbol='k*', size=12, color='k',
                         label_offset=(None, None), ax=None, **kwargs):
        """
        Add a label at a particular time along track.

        Parameters
        ----------
        label_time : str
            UTC time to use as start time for subsetting in datetime format.
            (e.g. 2014-08-20 12:30:00)
        label_text : str
            Optional label to use, if none given a time stamp is generated.
        add_symbol : bool
            True to add a symbol at the label_time, provided by symbol keyword.
        symbol : str
            Color-symbol combo per matplotlib.pyplot.plot
        size : int
            Font size used.
        color : str
            Matplotlib color
        label_offset : [float, float]
            Offset in lon, lat degrees for the bottom left corner of the label
            text relative to the point. A value of None will use 0.01 default.
        ax : Matplotlib axis instance
            Axis on which to plot.
        """

        # parse parameters
        ax = common._parse_ax(ax)

        # Only plot every nth barb from the barbspacing parameter
        lon_offset, lat_offset = label_offset
        if lon_offset is None:
            lon_offset = 0.01
        if lat_offset is None:
            lat_offset = 0.01

        # Get start and end times (this deals with subsets)
        dt_point = self._get_datetime(label_time)

        time = self.time['data'][(self.time['data'][:] >= dt_point)]
        lon = self.longitude['data'][(self.time['data'][:] >= dt_point)]
        lat = self.latitude['data'][(self.time['data'][:] >= dt_point)]

        xpos_mrkr, ypos_mrkr = self.basemap(lon[0], lat[0])
        xpos_text, ypos_text = self.basemap(
            lon[0] + lon_offset, lat[0] + lat_offset)

        # Check to make sure there are no missing lon/lat values and
        # only plot the valid points
        # This causes fatal error when trying to save figure later
        if np.logical_and(np.isfinite(xpos_mrkr), np.isfinite(ypos_mrkr)):
            # Plot the symbol
            if add_symbol:
                ax.plot(xpos_mrkr, ypos_mrkr, symbol, **kwargs)
            if label_text is None:
                label_text = time[0].strftime("%H:%M")

            # Attach the text
            ax.text(xpos_text, ypos_text, label_text,
                    fontsize=size, color=color)
        return

##################
#  Line methods  #
##################
    # Original Multiline color plotting found at
    # http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    # Modified for this class#

    def _colorline(self, x, y, z, data, cmap=None, norm=None,
                   linewidth=1.5, alpha=0.5, vmin=None, vmax=None):
        '''
        Plot a colored line with coordinates x and y.
        Optionally specify colors in the array z.
        Optionally specify a colormap, a norm function and a line width.

        Parameters
        ----------
        x : float array
            X coordinates to map.
        Y : float array
            Y coordinates to map.
        z : float array
            Variable data to plot.
        cmap : str
            String of colormap to use.
        norm : float array
            Normalized array for colormapping.
        linewidth : float
            Sets the width of the line drawn.
        alpha : float
            Sets the transparency of the line (0 = transparent, 1 = opaque).
        vmin : float
            Sets the minimum to scale colormap luminance.
        vmax : float
            Sets the maximum to scale colormap luminance.
        '''
        # Get the colormap to work with
        if cmap is None:
            cmap = 'spectral'

        if vmin is None:
            vmin = data.min()

        if vmax is None:
            vmax = data.max()

        # Use the vmin/vmax as a mask
        data = np.ma.masked_outside(data, vmin, vmax)

        # Set the normalization to the min and max
        if norm is None:
            norm = plt.Normalize(vmin, vmax)

        # Create list of line segments from x and y coordinates,
        # in the correct format for LineCollection:
        # an array of the form numlines x (points per line) x 2 (x & y) array
        points = np.array([x, y, z]).T.reshape(-1, 1, 3)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        lc = Line3DCollection(segments, cmap=cmap, norm=norm,
                              linewidth=linewidth, alpha=alpha)
        lc.set_array(data)

        return lc

    # NOTE: This is in here for legacy, RECOMMENDED NOT TO USE
    def draw_marker_track(self, cmap=None, msize=2, mstyle='.'):
        """
        Draws a track by plotting each marker and colored based on height
        This is REALLY slow for large number of points.
        Left in for reference.

        Highly recommended to not use this function.
        """
        # Color code the marker based on height
        def get_marker_color(altitude):
            if ((altitude > 25.) & (altitude <= 150.)):
                mcol = 'm'
                return mcol
            elif ((altitude > 150.) & (altitude <= 460.)):
                mcol = 'b'
                return mcol
            elif (altitude > 460.) & (altitude <= 915.):
                mcol = 'g'
                return mcol
            elif (altitude > 915.) & (altitude <= 1525.):
                mcol = 'c'
                return mcol
            elif (altitude > 1525.) & (altitude <= 2290.):
                mcol = 'y'
                return mcol
            elif (altitude > 2290.) & (altitude <= 3050.):
                mcol = 'r'
                return mcol
            elif (altitude > 3050.):
                mcol = 'k'
                return mcol
            else:
                mcol = 'w'
                return mcol

        # Get the colormap to work with
        if cmap is None:
            cmap = 'spectral'

        # Loop through each point.
        # Unfortunately as of Aug 19, 2014 basemap.scatter
        # was not working.  No display and I have no idea why!
        # This is SLOW and not recommended.
        for lon, lat, alt in zip(
                self.longitude['data'][:].compressed(),
                self.latitude['data'][:].compressed(),
                self.altitude['data'][:].compressed()):
            mcol = get_marker_color(alt)
            self.basemap.plot(lon, lat, color=mcol, latlon=True,
                              marker=mstyle, markersize=msize)

#################
#  Get methods  #
#################

    def _get_datetime(self, time_string, get_start=False, get_end=False):
        '''Get a start time as datetime instance for subsetting.'''
        # Check to see if time is subsetted
        if time_string is None:
                if get_start is True:
                    dt = self.time['data'][:].min()
                if get_end is True:
                    dt = self.time['data'][:].max()
        else:
            tStr = [time_string[0:4], time_string[5:7], time_string[8:10],
                    time_string[11:13], time_string[14:16],
                    time_string[17:19], '0']
            tInt = [int(s) for s in tStr]
            try:
                dt = datetime(tInt[0], tInt[1], tInt[2], tInt[3],
                              tInt[4], tInt[5], tInt[6])
            except:
                import warnings
                warnings.warn(common.DATE_STRING_FORMAT)
                return

        return dt

    def _get_start_datetime(self, start_time):
        '''
        Get a start time as datetime instance for subsetting.
        '''
        # Check to see if time is subsetted
        if start_time is None:
            dt_start = self.time['data'][:].min()
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

        return dt_start

    def _get_end_datetime(self, end_time):
        '''
        Get a start time as datetime instance for subsetting.
        '''
        # Check to see if the time is subsetted
        if end_time is None:
            dt_end = self.time['data'][:].max()
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

        return dt_end

    def _get_x_y_time_subset(self, start_time, end_time, return_time=False):
        '''
        Get a subsetted X and Y to control track length if input by user.
        '''
        # Check to see if time is subsetted
        dt_start = self._get_datetime(start_time, get_start=True)
        dt_end = self._get_datetime(end_time, get_end=True)

        x = self.x[(self.time['data'][:] >= dt_start) &
                   (self.time['data'][:] <= dt_end)]
        y = self.y[(self.time['data'][:] >= dt_start) &
                   (self.time['data'][:] <= dt_end)]

        if return_time:
            time = self.time['data'][(self.time['data'][:] >= dt_start) &
                                     (self.time['data'][:] <= dt_end)]

            return x, y, time
        else:

            return x, y

    def _get_lon_lat_time_subset(self, start_time, end_time):
        '''
        Get a subsetted Lon and Lat to control track length if input by user.
        '''
        # Check to see if time is subsetted
        dt_start = self._get_datetime(start_time, get_start=True)
        dt_end = self._get_datetime(end_time, get_end=True)

        lon = self.longitude['data'][(self.time['data'][:] >= dt_start) &
                                     (self.time['data'][:] <= dt_end)]
        lat = self.latitude['data'][(self.time['data'][:] >= dt_start) &
                                    (self.time['data'][:] <= dt_end)]

        return lon, lat

    def _get_var_time_subset(self, field, start_time, end_time):
        '''
        Get a subsetted time and Variable to control track
        length if input by user.
        '''
        # Check to see if time is subsetted
        dt_start = self._get_datetime(start_time, get_start=True)
        dt_end = self._get_datetime(end_time, get_end=True)
        var = self.flight_data[field]
        vsub = var['data'][(self.time['data'][:] >= dt_start) &
                           (self.time['data'][:] <= dt_end)]
        return vsub

    def _get_time_var_time_subset(self, field, start_time, end_time):
        '''
        Get a subsetted time and Variable to control track
        length if input by user.
        '''
        # Check to see if time is subsetted
        dt_start = self._get_datetime(start_time, get_start=True)
        dt_end = self._get_datetime(end_time, get_end=True)
        tsub = self.time['data'][(self.time['data'][:] >= dt_start) &
                                 (self.time['data'][:] <= dt_end)]
        var = self.flight_data[field]
        vsub = var['data'][(self.time['data'][:] >= dt_start) &
                           (self.time['data'][:] <= dt_end)]
        return tsub, vsub

    def _get_time_subset(self, start_time, end_time):
        '''
        Get a subsetted time to control track length if input by user.
        '''
        # Check to see if time is subsetted
        dt_start = self._get_datetime(start_time, get_start=True)
        dt_end = self._get_datetime(end_time, get_end=True)
        tsub = self.time['data'][(self.time['data'][:] >= dt_start) &
                                 (self.time['data'][:] <= dt_end)]
        return tsub

    def _get_lat_index(self, value, radar):
        ''' Calculate the exact index position within latitude array. '''
        # Find the spacing
        dp = radar['latitude']['data'][1] - radar['latitude']['data'][0]

        # Calculate the relative position
        pos = (value - radar['latitude']['data'][0]) / dp
        return pos

    def _get_lon_index(self, value, radar):
        ''' Calculate the exact index position within longitude array. '''
        # Find the spacing
        dp = radar['longitude']['data'][1] - radar['longitude']['data'][0]

        # Calculate the relative position
        pos = (value - radar['longitude']['data'][0]) / dp
        return pos

    def _get_radar_variable_dict_data(self, radar, field):
        ''' Get the variable from the fields dictionary. '''
        Var, data = radar['fields'][field], radar['fields'][field]['data'][:]
        return Var, data

    def _get_var_dict(self, field):
        ''' Get the variable from the flight dictionary. '''
        var, data = self.flight_data[field], self.flight_data[field]['data'][:]
        return var, data
