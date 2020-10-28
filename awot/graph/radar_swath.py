"""
awot.graph.radar_swath
===========================

A group of scripts create plots of horizontal swaths from airborne radar.

"""

from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import ticker
import numpy as np
import scipy.ndimage as scim

from .common import (find_nearest_indices, get_masked_data,
                     _check_basemap, create_basemap,
                     _get_start_datetime, _get_end_datetime,
                     _parse_ax_fig, _parse_ax)


class RadarSwathPlot(object):
    """Class to plot a horizontal radar image."""

    def __init__(self, radar, basemap=None,
                 lon_name=None, lat_name=None, height_name=None,
                 range_name=None, time_name=None,
                 resolution='l', projection='lcc', area_thresh=None,
                 ax=None):
        '''
        Parameters
        ----------
        radar : dict
            AWOT radar object.
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

        Initialize the class to create plots.
        '''

        # Save the airborne class to this class in case cross-section is passed
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
        if range_name is None:
            self.range = self.radar['range']
        else:
            self.range = self.radar[range_name]
        if time_name is None:
            self.time = self.radar['time']
        else:
            self.time = self.radar[time_name]

        if height_name is None:
            alt2D, range2D = np.meshgrid(
                                  self.radar['altitude']['data'],
                                  self.range['data'])
            self.height = {
                           'data': alt2D - range2D,
                           'units': self.radar['altitude']['units']
                           }
        else:
            self.height = self.radar[height_name]

        if basemap is not None:
            self.basemap = basemap
        else:
            # Create Basemap instance
            lon_0 = np.median(self.longitude['data'][:])
            lat_0 = np.median(self.latitude['data'][:])
            corners = [np.min(self.longitude['data'][:]),
                       np.min(self.latitude['data'][:]),
                       np.max(self.longitude['data'][:]),
                       np.max(self.latitude['data'][:])]
            self.basemap = create_basemap(
                             proj=projection,
                             lon_0=lon_0, lat_0=lat_0,
                             corners=corners, ax=ax,
                             resolution=resolution,
                             area_thresh=area_thresh)

####################
#   Plot modules  ##
####################

    def plot_at_range_height(self, field, range_height, mask_procedure=None,
                             mask_tuple=None, cminmax=(0., 60.),
                             clevs=25, vmin=15., vmax=60.,
                             clabel='dBZ', title=" ", title_size=20,
                             cmap='gist_ncar', color_bar=True, cb_pad="5%",
                             cb_loc='right', cb_tick_int=2,
                             start_time=None, end_time=None,
                             ax=None, fig=None):
        """
        Produce a plot of radar swath at the give range height.

        Parameters
        ----------
        field : str
            3-D variable (e.g. Reflectivity [dBZ]) to use in plot.
        range_height : float
            Height in kilometers at which to plot the horizontal field.
        mask_procedure : str
            String indicating how to apply mask via numpy, possibilities are:
            'less', 'less_equal', 'greater', 'greater_equal',
            'equal', 'inside', 'outside'.
        mask_tuple : (str, float[, float])
            Tuple containing the field name and value(s) below which to mask
            field prior to plotting, for example to mask all data where.
        cminmax : tuple
            (min,max) values for controur levels.
        clevs : int
            Number of contour levels.
        vmin : float
            Minimum contour value to display.
        vmax : float
            Maximum contour value to display.
        clabel : str
            Label for colorbar (e.g. units 'dBZ').
        title : str
            Plot title.
        title_size : int
            Font size of title to display.
        cmap : str
            Matplotlib color map to use.
        color_bar : boolean
            True to add colorbar, False does not.
        cb_pad : str
            Pad to move colorbar, in the form "5%", pos is to
            right for righthand location.
        cb_loc : str
            Location of colorbar, default is 'right', also available:
            'bottom', 'top', 'left'.
        cb_tick_int : int
            Interval to use for colorbar tick labels,
            higher number "thins" labels.
        start_time : str
            UTC time to use as start time for subsetting in datetime format.
            (e.g. 2014-08-20 12:30:00)
        end_time : str
            UTC time to use as an end time for subsetting in datetime format.
            (e.g. 2014-08-20 16:30:00)
        ax : Matplotlib axis instance
            Axis to plot.
            None will use the current axis.
        fig : Matplotlib figure instance
            Figure on which to add the plot.
            None will use the current figure.

        Notes
        -----
        Defaults are established during DYNAMO project analysis.
        """
        # parse parameters
        ax, fig = _parse_ax_fig(ax, fig)

        # Grab the variable dictionary of interest to plot
#        Var = self.fields[field]

        # Return masked or unmasked variable
#        Var, Data = self._get_variable_dict_data(field)
#        if mask_procedure is not None:
#            Data = get_masked_data(Data, mask_procedure, mask_tuple)

        # Return masked or unmasked variable
        # Subsetted if desired
        Var, tsub, Data = self._get_variable_dict_data_time_subset(
            field, start_time, end_time)
        if mask_procedure is not None:
            Data = get_masked_data(Data, mask_procedure, mask_tuple)

        # Find the closest vertical point
        Zind = find_nearest_indices(self.range['data'][:], range_height)
        print(f"Closest level:{self.range['data'][Zind]:4.1f}")

        # Convert lats/lons to 2D grid
        Ht2D, Lon2D = np.meshgrid(self.range['data'][:],
                                  self.longitude['data'][:])
        Ht2D, Lat2D = np.meshgrid(self.range['data'][:],
                                  self.longitude['data'][:])

        # Convert lats/lons to map projection coordinates

        x, y = self.basemap(Lon2D[:, Zind], Lat2D[:, Zind])

        # Create plots
        cs = self.basemap.pcolormesh(x, y, Data[:, Zind],
                                     vmin=vmin, vmax=vmax, cmap=cmap)

#    plt.colors.colormap.set_under('white')
        # Add Colorbar
        if color_bar:
            cbStr = "%s at %4.1f %s" % (Var['long_name'],
                                        self.height['data'][Zind],
                                        self.height['units'])
            cb = self.basemap.colorbar(
                cs, location=cb_loc, pad=cb_pad)
            cb.set_label(cbStr)
            # Set the number of ticks in the colorbar based upon
            # number of contours and the tick interval selected
            tick_locator = ticker.MaxNLocator(nbins=int(clevs / cb_tick_int))
            cb.locator = tick_locator
            cb.update_ticks()

        # Add title
        ax.set_title(title, fontsize=title_size)

        return

#     def overlay_wind_vector(self, height_level=2., vtrim=4, vlw=1.3, vhw=2.5,
#                             vscale=400, refVec=True, refU=10., refUposX=1.05,
#                             refUposY=1.015, qcolor='k', ax=None, fig=None):
#         """
#         Overlays a 2-D wind field at specified height onto map
#
#         Parameters
#         ----------
#         height_level : float
#             Height level for winds in horizontal plot.
#         vtrim : float
#             The number of vector arrows will be thinned by this factor.
#         vlw : float
#             Vector arrow linewidth.
#         vhw : float
#             Vector arrow headwidth.
#         vscale : int
#             Vector arrow scale (smaller = longer arrow).
#         refVec : boolean
#             True to attach a reference vector, False to not.
#         refU : float
#             Magnitude of reference vector [m/s].
#         refU_Xpos : float
#             X position of reference vector in axes-relative coordinates.
#         refU_Ypos : float
#             Y position of reference vector in axes-relative coordinates.
#         ax : Matplotlib axis instance
#             Axis to plot. None will use the current axis.
#         fig : Matplotlib figure instance
#             Figure to add the plot. None will use the current figure.
#
#         Notes
#         -----
#         U  is Along aircraft longitudinal axis wind.
#         V  is Perpendicular aircraft longitudinal axis wind.
#         """
#         # parse parameters
#         ax, fig = _parse_ax_fig(ax, fig)
#
#         # Find the closest vertical point for desired wind field
#         Htind = find_nearest_indices(self.height['data'][:], height_level)
#
#         # transform to porjection grid
#         U = self.fields['Uwind']['data'][Htind, :, :]
#         V = self.fields['Vwind']['data'][Htind, :, :]
#         lon = self.longitude['data'][:]
#         lat = self.latitude['data'][:]
#         uproj, vproj, xx, yy = self.basemap.transform_vector(
#             U, V, lon, lat, len(lon) / vtrim, len(lat) / vtrim,
#             returnxy=True, masked=True)
#
#         # Overplot the vectors
#         Q = self.basemap.quiver(xx, yy, uproj, vproj, scale=vscale,
#                                 headwidth=vhw, linewidths=vlw, color=qcolor)
#
#         # Make a quiver key to attach to figure.
#         qkLab = str("%4.1f" % refU) + 'm/s at ' + \
#             str("%4.1f" % height_level) + ' km'
#         # , fontproperties={'weight': 'bold'})
#         qk = ax.quiverkey(Q, refUposX, refUposY, refU, qkLab)
#
#         return
#
#     def plot_point(self, lon, lat, symbol='ro', label_text=None,
#                    label_offset=(None, None), **kwargs):
#         """
#         Plot a point on the current map.
#
#         Additional arguments are passed to basemap.plot.
#
#         Parameters
#         ----------
#         lon : float
#             Longitude of point to plot.
#         lat : float
#             Latitude of point to plot.
#         symbol : str
#             Matplotlib compatible string which specified the symbol of the
#             point.
#         label_text : str, optional.
#             Text to label symbol with.  If None no label will be added.
#         label_offset : [float, float]
#             Offset in lon, lat degrees bottom left corner of the label
#             text relative to the point. None will use 0.01 default.
#         """
#         lon_offset, lat_offset = label_offset
#         if lon_offset is None:
#             lon_offset = 0.01
#         if lat_offset is None:
#             lat_offset = 0.01
#
#         # Plot the symbol
#         self.basemap.plot(lon, lat, symbol, latlon=True, **kwargs)
#         # Attach the text
#         if label_text is not None:
#             # basemap does not have a text method so we must determine
#             # the x and y points and plot them on the basemap's axis.
#             x_text, y_text = self.basemap(lon + lon_offset,
#                                           lat + lat_offset)
#             self.basemap.ax.text(x_text, y_text, label_text)
#
#     def plot_line_geo(self, line_lons, line_lats,
#                       line_style='r-', lw=3, alpha=0.2,
#                       label0=True, label_offset=(0.01, 0.01),
#                       ax=None, fig=None, **kwargs):
#         """
#         Plot a line segments on the current map given values in lat and lon.
#
#         Additional arguments are passed to basemap.plot.
#
#         Parameters
#         ----------
#         line_lons : array
#             Longitude of line segment to plot.
#         line_lats : array
#             Latitude of line segment to plot.
#         line_style : str
#             Matplotlib compatible string which specifies the line style.
#         lw : int
#             LInewidth.
#         alpha : float
#             Transparency, 0: transparent, 1: opaque.
#         label0 : bool
#             True to label the starting point of line, False is no label.
#         label_offset : [float, float]
#             Offset in lon, lat degrees for bottom left corner of the label
#             text relative to the point. None will use 0.01 default.
#         ax : Matplotlib axis instance
#             Axis to plot. None will use the current axis.
#         fig : Matplotlib figure instance
#             Matplotlibl figure instance to add plot.
#             None will use the current figure.
#         """
#         # parse parameters
#         if ax is None:
#             ax = _parse_ax(ax)
#
#         x, y = self.basemap(line_lons, line_lats)
#         self.basemap.plot(x, y, line_style, lw=lw, alpha=alpha, ax=ax)
#
#         # Overplot the 0 point
#         if label0:
#             self.plot_point(line_lons[0], line_lats[0], 'ko',
#                             label_text='0', label_offset=label_offset,
#                             markersize=0.5)

########################
#   3-D plot methods  ##
########################

#     def grid_3d(self, surf_field, surf_min=-5., surf_max=5.,
#                    surf_cmap='RdBu_r', rstride=5, cstride=5,
#                    plot_contour=False, cont_field=None, range_height=2.,
#                    cminmax=(0., 60.), clevs=25, vmin=15., vmax=60.,
#                    cmap='gist_ncar', alpha=0., zlims=(-5., 5.), dlat=1.,
#                    dlon=1., plot_track=True, title=" ", fig=None, ax=None):
#         """
#         Wrapper to call the 3D plotting function for backwards compatability.
#         """
#         r3d = Radar3DPlot(self.radar, basemap=self.basemap)
#
#         r3d.DPJgrid_3d(self.airborne, surf_field, surf_min=surf_min,
#                        surf_max=surf_max, surf_cmap=surf_cmap,
#                        rstride=rstride, cstride=cstride,
#                        plot_contour=plot_contour, cont_field=cont_field,
#                        range_height=range_height,
#                        cminmax=cminmax, clevs=clevs, vmin=vmin, vmax=vmax,
#                        cmap=cmap, alpha=alpha,
#                        zlims=zlims, dlat=dlat, dlon=dlon,
#                        plot_track=plot_track,
#                        title=title, fig=fig, ax=ax)

#############################
#   Vertical plot methods  ##
#############################

    def plot_cross_section(self, field, start_pt, end_pt, xs_length=500,
                           mask_procedure=None, mask_tuple=None,
                           title=" ", title_size=20, cminmax=(0., 60.),
                           clevs=25, vmin=15., vmax=60., cmap='gist_ncar',
                           clabel='dBZ', color_bar=True, cb_pad=.05,
                           cb_orient='vertical', cb_tick_int=2,
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
        cminmax : tuple
            (min,max) values for controur levels.
        clevs : integer
            Number of contour levels.
        vmin : float
            Minimum contour value to display.
        vmax : float
            Maximum contour value to display.
        clabel : str
            Label for colorbar (e.g. units 'dBZ').
        title : str
            Plot title.
        title_size : int
            Font size of title to display.
        cmap : string
            Matplotlib color map to use
        color_bar : bool
            True to add colorbar, False does not.
        cb_pad : str
            Pad to move colorbar, in the form "5%", pos is to
            right for righthand location.
        cb_loc : str
            Location of colorbar, default is 'right', also available:
            'bottom', 'top', 'left'.
        cb_tick_int : int
            Interval to use for colorbar tick labels,
            higher number "thins" labels.
        ax : Matplotlib axis instance
            Axis to plot. None will use the current axis.
        fig : Matplotlib figure instance
            Figure on which to add the plot. None will use the current figure.
        '''
        rvp = RadarVerticalPlot(self.radar, basemap=self.basemap)

        rvp.plot_cross_section(
            field, start_pt, end_pt, xs_length=xs_length,
            mask_procedure=mask_procedure, mask_tuple=mask_tuple,
            title=title, title_size=title_size,
            cminmax=cminmax, clevs=clevs, vmin=vmin, vmax=vmax,
            cmap=cmap, clabel=clabel, color_bar=color_bar, cb_pad=cb_pad,
            cb_orient=cb_orient, cb_tick_int=cb_tick_int, ax=ax, fig=fig)

###################
#   Get methods  ##
###################

    def _get_variable_dict(self, field):
        '''Get the variable from the fields dictionary'''
        Var = self.fields[field]
        return Var

    def _get_variable_dict_data(self, field):
        '''Get the variable from the fields dictionary.'''
        Var, data = self.fields[field], self.fields[field]['data'][:]
        return Var, data

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

    def _get_variable_dict_data_time_subset(self, field, start_time, end_time):
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
