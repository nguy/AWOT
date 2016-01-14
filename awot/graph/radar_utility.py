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

from matplotlib.colors import from_levels_and_colors

from . import common

class RadarUtilityPlot(object):
    """Create vertical plot of radar data."""

    def __init__(self, radar, basemap=None,
                lon_name=None, lat_name=None, height_name=None,
                time_name=None):
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

        # See what the field shape looks like for inspection
        # of time and height arrays later
        fieldshape = self.radar['fields'][self.radar['fields'].keys()[0]]['data'].shape

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

        # Numpy meshgrid can do N-Dimensional. May be worth exploring,
        # BUT would need the indices variables.
        self.heightfield = self.height.copy()
        if len(self.heightfield['data'].shape) == 1:
            self.heightfield['data'] = np.resize(self.height['data'][:], fieldshape)

        # Attempt to pull in time if found
        if time_name is None:
            try:
                self.time = self.radar['time']
                self.timefield = self.time.copy()
                # This step is slow converting back-and-forth between datenum
                # instances, but Numpy does not support meshgrid with datenum dtype
                if len(self.timefield['data'].shape) == 1:
                    timenum = date2num(self.time['data'][:], self.time['units'])
                    self.timefield['data'] = num2date(
                          np.resize(timenum, fieldshape), self.time['units'])
            except:
                self.time = None
                self.timefield = None
        else:
            self.time = self.radar[time_name]

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
        ax = common._parse_ax(ax)
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

        bifreq, xedges, yedges = np.histogram2d(
                 xarr.ravel(), yarr.ravel(), bins=(binsx, binsy), normed=True)
        X, Y = np.meshgrid(xedges, yedges)

        cb_label = "Frequency"
        if plot_percent:
            bifreq = bifreq * 100.
            cb_label = cb_label + " (%)"

        if mask_below is not None:
            bifreq = np.ma.masked_where(bifreq < mask_below, bifreq)
        # Plot the data
        p = ax.pcolormesh(X, Y, bifreq.T, vmin=vmin, vmax=vmax, cmap=cmap)

        # Set the axes
        common._set_axes(ax, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max,
                  title=title, titleFontSize=titleFontSize,
                  xlab=xlab, ylab=ylab, xpad=xpad, ypad=ypad,
                  xlabFontSize=xlabFontSize, ylabFontSize=ylabFontSize)

        if plot_colorbar:
            cb = common.add_colorbar(ax, p, orientation=cb_orient, pad=cb_pad,
                 label=cb_label, fontsize=cb_fontsize,
                 ticklabel_size=cb_ticklabel_size,
                 clevs=cb_levs, tick_interval=cb_tick_int)

        # Create a dictionary to return
        bivar = {'frequency' : bifreq.T,
                 'frequency_label' : cb_label,
                 'xaxis' : X,
                 'yaxis' : Y}
        return bivar

    def plot_cfad(self, field, height_axis=1,
                  xbinsminmax=None, nbinsx=50,
                  points_thresh_fraction=0.1,
                  start_time=None, end_time=None,
                  vmin=None, vmax=None, cmap=None,
                  discrete_levels=True, levels=None,
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
                  quantiles=None,
                  qcolor='k', qlabels_on=True,
                  qlabel_color=None, qlabel_size=None,
                  ax=None, fig=None):
        """
        Create a frequency by altitude distribution plot of two variables.
        This is the traditional method of calculating a frequency distribution
        at each height of input array by iterating through the height array
        and input data array.

        NOTE: This routine only works with Cartesian data.

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
        points_thresh_fraction : float
            The fraction of points that must be present for the
            CFAD to be calculated. Following Yuter and Houzed 1995,
            the default values is 0.1 (10%) of potential data coverage
            is required. This threshold removes anomolous results when
            a small number of points is present.
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
        discrete_levels : boolean
            True to create discrete levels for colormap. False will apply
            the default linear normalization of luminance for colorbar scale.
        levels : array
            An list of levels to be used for display. If chosen discrete
            color will be used in the colorbar instead of a linear luminance
            mapping.
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
        ax = common._parse_ax(ax)
        # Snag the data from requested field
        xarr = self._get_fields_variable_dict_data_time_subset(
            field, start_time, end_time)

        if xbinsminmax is None:
            xbinsminmax = (np.ma.min(xarr), np.ma.max(xarr))
        binsx = np.linspace(xbinsminmax[0], xbinsminmax[1], nbinsx, endpoint=True)

        cb_label = "Frequency"
        percent = False

        if plot_percent:
            cb_label = cb_label + " (%)"

        # Reshape the array so that height axis is first dimension
        if height_axis != 0:
            ht = np.rollaxis(self.heightfield['data'], height_axis)
            xarr = np.rollaxis(xarr, height_axis)
        else:
            ht = self.heightfield['data'].copy()

        # Create CFAD array to fill
        nh = xarr.shape[0]
        bin_pts = np.empty((nh, len(binsx)-1))
        bin_perc = np.empty((nh, len(binsx)-1))

        for nn in range(nh):
#            CFAD[nn, :], bin_edges = np.histogram(
#                   xarr[nn, ...], bins=binsx, density=plot_percent)

            # Check data for good points
            condition = np.logical_or(
                         np.isfinite(xarr[nn, ...].ravel()),
                         xarr.mask[nn, ...].ravel() == 0)
            # Sort the good data from low to high values
            array = np.sort(xarr[nn, ...].ravel()[condition])
            # Calculate the fraction of points out of possible total
            ptsfrac = float(len(array))/float(len(xarr[nn, ...].ravel()))
            if ptsfrac > points_thresh_fraction:
                bin_pts[nn, :], bin_edges  = np.histogram(
                    array,  bins=binsx, density=False)
            bin_perc[nn, :] = bin_pts[nn, :] / bin_pts[nn, :].sum() * 100.
            del(ptsfrac, array, condition)

        # Mask any invalid or negative numbers
        bin_pts = np.ma.masked_invalid(bin_pts)
        bin_pts = np.ma.masked_less(bin_pts, 0.)
        bin_perc = np.ma.masked_invalid(bin_perc)
        bin_perc = np.ma.masked_less(bin_perc, 0.)

        if plot_percent:
            CFAD = bin_perc
        else:
            CFAD = bin_pts

        X, Y = np.meshgrid(binsx, self.height['data'][:])
        if mask_below is not None:
            CFAD = np.ma.masked_where(CFAD < mask_below, CFAD)

        # Set the axes
        common._set_axes(ax, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max,
                  title=title, titleFontSize=titleFontSize,
                  xlab=xlab, ylab=ylab, xpad=xpad, ypad=ypad,
                  xlabFontSize=xlabFontSize, ylabFontSize=ylabFontSize)

        # Plot the data
        norm, levpos, colors = None, None, None
        cm = plt.get_cmap(cmap)
        if discrete_levels:
            # Default to these levels if none chosen
            if levels is None:
                levels = [.1, .5, 1, 2, 5, 7, 10, 15, 25]
            # Get the colormap and calculate data spaced by number of levels

            levpos = np.rint(np.squeeze([np.linspace(0, 255, len(levels))])).astype(int)
            colors = cm(levpos)
            # Convert levels to colormap values
            cmap, norm = from_levels_and_colors(levels, colors, extend='max')
        p = ax.pcolormesh(X, Y, CFAD, vmin=vmin, vmax=vmax, norm=norm, cmap=cmap)

        if plot_colorbar:
            cb = common.add_colorbar(ax, p, orientation=cb_orient, pad=cb_pad,
                 label=cb_label, fontsize=cb_fontsize,
                 ticklabel_size=cb_ticklabel_size,
                 clevs=cb_levs, tick_interval=cb_tick_int)

        if quantiles is not None:
            qArr = self.plot_quantiles(field, quantiles=quantiles,
                                       height_axis=height_axis,
                                       start_time=start_time, end_time=end_time,
                                       qcolor=qcolor, qlabels_on=qlabels_on,
                                       qlabel_color=qlabel_color,
                                       qlabel_size=qlabel_size,
                                       setup_axes=False, ax=ax)


        # Clean up any potentially lingering variables
        del(norm, levels, cm, levpos, colors, CFAD)
        # Populate the dictionary
        cfad_dict = {'frequency_points' : bin_pts,
                     'frequency_percent' : bin_perc,
                     'xaxis' : X,
                     'yaxis' : Y
                     }
        return cfad_dict

    def plot_quantiles(self, field, quantiles=None, height_axis=1,
                  start_time=None, end_time=None,
                  qcolor='k', qlabels_on=True,
                  qlabel_color=None, qlabel_size=None,
                  mask_above_height=None, mask_below_height=None,
                  mask_between_height=None,
                  x_min=None, x_max=None,
                  y_min=None, y_max=None,
                  xlab=None, xlabFontSize=None, xpad=None,
                  ylab=None, ylabFontSize=None, ypad=None,
                  title=None, titleFontSize=None,
                  setup_axes=True,
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
        quantiles : list
            A list of percentage values for quantile calculations.
        height_axis : int
            The dimension to perform quantile calculation over (non-height).
        start_time : str
            UTC time to use as start time for subsetting in datetime format.
            (e.g. 2014-08-20 12:30:00)
        end_time : str
            UTC time to use as an end time for subsetting in datetime format.
            (e.g. 2014-08-20 16:30:00)
        qcolor : str
            Color to use for quantile plot lines.
        qlabels_on : boolean
            True to print labels of quantiles on plot. False for no labels.
        qlabel_color : str
            Color of quantile labels if activated. Default black.
        qlabel_size : int
            Size of the quantile labels.
        x_min : float
            Minimum value for X-axis.
        x_max : float
            Maximum value for X-axis.
        y_min : float
            Minimum value for Y-axis.
        y_max : float
            Maximum value for Y-axis.
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
        setup_axes : bool
            If True, the axes will be set up using either keyword information
            or attempted to determine automatically.
            If False, existing axes will be used, either passed by keyword ax
            or current working axis.
        ax : Matplotlib axis instance
            Axis to plot. None will use the current axis.
        fig : Matplotlib figure instance
            Figure on which to add the plot. None will use the current figure.
        """
        # parse parameters
        ax = common._parse_ax(ax)
        # Snag the data from requested field
        xarr = self._get_fields_variable_dict_data_time_subset(
            field, start_time, end_time)

        if quantiles is None:
            quantiles = [5, 10, 25, 50, 75, 90]

        # Reshape the array so that height axis is first dimension
        if height_axis != 0:
            ht = np.rollaxis(self.heightfield['data'], height_axis)
            xarr = np.rollaxis(xarr, height_axis)
        else:
            ht = self.heightfield['data'].copy()

        # Apply mask to altitudes if indicated
        ### NG 	THESE ARE NOT WORKING PROPERLY DON'T KNOW WHY ###
        if mask_above_height is not None:
            xarr = np.ma.masked_where(ht > mask_above_height, xarr)
        if mask_below_height is not None:
            xarr = np.ma.masked_where(ht < mask_below_height, xarr)
        if mask_between_height is not None and np.shape(mask_between_height) != 2:
            condc = ((ht < mask_between_height[0]) &
                     (ht > mask_between_height[1]))
            xarr = np.ma.masked_where(condc, xarr)

         # Check data for good points
        xarr = np.ma.masked_where(~(np.isfinite(xarr)), xarr)

        # Calculate the quantile profiles
        qArr = self.calc_quantiles(xarr, self.height['data'][:], quantiles)

        # Set the axes
        if setup_axes:
            common._set_axes(ax, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max,
                      title=title, titleFontSize=titleFontSize,
                      xlab=xlab, ylab=ylab, xpad=xpad, ypad=ypad,
                      xlabFontSize=xlabFontSize, ylabFontSize=ylabFontSize)

        # Plot the data
        self.add_quantiles_to_axis(ax, qArr, qcolor, qlabels_on)
        return qArr

    def add_quantiles_to_axis(self, ax, qArr, qcolor, qlabels_on, qlabel_size=10,
                              qlabel_color='k'):
        for num in range(len(qArr['quantiles'])):
            p = ax.plot(qArr['profiles'][:, num], qArr['yaxis'][:], color=qcolor)
            ytextloc = self.heightfield['data'][:, 0].max() * 1.05
            if qlabels_on:
                ax.text(qArr['profiles'][-1, num], ytextloc,
                        str(qArr['quantiles'][num]), color=qlabel_color,
                        size=qlabel_size)

###########################
#   Calculation methods   #
###########################
    def calc_cfad_quantile(self, data, binsx, points_thresh_fraction, quantiles=None):
        # Create array to fill
        nh = data.shape[0]
        bin_pts = np.empty((nh, len(binsx)-1))
        bin_perc = np.empty((nh, len(binsx)-1))

        if quantiles is not None:
            pass

        for nn in range(nh):
            # Check data for good points
            condition = np.logical_or(
                         np.isfinite(data[nn, ...].ravel()),
                         data.mask[nn, ...].ravel() == 0)
            # Sort the good data from low to high values
            data = np.sort(data[nn, ...].ravel()[condition])
            # Calculate the fraction of points out of possible total
            ptsfrac = float(len(data))/float(len(data[nn, ...].ravel()))
            if ptsfrac > points_thresh_fraction:
                bin_pts[nn, :], bin_edges  = np.histogram(
                    data,  bins=binsx, density=False)
            bin_perc[nn, :] = bin_pts[nn, :] / bin_pts[nn, :].sum() * 100.
            del(ptsfrac, data, condition)

        # Mask any invalid or negative numbers
        bin_pts = np.ma.masked_invalid(bin_pts)
        bin_pts = np.ma.masked_less(bin_pts, 0.)
        bin_perc = np.ma.masked_invalid(bin_perc)
        bin_perc = np.ma.masked_less(bin_perc, 0.)
        cfad_dict = {'frequency_points' : bin_pts,
                     'frequency_percent' : bin_perc,
                     'xaxis' : X,
                     'yaxis' : Y
                     }
        return cfad_dict

    def calc_quantiles(self, xarr, height, quantiles):
        # Create array to fill
        nh = xarr.shape[0]
        qArr = np.empty((nh, len(quantiles)))
        for nn in range(nh):
#             # Sort the good data from low to high values
#             condition = (xarr.mask[nn, ...] == True)
#             data = np.sort(xarr[nn, ...][condition].ravel())
            data = np.sort(xarr[nn, ...].ravel())
#             inds = np.rint(np.array(quantiles)/100. * len(data)).astype(int)
#             qArr[nn, :] = data[inds]
            qArr[nn, :] = mstats.mquantiles(data, prob=[x / 100. for x in quantiles], axis=None)
        quant_dict = {'quantiles' : quantiles,
                      'profiles' : qArr,
                      'yaxis' : height
                      }
        return quant_dict

###################
#   Get methods   #
###################

    def _get_fields_variable_dict_data_time_subset(self, field, start_time, end_time):
        '''
        Get the variable from the fields dictionary.
        Subset the time when in time series format.
        '''
        data = self.fields[field]['data'][:]

        # Check to see if time is subsetted
        if self.time != None:
            dt_start = common._get_start_datetime(self.time, start_time)
            dt_end = common._get_end_datetime(self.time, end_time)
            tsub = self.time['data'][(self.time['data'] >= dt_start) &
                                     (self.time['data'] <= dt_end)]
            datasub = data[(self.time['data'] >= dt_start) &
                           (self.time['data'] <= dt_end)]
        else:
            datasub = data
        return datasub

    def _get_variable_dict_data_time_subset(self, field, start_time, end_time):
        '''
        Get the variable from the fields dictionary.
        Subset the time when in time series format.
        '''
        data = self.radar[field]['data'][:]

        # Check to see if time is subsetted
        if self.time != None:
            dt_start = common._get_start_datetime(self.time, start_time)
            dt_end = common._get_end_datetime(self.time, end_time)
            tsub = self.time['data'][(self.time['data'] >= dt_start) &
                                     (self.time['data'] <= dt_end)]
            datasub = data[(self.time['data'] >= dt_start) &
                           (self.time['data'] <= dt_end)]
        else:
            datasub = data
        return datasub

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
        if self.time != None:
            dt_start = common._get_start_datetime(self.time, start_time)
            dt_end = common._get_end_datetime(self.time, end_time)

            data1sub = data1[(self.timefield['data'] >= dt_start) &
                           (self.timefield['data'] <= dt_end)]

            data2sub = data2[(self.timefield['data'] >= dt_start) &
                           (self.timefield['data'] <= dt_end)]
        else:
            data1sub, data2sub = data1, data2
        return data1sub, data2sub

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
