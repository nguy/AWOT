"""
awot.graph.flight_level
=========================

A group of scripts to create plots of flight level data. 

Created by Nick Guy.

TODO:
Fix the time_stamp module to work with savefig.  Okay on screen output,
but for some reason it bombs when trying to save.

Improve time series variable plotting?
"""
# HISTORY::
# 20 Aug 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
#-------------------------------------------------------------------
# Load the needed packages
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from datetime import datetime

from .common import plot_date_ts
#===============================================================
# BEGIN FUNCTIONS
#===============================================================
class FlightLevel(object):
    """
    To create a FlightLevel instance:
    new_instance = FlightLevel() or new_instance = FlightLevel(AirborneInstance)
    
    Notable attributes
    ------------------
    
    A basemap call can be passed if created externally.
    """
    def __init__(self, airborne, basemap=None):
        '''Intitialize the class to create plots'''
#        self.basemap = airborne.basemap
        self.longitude = airborne.flight_data['longitude']
        self.latitude = airborne.flight_data['latitude']
        self.altitude = airborne.flight_data['altitude']
        self.flight_number = airborne.flight_data['flight_number']
        self.project = airborne.flight_data['project']
        self.platform = airborne.flight_data['platform']
        self.Uwind = airborne.flight_data['Uwind']
        self.Vwind = airborne.flight_data['Vwind']
        self.time = airborne.flight_data['time']
        self.flight_data = airborne.flight_data
        
        if basemap != None:
            self.basemap = basemap
        try:
            self.basemap = airborne.basemap
        except:
            print "Warning: No basemap instance"
            
        # Calculate x,y map position coordinates
        self.x, self.y = self.basemap(self.longitude, self.latitude)        
       
    ###############
    # Track plots #
    ###############

    def plot_trackmap(self,
                   color_by_altitude=False, track_cmap='spectral', 
                   track_color='k', track_lw=1.5, alpha =1.,
                   start_time=None, end_time=None,
                   
                   min_altitude=None, max_altitude=None,
                   add_cb=True, cbloc='right', cbsize='3%', cbpad='10%',

                   addlegend=False, legLoc='lower right', legLab=None,
                   addtitle=False, title=None,
                   ax=None, fig=None, **kwargs):
        """
        Create a plot of the aircraft flight track, with optional enhancements.
        Note that the legend only works with single color track at this time.
    
        Parameters::
        ----------
        color_by_altitude : bool
            True results in flight track changing color with altitude
            False displays track in a single color
        track_cmap : colormap
            If color_by_altitude True, it will use this colormap
        track_color : string
            If color_by_altitude False, this color (see matplotlib) is used for track
        track_lw : float or int
            Line width to use for track
        alpha : float
            Alpha value for transparency (0 = transparent, 1. = solid)
        start_time : string
            UTC time to use as start time for subsetting in datetime format
            (e.g. 2014-08-20 12:30:00)
        end_time : string
            UTC time to use as an end time for subsetting in datetime format
            (e.g. 2014-08-20 16:30:00)
                
        min_altitude : float
            Minimum altitude to use in mapping color
        max_altitude : float
            Maximum altitude to use in mapping color
        add_cb : boolean
            True plots a colorbar, Fasle does not
        cbloc : 'str'
            Location of colorbar.  
        cbsize : 'str'
            Size of colorbar 'N%'
        cbpad : 'str'
            Pad between parent axes and colorbar axes in same units as size. 'N%'
            
        addlegend : bool
            Defaults to No legend drawn
        legLoc : str
            See matplotlib legend object documentation
        addtitle : bool
            Defaults to no title for plot
        title : str
            See matplotlib axis object documentation
                
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.
                
        **kwargs will pass specific arguments to subprograms, 
          see Basemap and Matplotlib.
        """
#        plt.close()
        
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
        
        # Check inputs
        if legLab is None:
            legLab = self.flight_number
            
        if min_altitude is None:
            min_altitude = self.altitude.min()
        
        if max_altitude is None:
            max_altitude = self.altitude.max()
            
        # Get start and end times (this deals with subsets)
        dt_start = self._get_start_datetime(start_time)
        dt_end = self._get_end_datetime(end_time)
        
        # Subset the data (will use min/max if None given)
        xSub, ySub = self._get_x_y_time_subset(start_time, end_time, return_time=False)
        lonSub, latSub = self._get_lon_lat_time_subset(start_time, end_time)
        timeSub, VarSub = self._get_time_var_time_subset('altitude', start_time, end_time)
        
        # Clean up the masked data for plotting
        all_mask = latSub.mask + lonSub.mask + VarSub.mask
        lonmask = np.ma.array(lonSub, mask=all_mask).compressed()
        latmask = np.ma.array(latSub, mask=all_mask).compressed()
        varmask = np.ma.array(VarSub, mask=all_mask).compressed()
        
        xmask, ymask = self.basemap(lonmask, latmask)

        # Plot the track either coloring by altitude or as single color  
        if color_by_altitude:
            lc = self._colorline(xmask, ymask, z=varmask,
                            cmap=track_cmap, linewidth=track_lw, alpha=alpha,
                            vmin=min_altitude, vmax=max_altitude)
    
            ax.add_collection(lc)
                            
            if add_cb:
                cb = self.basemap.colorbar(lc, location=cbloc, size=cbsize, pad=cbpad)
                cb.set_label('Altitude (m)')
    
        else:
            p = ax.plot(xmask, ymask,
                   color=track_color, lw=track_lw, alpha=alpha, label=legLab)
                   
        if addlegend:
            self.basemap.ax.legend(loc=legLoc, fontsize=11, frameon=True,  
                       handlelength=2, handletextpad=1, labelspacing=.2)
                       
        if addtitle:
            if title is None:
                title = self.project +' '+self.platform
            self.basemap.ax.set_title(title)
            
    ####### 
       
    def plot_trackmap_variable(self, 
                   field=None, track_lw=1.5, alpha=1.,
                   start_time=None, end_time=None,
                   
                   track_cmap='jet',
                   min_value=None, max_value=None,
                   add_cb=True, cbloc='right', cbsize='3%', cbpad='10%', 
                   cblabel='Temperature (C)',
                   
                   addlegend=False, legLoc='lower right', legLab=None,
                   addtitle=False, title=None,
                   ax=None, fig=None, **kwargs):
        """
        Create a plot of the aircraft flight track, with optional enhancements.
        Note that the legend only works with single color track at this time.
    
        Parameters::
        ----------
        field : float
            Variable passed to plot, defaults to temperature
        track_lw : float or int
            Line width to use for track
        alpha : float
            Alpha value for transparency (0 = transparent, 1. = solid)
        start_time : string
            UTC time to use as start time for subsetting in datetime format
            (e.g. 2014-08-20 12:30:00)
        end_time : string
            UTC time to use as an end time for subsetting in datetime format
            (e.g. 2014-08-20 16:30:00)
            
        track_cmap : colormap
            Colormap to use for coding values
        min_value : float
            Minimum value to use in mapping color
        max_value : float
            Maximum value to use in mapping color
        add_cb : boolean
            True plots a colorbar, Fasle does not
        cbloc : 'str'
            Location of colorbar.  
        cbsize : 'str'
            Size of colorbar 'N%'
        cbpad : 'str'
            Pad between parent axes and colorbar axes in same units as size. 'N%'
            
        addlegend : bool
            Defaults to No legend drawn
        legLoc : str
            See matplotlib legend object documentation
        addtitle : bool
            Defaults to no title for plot
        title : str
            See matplotlib axis object documentation
            
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.
                
        **kwargs will pass specific arguments to subprograms, 
          see Basemap and Matplotlib.
        OUTPUT::
            p : plot instance
                Plot
        """
        
        # Pull out the variable to plot
        if field is None:
            field = 'temperature'
        
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
            
        # Get start and end times (this deals with subsets)
        dt_start = self._get_start_datetime(start_time)
        dt_end = self._get_end_datetime(end_time)
        
        # Subset the data (will use min/max if None given)
        xSub, ySub = self._get_x_y_time_subset(start_time, end_time, return_time=False)
        lonSub, latSub = self._get_lon_lat_time_subset(start_time, end_time)
        timeSub, VarSub = self._get_time_var_time_subset(field, start_time, end_time)
        
        # Clean up the masked data for plotting
        all_mask = latSub.mask + lonSub.mask + VarSub.mask
        lonmask = np.ma.array(lonSub, mask=all_mask).compressed()
        latmask = np.ma.array(latSub, mask=all_mask).compressed()
        varmask = np.ma.array(VarSub, mask=all_mask).compressed()
        
        xmask, ymask = self.basemap(lonmask, latmask)
        
        # Check inputs
        if min_value is None:
            min_value = varmask.min()
         
        if max_value is None:
            max_value = varmask.max()
        
        if legLab is None:
            legLab = self.flight_number
           
        # Plot the track
        lc = self._colorline(xmask, ymask, z=varmask, cmap=track_cmap,
#        lc = self._colorline(z=varm, cmap=track_cmap, 
                        linewidth=track_lw, alpha=alpha,
                        vmin=min_value, vmax=max_value)
    
        ax.add_collection(lc)
        
        if add_cb:                    
            cb = self.basemap.colorbar(lc, location=cbloc, size=cbsize, pad=cbpad)
            cb.set_label(cblabel)
                   
        if addlegend:
            self.basemap.ax.legend(loc=legLoc, fontsize=11, frameon=True,  
                       handlelength=2, handletextpad=1, labelspacing=.2)
                       
        if addtitle:
            if title is None:
                title = self.project +' '+self.platform
            ax.set_title(title)
            
    ##########################
    # Basemap add-on methods #
    ##########################
    
    def draw_boundary(self, **kwargs):
        """
        Draw a boundary around the map
        Parameters::
        ----------
           See basemap documentation
        """
        
        self.basemap.drawmapboundary(**kwargs)
    
    def draw_scale(self, location='lower left',
                   lon=None, lat=None, lon0=None, lat0=None, length=None, **kwargs):
        """
        Draw a reference scale on the map
        Parameters::
        ----------
        location : str
           Where to put the scale, default is 'lower left'
           Options are 'lower left', 'upper left', 'lower right', 'upper right',
            'upper middle', 'lower middle'
        
        All other parameters, See basemap documentation
        
        Note that the location parameter can be overridden by setting 
        the lon and lat variables explicitly
        """
        
        # parse parameters
#        ax = self._parse_ax(ax)
        
        # Check on where to put the label, defaults to the else statement
        # that calculates 'lower left'
        if (location.lower() == 'upper left') or (location.lower() == 'upper_left'):
            xScaleOffset = (self.basemap.urcrnrlon - self.basemap.llcrnrlon) * 0.10
            yScaleOffset = (self.basemap.urcrnrlat - self.basemap.llcrnrlat) * 0.90
        elif (location.lower() == 'upper middle') or (location.lower() == 'upper_middle'):
            xScaleOffset = (self.basemap.urcrnrlon - self.basemap.llcrnrlon) * 0.50
            yScaleOffset = (self.basemap.urcrnrlat - self.basemap.llcrnrlat) * 0.90
        elif (location.lower() == 'upper right') or (location.lower() == 'upper_right'):
            xScaleOffset = (self.basemap.urcrnrlon - self.basemap.llcrnrlon) * 0.90
            yScaleOffset = (self.basemap.urcrnrlat - self.basemap.llcrnrlat) * 0.90
        elif (location.lower() == 'lower right') or (location.lower() == 'lower_right'):
            xScaleOffset = (self.basemap.urcrnrlon - self.basemap.llcrnrlon) * 0.90
            yScaleOffset = (self.basemap.urcrnrlat - self.basemap.llcrnrlat) * 0.10
        elif (location.lower() == 'lower middle') or (location.lower() == 'lower_middle'):
            xScaleOffset = (self.basemap.urcrnrlon - self.basemap.llcrnrlon) * 0.50
            yScaleOffset = (self.basemap.urcrnrlat - self.basemap.llcrnrlat) * 0.10
        else :
            xScaleOffset = (self.basemap.urcrnrlon - self.basemap.llcrnrlon) * 0.10
            yScaleOffset = (self.basemap.urcrnrlat - self.basemap.llcrnrlat) * 0.10
                        
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
                                  length=length, barstyle='fancy', labelstyle='simple',
                                  **kwargs)
    
    def draw_barbs(self, barbspacing=10, barbcolor='k', flagcolor='k', 
                   blength=8, lw=0.5, **kwargs):
        """
        Draw a reference scale on the map
        Parameters::
        ----------
        barbspacing : int
           Plot every nth via the barbspacing provided.
        See basemap documentation
        """
                           
        # Only plot every nth barb from the barbspacing parameter
        self.basemap.barbs(self.longitude[::barbspacing], self.latitude[::barbspacing],
                           self.Uwind[::barbspacing], self.Vwind[::barbspacing],  
                           latlon=True, barbcolor=barbcolor, flagcolor=flagcolor, 
                           linewidth=lw, **kwargs)
    
    def plot_point(self, lon, lat, symbol='ro', label_text=None,
                   label_offset=(None, None), **kwargs):
        """
        Plot a point on the current map.

        Additional arguments are passed to basemap.plot.

        Parameters::
        ----------
        lon : float
            Longitude of point to plot.
        lat : float
            Latitude of point to plot.
        symbol : str
            Matplotlib compatible string which specified the symbol of the
            point.
        label_text : str, optional.
            Text to label symbol with.  If None no label will be added.
        label_offset : [float, float]
            Offset in lon, lat degrees for the bottom left corner of the label
            text relative to the point. A value of None will use 0.01 default.
        """
        
        lon_offset, lat_offset = label_offset
        if lon_offset is None:
            lon_offset = 0.01
        if lat_offset is None:
            lat_offset = 0.01
            
        # Plot the symbol
        self.basemap.plot(lon, lat, symbol, latlon=True, **kwargs)
        # Attach the text
        if label_text is not None:
            # basemap does not have a text method so we must determine
            # the x and y points and plot them on the basemap's axis.
            x_text, y_text = self.basemap(lon + lon_offset,
                                          lat + lat_offset)
            self.basemap.ax.text(x_text, y_text, label_text)
            
    def plot_line_geo(self, line_lons, line_lats, line_style='r-', **kwargs):
        """
        Plot a line segments on the current map given values in lat and lon.

        Additional arguments are passed to basemap.plot.

        Parameters::
        ----------
        line_lons : array
            Longitude of line segment to plot.
        line_lats : array
            Latitude of line segment to plot.
        line_style : str
            Matplotlib compatible string which specifies the line style.
        """
        
        self.basemap.plot(line_lons, line_lats, line_style, latlon=True,
                          **kwargs)
  
    def time_stamps(self, labelspacing=1800, symbol='k*',
                    size=12, color='k',
                    label_offset=(None, None), 
                   start_time=None, end_time=None,
                    ax=None, **kwargs):
        """
        Add time stamps along track according to labelspacing.
        
        Parameters::
        ----------
        labelspacing : int
            Plot every nth via the labelpacing provided. 
            The default assumes seconds and labels every 30 minutes.
        symbol : str
            Color-symbol combo per matplotlib.pyplot.plot
        size : int
            Font size used
        color : str
            Matplotlib color
        label_offset : [float, float]
            Offset in lon, lat degrees for the bottom left corner of the label
            text relative to the point. A value of None will use 0.01 default.
        start_time : string
            UTC time to use as start time for subsetting in datetime format
            (e.g. 2014-08-20 12:30:00)
        end_time : string
            UTC time to use as an end time for subsetting in datetime format
            (e.g. 2014-08-20 16:30:00)
        ax : Axis
            Optional Axes instance
        """
        
        # parse parameters
        ax = self._parse_ax(ax)
                           
        # Only plot every nth barb from the barbspacing parameter
        lon_offset, lat_offset = label_offset
        if lon_offset is None:
            lon_offset = 0.01
        if lat_offset is None:
            lat_offset = 0.01
            
        # Get start and end times (this deals with subsets)
        dt_start = self._get_start_datetime(start_time)
        dt_end = self._get_end_datetime(end_time)
   
        time = self.time[(self.time >= dt_start) & (self.time <= dt_end)]
        lon = self.longitude[(self.time >= dt_start) & (self.time <= dt_end)]
        lat = self.latitude[(self.time >= dt_start) & (self.time <= dt_end)]
            
        # Find the number of points to plot
        labeltimes = time[::labelspacing]
        lons = lon[::labelspacing]
        lats = lat[::labelspacing]
        xpos_star, ypos_star = self.basemap(lons, lats)
        xpos_text, ypos_text = self.basemap(lons + lon_offset, lats + lat_offset)
        
        for nn in range(len(labeltimes)):
            # Plot the symbol
            ax.plot(xpos_star[nn], ypos_star[nn], symbol, **kwargs)
            label_text = labeltimes[nn].strftime("%H:%M")

            # Attach the text
            ax.text(xpos_text[nn], ypos_text[nn], label_text, fontsize=size, color=color)
            
    #######################
    # Time Series methods #
    #######################
    
    def plot_timeseries(self, field, colF='ko', msize=1.5, lw=2,
                dForm='%H:%M',tz=None, xdate=True,
                date_MinTicker='minute',
                other_MajTicks=None, other_MinTicks=None,
                
                start_time=None, end_time=None,
                
                title=None,
                xlab=' ', xlabFontSize=16, xpad=7,
                ylab=' ', ylabFontSize=16, ypad=7,
                ax=None):
        """
        Wrapper function to produce a time series plot of variable indicated
        Parameters::
        ----------
        field : float
            Variable to plot as time series
        colF : str
            Color and marker shortcut (see python documentation)
        msize : float
            Marker size
        lw : float
            Linewidth to use with line
        dForm : str
            Format of the time string for x-axis labels
        tz : str
            Time zone info to use when creating axis labels (see datetime)
        xdate : boolean
            True to use X-axis as date axis, false implies Y-axis is date axis
        date_MinTicker : str
            Sting to set minor ticks of date axis,
            'second','minute','hour','day' supported
        other_MajTicks : float
            Values for major tickmark spacing, non-date axis
        other_MinTicks : float
            Values for minor tickmark spacing, non-date axis
    
        start_time : string
            UTC time to use as start time for subsetting in datetime format
            (e.g. 2014-08-20 12:30:00)
        end_time : string
            UTC time to use as an end time for subsetting in datetime format
            (e.g. 2014-08-20 16:30:00)
    
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
        ax : Axes instance
            Optional axes instance to plot the graph
        """
        # parse parameters
        ax = self._parse_ax(ax)
        # Get the subsetted data
        tSub, VarSub = self._get_time_var_time_subset(field, 
                                             start_time=start_time, end_time=end_time)
        
        # Plot the time series
        ts = plot_date_ts(tSub, VarSub, 
                colF=colF, msize=msize, lw=lw,
                dForm=dForm,tz=tz, xdate=xdate, 
                date_MinTicker=date_MinTicker,
                other_MajTicks=other_MajTicks, other_MinTicks=other_MinTicks,
                title=title,
                xlab=xlab, xlabFontSize=xlabFontSize, xpad=xpad,
                ylab=ylab, ylabFontSize=ylabFontSize, ypad=ypad,
                ax=ax)
                            
    ######################
    # Line methods  #
    ######################
    # Original Multiline color plotting found at 
    # http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    # Modified for this class# 
    
    def _colorline(self, x, y, z=None, cmap=None, norm=None, linewidth=1.5, alpha=0.5,
                   vmin=None, vmax=None):
        '''
        Plot a colored line with coordinates x and y
        Optionally specify colors in the array z
        Optionally specify a colormap, a norm function and a line width

        Parameters::
        ----------
        x : float array
            X coordinates to map
        Y : float array
            Y coordinates to map
        z : float array
            Variable data to plot
        cmap : str
            String of colormap to use
        norm : float array
            Normalized array for colormapping
        linewidth : float
            Sets the width of the line drawn
        alpha : float
            Sets the transparency of the line (0 = transparent, 1 = opaque)
        vmin : float
            Sets the minimum to scale colormap luminance
        vmax : float
            Sets the maximum to scale colormap luminance
        '''
                          
        # Get the colormap to work with
        if cmap is None:
            cmap = 'spectral'
        
        # Set the data to be plotted:
        if z is None:
            data = self.altitude.copy()
        else:
            data = z.copy()
            
        if vmin is None:
            vmin = data.min()
        
        if vmax is None:
            vmax = data.max()
        
        data = np.ma.masked_outside(data, vmin, vmax)
           
        # Set the normalization to the min and max
        if norm is None:
            norm = plt.Normalize(vmin, vmax)

        # Create list of line segments from x and y coordinates, 
        # in the correct format for LineCollection:
        # an array of the form   numlines x (points per line) x 2 (x and y) array
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
        lc = LineCollection(segments, cmap=cmap, norm=norm, 
                            linewidth=linewidth, alpha=alpha)
        lc.set_array(data)

        return lc
        
    # NOTE: This is in here for legacy, RECOMMENDED NOT TO USE 
    def draw_marker_track(self,cmap=None,msize=2,mstyle='.'):
        """
        Draws a track by plotting each marker and colored based on height
        This is REALLY slow for large number of points.  
        Left in for reference.
        
        Highly recommended to not use this function
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

        # Loop through each point.  Unfortunately as of Aug 19, 2014 basemap.scatter
        # was not working.  No display and I have no idea why!  This is SLOW and 
        # not recommended.
        for lon, lat, alt in zip(self.longitude.compressed(), 
                             self.latitude.compressed(), self.altitude.compressed()):
            mcol = get_marker_color(alt)
            self.basemap.plot(lon, lat, color=mcol, latlon=True,
                            marker=mstyle, markersize=msize)             

    ####################
    # Get methods #
    ####################     
                            
    def _get_start_datetime(self, start_time):
        '''
        Get a start time as datetime instance for subsetting.
        '''
        # Check to see if time is subsetted
        if start_time is None:
            dt_start = self.time.min()
        else:
            startStr = [start_time[0:4],start_time[5:7],start_time[8:10],
                        start_time[11:13],start_time[14:16],start_time[17:19],'0']
            startInt = [ int(s) for s in startStr ]
            try:
                dt_start = datetime(startInt[0],startInt[1],startInt[2],startInt[3],
                                    startInt[4],startInt[5],startInt[6])
            except:
                print "Check the format of date string (e.g. '2014-08-20 12:30:00')"
                return
                
        return dt_start    
                            
    def _get_end_datetime(self, end_time):
        '''
        Get a start time as datetime instance for subsetting.
        '''
        #Check to see if the time is subsetted
        if end_time is None:
            dt_end = self.time.max()
        else:
            endStr = [end_time[0:4],end_time[5:7],end_time[8:10],
                        end_time[11:13],end_time[14:16],end_time[17:19],'0']
            endInt = [ int(s) for s in endStr ]           
            try:
                dt_end = datetime(endInt[0],endInt[1],endInt[2],endInt[3],
                                    endInt[4],endInt[5],endInt[6])
            except:
                print "Check the format of date string (e.g. '2014-08-20 12:30:00')"
                return
                    
        return dt_end
                      
    def _get_x_y_time_subset(self, start_time, end_time, return_time=False):
        '''
        Get a subsetted X and Y to control track length if input by user.
        '''
        # Check to see if time is subsetted
        dt_start = self._get_start_datetime(start_time)
        dt_end = self._get_end_datetime(end_time)
        
        x = self.x[(self.time >= dt_start) & (self.time <= dt_end)]
        y = self.y[(self.time >= dt_start) & (self.time <= dt_end)]
        
        if return_time:
            time = self.time[(self.time >= dt_start) & (self.time <= dt_end)]
            
            return x, y, time
        else:
        
            return x, y
                      
    def _get_lon_lat_time_subset(self, start_time, end_time):
        '''
        Get a subsetted Lon and Lat to control track length if input by user.
        '''
        # Check to see if time is subsetted
        dt_start = self._get_start_datetime(start_time)
        dt_end = self._get_end_datetime(end_time)
        
        lon = self.longitude[(self.time >= dt_start) & (self.time <= dt_end)]
        lat = self.latitude[(self.time >= dt_start) & (self.time <= dt_end)]
        
        return lon, lat
        
    def _get_time_var_time_subset(self, field, start_time, end_time):
        '''
        Get a subsetted time and Variable to control track length if input by user.
        '''
        # Check to see if time is subsetted
        dt_start = self._get_start_datetime(start_time)
        dt_end = self._get_end_datetime(end_time)
        
        tsub = self.time[(self.time >= dt_start) & (self.time <= dt_end)]
        var = self.flight_data[field]
        vsub = var[(self.time >= dt_start) & (self.time <= dt_end)]
        
        return tsub, vsub
        
        
    ####################
    # Save methods #
    ####################     
                           
    def save_figure(self, figName='awot_plot', figType="png"):
        '''Save the current plot
        
        Parameters::
        ------------
        figName : str
            Figure name
        figType : str
            Figure format, default to .png
        
        '''
        plt.savefig(figName+'.'+figType, format=figType)

    ####################
    # Parseing methods #
    ####################
    
    def _parse_ax_fig(self, ax, fig):
        """ Parse and return ax and fig parameters. """
        if ax is None:
            ax = plt.gca()
        if fig is None:
            fig = plt.gcf()
        return ax, fig
        
    def _parse_ax(self, ax):
        """ Parse and return ax and fig parameters. """
        if ax is None:
            ax = plt.gca()
        return ax
        
    #################
    # Check methods #
    #################       
    def _check_basemap(self):
        """ Check that basemap is not None, raise ValueError if it is. """
        if self.basemap is None:
            raise ValueError('no basemap plotted')