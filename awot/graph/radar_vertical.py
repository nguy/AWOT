"""
awot.graph.radar_vertical
=========================

A group of scripts create various plots of native coordinate data 
collected by the NOAA P-3 tail Doppler radar. 

Author:
    06 Mar 2014 - Created by Nick Guy, NOAA/NSSL/WRDD, NRC.
    10 Sep 2014 - Refactored to a class sctructure to fit AWOT flow

"""
# FUNCTIONS::
# polar_sweep - Plot polar coordinate data on polar coordinate axis
# polar_sweep_grid - Plotting transformed data to Cartesian output
# sweep_to_Cart - Polar coordinates transformed to Cartesian
# sweep_aircraft_relative - Polar coord data transformed to aircraft-relative frame
# sweep_track_relative - Polar coord data transformed to track-relative frame
# sweep_earth_relative - Polar coord data transformed to earth-relative frame
#-------------------------------------------------------------------
# Load the needed packages
import numpy as np
import matplotlib.pyplot as plt
from  matplotlib.dates import DateFormatter, date2num
import scipy.ndimage as scim

from .common import plot_polar_contour, get_masked_data
from .common import _get_start_datetime, _get_end_datetime
from .common import contour_date_ts
from .coord_transform import radar_coords_to_cart_track_relative, \
       radar_coords_to_cart_earth_relative, radar_coords_to_cart_aircraft_relative
       
       
# Define various constants that may be used for calculations
RE = 6371.  # Earth radius (average)
#===============================================================
# BEGIN FUNCTIONS
#**===============================================================
########################
# Vertical Sweep Class #
########################

class RadarSweepPlot(object):
    """
    To create a RadarPlot instance:
    new_instance = RadarSweepPlot() or new_instance = RadarSweepPlot(AirborneInstance)
    
    Notable attributes
    ------------------
    
    A basemap call can be passed if created externally.
    """
    def __init__(self, airborne, basemap=None, instrument=None):
        '''Intitialize the class to create plots'''
        # Check the instrument to see how to import airborne class
        if instrument is None:
            print "Trying tail Doppler radar, please specify instrument type!"
            instrument = 'tdr_sweep'
        elif instrument == 'tdr_sweep':
            radar_data = airborne.tdr_sweep_radar_data
        elif instrument == 'tdr_grid':
            print "Use the radar_grid library (RadarRadarHorizontalPlotPlot class)"
            return
            
        # Now initialize the RadarHorizontalPlot Class
        self.longitude = radar_data.longitude
        self.latitude = radar_data.latitude
        self.altitude = radar_data.altitude
        self.fields = radar_data.fields
        self.range =radar_data.range
        self.rotation = radar_data.rotation
        self.drift = radar_data.drift
        self.heading = radar_data.heading
        self.pitch = radar_data.pitch
        self.roll = radar_data.roll
        self.tilt = radar_data.tilt
        
        if basemap != None:
            self.basemap = basemap
        try:
            self.basemap = airborne.basemap
        except:
            print "Warning: No basemap instance"
    
    ############################
    # Plotting generate method #
    ############################

    def plot_to_grid(self, Xcoord, Ycoord, Values,
               vmin=-24., vmax=64., cmap='jet',
               xlims=None, ylims=None,
               xlab=None, ylab=None, title=None,
               grid_on=True, plot_in_km=True,
               cb_flag=True, cb_label=None,
               cb_orient=None, cb_pad=None, cb_width=None,
               ax=None, fig=None):
               
        """Plot a sweep of native (polar) projected on a Cartesian plane.
        This method does the actual plotting and is called from the various other
        projection methods.
        
        Parameters::
        ----------
        Xcoord : float array
            X-axis array to use for plot
        Ycoord : float array
            Y-axis array to use for plot
        Values : float array
            2D Data array to plot (Ydims,Xdims)
        vmin : float
            Minimum value to display
        vmax : float
            Maximum value to display
        cmap : string
            Matplotlib color map to use
            
        xlims : tuple
            X-axis limits (min,max)
        ylims : tuple
            Y-axis limits (min,max)
        xlab : str
            X-axis label
        ylab : str
            Y-axis label
        title : str
            Plot title, if None then it is omitted
            
        grid_on : boolean
            True will add a grid to plot, False leaves off
        plot_in_km : boolean
            If True (default) the X, Y coordinates are converted to km instead of m
            
        cb_flag : boolean
            True to add colorbar, False does not
        cb_orient : str
            Location of colorbar if True, either 'horizontal' or 'vertical'
        cb_pad : str
            Pad to move colorbar, fraction in the form ".05", 
            pos is to right for righthand location
        cb_label : str
            Colorbar label
                
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the plot to. None will use the current figure.
        
        Notes::
        -----
        Plotting convention is to project the data onto a 2D surface 
        looking from the back of aircraft forward.
        """
    
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
        
        if plot_in_km:
            Xcoord = Xcoord / 1000.
            Ycoord = Ycoord / 1000.
            xunit = 'km'
            yunit = 'km'
        else:
            xunit = 'm'
            yunit = 'm'
    
        # Plot the data
        p = ax.pcolormesh(Xcoord, Ycoord, Values, cmap=cmap, vmin=vmin, vmax=vmax)
               
        # Set the title
        if title == None:
            pass
        else:
            ax.set_title(title)
            
        # Set the axes limits
        if xlims == None:
            xlims = (-1.05 * Xcoord.max(), 1.05 * Xcoord.max())
        ax.set_xlim(xlims)
        if ylims == None:
            ylims = (-10.,30.)
        ax.set_ylim(ylims)

        print xlims, ylims
        # Set the axes labels
        if xlab == None:
            pass
        else:
            ax.set_xlabel(xlab+' ('+xunit+')')
        if ylab == None:
            pass
        else:
            ax.set_ylabel(ylab+' ('+yunit+')')
        
        # Check if grid lines are desired
        if grid_on:
            ax.grid()
       
        # Set the colorbar if desired
        if cb_width is None:
            cb_width = .10
        if cb_orient is None:
            cb_orient = 'horizontal'
        if cb_flag:
            if cb_pad is None:
                if cb_orient == 'horizontal':
                    cb_pad = .15
                elif cb_orient == 'vertical':
                    cb_pad = .05
                
            cb = plt.colorbar(mappable=p, orientation=cb_orient, pad=cb_pad,
                              fraction=cb_width)
        
            if cb_label == None:
                pass
            else:
                cb.set_label(cb_label)
            
        return p

    ######################
    # Sweep plot modules #
    ######################
    
    def plot_polar_sweep(self, field, nlevs=30, mask_procedure=None, mask_tuple=None,
               vmin=None, vmax=None, cmap='gist_ncar',
               rng_rings=None, rot_angs=None,
               cb_flag=True, cb_orient=None, cb_label=None,
               title=None, ax=None, fig=None):
        """
        Plot a sweep of native (polar) coordinate radar data on polar format plot
        
        Parameters::
        ----------
        field : str
            Variable name (e.g. 'dBZ') to use in plot
        nlevs : int
            Number of contour levels to plot
        mask_procedure : str
            String indicating how to apply mask via numpy, possibilities are:
            'less', 'less_equal', 'greater', 'greater_equal', 'equal', 'inside', 'outside'
        mask_tuple : (str, float[, float])
            Tuple containing the field name and value(s) below which to mask
            field prior to plotting, for example to mask all data where
        vmin : float
            Minimum value to display
        vmax : float
            Maximum value to display
        cmap : string
            Matplotlib color map to use
            
        cb_flag : boolean
            True turns on colorbar, False no colorbar
        cb_orient : str
            Orientation of colorbar ('vertical' or 'horizontal')
        cb_label : str
            Colorbar label, None will use a default label generated from the
            field information.
        title : str
            Title to label plot with, if none given it is omitted
            
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the plot to. None will use the current figure.

        Notes::
        -----
        Variables used to calculate polar coordinate positions:
        Rotation angle with respect to instrument [degrees]
        Range along ray [m]
        
        Defaults are established during DYNAMO project analysis
        """
    
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
        
        # Get variable
        Var, Data = self._get_variable_dict_data(field)
        
        if mask_procedure != None:
            Data = get_masked_data(Data, mask_procedure, mask_tuple)
        
        if vmin is None:
            vmin = Data.min()
            
        if vmax is None:
            vmax = Data.max()
        
        # Set variable for calculation
        range, rotation = self._get_2D_var(self.rotation['data'][:])
            
        # Plot the polar coordinate radar data
        p = plot_polar_contour(Data, rotation, range,
                           nlevs=nlevs, vmin=vmin, vmax=vmax, cmap=cmap)

        # Set the range and turn grid on
        ax.set_rmax(1.05 * range.max())
        ax.grid(True)
               
        # Set the title
        if title == None:
            pass
        else:
            ax.set_title(title)
        
        # Set the range rings if set
        if rng_rings == None:
            pass
        else:
           plt.rgrids(rng_rings)
        
        # Set the rotation angles (theta) if set
        if rot_angs == None:
            pass
        else:
           plt.thetagrids(rot_angs)
       
        # Set the colorbar if desired
        if cb_flag:
            cb = plt.colorbar(mappable=p,orientation=cb_orient)
            if cb_label == None:
                cb_label = Var['long_name'] +' (' + Var['units'] + ')'
            else:
                cb.set_label(cb_label)
    
        return
        
    ##########
    
    def plot_aircraft_relative(self, field, mask_procedure=None, mask_tuple=None,
               vmin=-24., vmax=64., cmap=None,
               xlims=None, ylims=None,
               xlab='Distance from aircraft', ylab='Altitude', title=None,
               grid_on=True, plot_in_km=True,
               cb_flag=True, cb_orient=None, cb_pad=None, cb_label=None,
               ax=None, fig=None):
        """
        Project native (polar) coordinate radar sweep data onto an
        aircraft-relative Cartesian coordinate grid.
        See Lee et al. (1994) Journal of Atmospheric and Oceanic Technology 
        for methodology and definitions.
        
        Parameters::
        ----------
        field : str
            Variable name (e.g. 'dBZ') to use in plot
        mask_procedure : str
            String indicating how to apply mask via numpy, possibilities are:
            'less', 'less_equal', 'greater', 'greater_equal', 'equal', 'inside', 'outside'
        mask_tuple : (str, float[, float])
            Tuple containing the field name and value(s) below which to mask
            field prior to plotting, for example to mask all data where
        data_proj : str
            Which direction the data is collected [ 'fore' or 'aft' ]
            Needed to display the data in forward-looking convention
        vmin : float
            Minimum value to display
        vmax : float
            Maximum value to display
        cmap : string
            Matplotlib color map to use
            
        xlims : tuple
            X-axis limits (min,max)
        ylims : tuple
            Y-axis limits (min,max)
        xlab : str
            X-axis label
        ylab : str
            Y-axis label
        title : str
            Plot title, if None then it is omitted
            
        grid_on : boolean
            True will add a grid to plot, False leaves off
        plot_in_km : boolean
            If True (default) the X, Y coordinates are converted to km instead of m
            
        cb_flag : boolean
            True to add colorbar, False does not
        cb_orient : str
            Location of colorbar if True, either 'horizontal' or 'vertical'
        cb_pad : str
            Pad to move colorbar, fraction in the form ".05", 
            pos is to right for righthand location
        cb_label : str
            Colorbar label
        
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the plot to. None will use the current figure.
                     
        Notes::
        ------
        Variables used in coordinate transformation calculations
        Rotation angle with respect to instrument [degrees]
        Range along ray [m]
        Tilt angle with respect to platform [degrees]
        
        Plotting convention is to project the data onto a 2D surface 
        looking from the back of aircraft forward.
        
        X,Y,Z coordinates are a direct projection from polar to Cartesian coordinates.
        
        This mapping does NOT take into account corrections for 
        roll, pitch, or drift of the aircraft.
        """
    
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
        
        # Get variable
        Var, Data = self._get_variable_dict_data(field)
        
        if mask_procedure != None:
            Data = get_masked_data(Data, mask_procedure, mask_tuple)
        
        if vmin is None:
            vmin = Data.min()
            
        if vmax is None:
            vmax = Data.max()
        
        # Set variable for calculation
        range, rotation = self._get_2D_var(self.rotation['data'][:])
        range, tilt = self._get_2D_var(self.tilt['data'][:])
            
        X, Y, Z = radar_coords_to_cart_aircraft_relative(range / 1000.0, rotation, tilt)
        
        if cb_label == None:
            cb_label = field +' (' + Var['units'] + ')'

        
        p = self.plot_to_grid(X, Z, Data, 
                         vmin=vmin, vmax=vmax, cmap=cmap,
                         xlims=xlims, ylims=ylims,
                         xlab=xlab, ylab=ylab, title=title,
                         grid_on=grid_on,
                         cb_flag=cb_flag, cb_orient=cb_orient, 
                         cb_pad=cb_pad, cb_label=cb_label,
                         ax=ax, fig=fig)
        return p
        
    ##########

    def plot_track_relative(self, field, mask_procedure=None, mask_tuple=None,
               vmin=-24., vmax=64., cmap=None,
               xlims=None, ylims=None,
               xlab='Distance from track', ylab='Altitude', title=None,
               grid_on=True, plot_in_km=True,
               cb_flag=True, cb_orient=None, cb_pad=None, cb_label=None,
               ax=None, fig=None):

        """
        Project native (polar) coordinate radar sweep data onto 
        track-relative Cartesian coordinate grid.
        See Lee et al. (1994) Journal of Atmospheric and Oceanic Technology 
        for methodology and definitions.
        
        Parameters::
        ----------
        field : str
            Variable name (e.g. 'dBZ') to use in plot
        mask_procedure : str
            String indicating how to apply mask via numpy, possibilities are:
            'less', 'less_equal', 'greater', 'greater_equal', 'equal', 'inside', 'outside'
        mask_tuple : (str, float[, float])
            Tuple containing the field name and value(s) below which to mask
            field prior to plotting, for example to mask all data where
        data_proj : str
            Which direction the data is collected [ 'fore' or 'aft' ]
            Needed to display the data in forward-looking convention
        vmin : float
            Minimum value to display
        vmax : float
            Maximum value to display
        cmap : string
            Matplotlib color map to use
            
        xlims : tuple
            X-axis limits (min,max)
        ylims : tuple
            Y-axis limits (min,max)
        xlab : str
            X-axis label
        ylab : str
            Y-axis label
        title : str
            Plot title, if None then it is omitted
            
        grid_on : boolean
            True will add a grid to plot, False leaves off
        plot_in_km : boolean
            If True (default) the X, Y coordinates are converted to km instead of m
            
        cb_flag : boolean
            True to add colorbar, False does not
        cb_orient : str
            Location of colorbar if True, either 'horizontal' or 'vertical'
        cb_pad : str
            Pad to move colorbar, fraction in the form ".05", 
            pos is to right for righthand location
        cb_label : str
            Colorbar label
        
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the plot to. None will use the current figure.
                     
        Notes::
        ------
        Variables used in coordinate transformation calculations
        Rotation angle with respect to instrument [degrees]
        Range along ray [m]
        Tilt angle with respect to platform [degrees]
        Roll angle of aircraft [degrees], right-wing down positive
        Drift angle [degrees], between heading and track
        Pitch angle [degrees], nose up positive
        
        Plotting convention is to project the data onto a 2D surface 
        looking from the back of aircraft forward.  
        
        X,Y,Z coordinates are a rotation of the data about the aircraft track, 
        following a direct projection from polar to Cartesian coordinates.
        
        This mapping corrects for roll, pitch, and drift of the aircraft.
        
        This is considered a leveled, heading-relative coordinate system.
        """
    
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
        
        # Get variable
        Var, Data = self._get_variable_dict_data(field)
        
        if mask_procedure != None:
            Data = get_masked_data(Data, mask_procedure, mask_tuple)
        
        if vmin is None:
            vmin = Data.min()
            
        if vmax is None:
            vmax = Data.max()
        
        # Set variable for calculation
        range, rotation = self._get_2D_var(self.rotation['data'][:])
        range, tilt = self._get_2D_var(self.tilt['data'][:])
        range, roll = self._get_2D_var(self.roll['data'][:])
        range, drift = self._get_2D_var(self.drift['data'][:])
        range, pitch = self._get_2D_var(self.pitch['data'][:])
            
        X, Y, Z = radar_coords_to_cart_track_relative(range / 1000.0, rotation, roll, 
                                                      drift, tilt, pitch)
        
        if cb_label == None:
            cb_label = field +' (' + Var['units'] + ')'
        
        p = self.plot_to_grid(X, Z, Data, 
                         vmin=vmin, vmax=vmax, cmap=cmap,
                         xlims=xlims, ylims=ylims,
                         xlab=xlab, ylab=ylab, title=title,
                         grid_on=grid_on,
                         cb_flag=cb_flag, cb_orient=cb_orient, 
                         cb_pad=cb_pad, cb_label=cb_label,
                         ax=ax, fig=fig)
        return p
        
    ##########
        
    def plot_earth_relative(self, field, mask_procedure=None, mask_tuple=None,
               vmin=-24., vmax=64., cmap=None,
               xlims=None, ylims=None,
               xlab='Distance from aircraft', ylab='Altitude', title=None,
               grid_on=True, plot_in_km=True,
               cb_flag=True, cb_orient='horizontal', cb_pad=None, cb_label=None,
               ax=None, fig=None):
        """
        Project native (polar) coordinate radar sweep data onto 
        earth-relative Cartesian coordinate grid.
       See Lee et al. (1994) Journal of Atmospheric and Oceanic Technology 
       for methodology and definitions.
        
        Parameters::
        ----------
        field : str
            Variable name (e.g. 'dBZ') to use in plot
        mask_procedure : str
            String indicating how to apply mask via numpy, possibilities are:
            'less', 'less_equal', 'greater', 'greater_equal', 'equal', 'inside', 'outside'
        mask_tuple : (str, float[, float])
            Tuple containing the field name and value(s) below which to mask
            field prior to plotting, for example to mask all data where
        data_proj : str
            Which direction the data is collected [ 'fore' or 'aft' ]
            Needed to display the data in forward-looking convention
        vmin : float
            Minimum value to display
        vmax : float
            Maximum value to display
        cmap : string
            Matplotlib color map to use
            
        xlims : tuple
            X-axis limits (min,max)
        ylims : tuple
            Y-axis limits (min,max)
        xlab : str
            X-axis label
        ylab : str
            Y-axis label
        title : str
            Plot title, if None then it is omitted
            
        grid_on : boolean
            True will add a grid to plot, False leaves off
        plot_in_km : boolean
            If True (default) the X, Y coordinates are converted to km instead of m
            
        cb_flag : boolean
            True to add colorbar, False does not
        cb_orient : str
            Location of colorbar if True, either 'horizontal' or 'vertical'
        cb_pad : str
            Pad to move colorbar, fraction in the form ".05", 
            pos is to right for righthand location
        cb_label : str
            Colorbar label
        
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the plot to. None will use the current figure.
                     
        Notes::
        ------
        ------
        Variables used in coordinate transformation calculations
        Rotation angle with respect to instrument [degrees]
        Range along ray [m]
        Tilt angle with respect to platform [degrees]
        Roll angle of aircraft [degrees], right-wing down positive
        Heading angle [degrees], clockwise from North
        Pitch angle [degrees], nose up positive
        
        Plotting convention is to project the data onto a 2D surface 
        looking from the back of aircraft forward.  
        
        X,Y,Z coordinates are a rotation of the data about an 
        earth-relative azimuth, following a direct projection from polar to 
        Cartesian coordinates.
        
        This mapping corrects for roll, pitch, and drift of the aircraft.
        
        This is considered a leveled, heading-relative coordinate system.
        """
    
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
        
        # Get variable
        Var, Data = self._get_variable_dict_data(field)
        
        if mask_procedure != None:
            Data = get_masked_data(Data, mask_procedure, mask_tuple)
        
        if vmin is None:
            vmin = Data.min()
            
        if vmax is None:
            vmax = Data.max()
        
        # Set variable for calculation
        range, rotation = self._get_2D_var(self.rotation['data'][:])
        range, tilt = self._get_2D_var(self.tilt['data'][:])
        range, roll = self._get_2D_var(self.roll['data'][:])
        range, pitch = self._get_2D_var(self.pitch['data'][:])
        range, heading = self._get_2D_var(self.heading['data'][:])
            
        X, Y, Z = radar_coords_to_cart_earth_relative(range / 1000.0, rotation, roll, 
                                                      heading, tilt, pitch)
        
        if cb_label == None:
            cb_label = field +' (' + Var['units'] + ')'
        
        p = self.plot_to_grid(X, Z, Data, 
                         vmin=vmin, vmax=vmax, cmap=cmap,
                         xlims=xlims, ylims=ylims,
                         xlab=xlab, ylab=ylab, title=title,
                         grid_on=grid_on,
                         cb_flag=cb_flag, cb_orient=cb_orient, 
                         cb_pad=cb_pad, cb_label=cb_label,
                         ax=ax, fig=fig)
        return p
        
    ####################
    # Get methods #
    #################### 
    
    def _get_variable_dict(self, field):
        '''Get the variable from the fields dictionary'''
        Var = self.fields[field]
       
        return Var
    
    def _get_variable_dict_data(self, field):
        '''Get the variable from the fields dictionary'''
        Var, data = self.fields[field], self.fields[field]['data'][:]
       
        return Var, data
        
    def _get_2D_var(self, Var):
        '''Return a 2D variable 
        Useful for coordinate transformations
        '''
        rg2D, Var2D = np.meshgrid(self.range['data'][:], Var)
        
        return rg2D, Var2D

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
        
    def _parse_fig(self, fig):
        """ Parse and return ax and fig parameters. """
        if fig is None:
            fig = plt.gcf()
        return fig
        

#######################
# Vertical Plot Class #
#######################

class RadarVerticalPlot(object):
    """
    To create a RadarHorizontalPlot instance:
    new_instance = RadarHorizontalPlot() or new_instance = RadarHorizontalPlot(AirborneInstance)
    
    Notable attributes
    ------------------
    
    A basemap call can be passed if created externally.
    """
    def __init__(self, airborne, basemap=None, instrument=None):
        '''Intitialize the class to create plots'''
        # Check the instrument to see how to import airborne class
        if instrument is None:
            print "Trying tail Doppler radar, please specify instrument type!"
            instrument = 'tdr_grid'
        elif instrument == 'tdr_grid':
            radar_data = airborne.tdr_radar_data
        elif instrument == 'lf':
            radar_data = airborne.lf_radar_data
        elif instrument == 'rasta':
            radar_data = airborne.rasta_radar_data
        elif instrument == 'ground':
            print "Connected using PyArt"
            radar_data = airborne.ground_radar_data
            return
        elif instrument == 'tdr_sweep':
            print "Use the radar_sweep library (RadarSweepPlot class)"
            return
            
        # Now initialize the RadarHorizontalPlot Class
        self.longitude = radar_data['longitude']
        self.latitude = radar_data['latitude']
        self.height = radar_data['height']
        self.fields = radar_data['fields']
        
        # Attempt to pull in time if found
        try:
            self.time = radar_data['time']
        except:
            print "Warning: No time variable found"
        
        if basemap != None:
            self.basemap = basemap
        try:
            self.basemap = airborne.basemap
        except:
            print "Warning: No basemap instance"
            
        # Save the airborne class to this class in case cross-section is passed
        self.airborne = airborne
            
    #########################
    # Vertical plot methods #
    #########################

    def plot_cross_section(self, field, start_pt, end_pt, xs_length=500,
                           mask_procedure=None, mask_tuple=None,
                           title=" ", title_size=20,
                           cminmax=(0.,60.), clevs=25, vmin=15., vmax=60.,
                           cmap='gist_ncar', clabel='dBZ',
                           color_bar=True, cb_pad=.05, cb_orient='vertical', cb_tick_int=2,
                           ax=None, fig=None):
        '''
        Plot a cross-section between two points
        
        Parameters::
        ----------
        field : str
            3-D variable (e.g. Reflectivity [dBZ]) to use in plot
        start_pt, end_pt : tuple
            (lat, lon) Tuple of start, end points for cross-section
        xs_length : int
            Number of to use for the cross section
        mask_procedure : str
            String indicating how to apply mask via numpy, possibilities are:
            'less', 'less_equal', 'greater', 'greater_equal', 'equal', 'inside', 'outside'
        mask_tuple : (str, float[, float])
            Tuple containing the field name and value(s) below which to mask
            field prior to plotting, for example to mask all data where
        cminmax : tuple
            (min,max) values for controur levels
        clevs : integer
            Number of contour levels
        vmin : float
            Minimum contour value to display
        vmax : float
            Maximum contour value to display
        clabel : string
            Label for colorbar (e.g. units 'dBZ')

        title : string
            Plot title
        title_size : int
            Font size of title to display
            
        cmap : string
            Matplotlib color map to use
        color_bar : boolean
            True to add colorbar, False does not
        cb_pad : str
            Pad to move colorbar, in the form "5%", pos is to right for righthand location
        cb_loc : str
            Location of colorbar, default is 'right', also available: 
            'bottom', 'top', 'left'
        cb_tick_int : int
            Interval to use for colorbar tick labels, higher number "thins" labels
                
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the plot to. None will use the current figure.
        '''
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
            
        # Return masked or unmasked variable
        Var, Data = self._get_variable_dict_data(field)
        if mask_procedure != None:
            Data = get_masked_data(Data, mask_procedure, mask_tuple)
            
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
            xs_data[:,nlev] = scim.map_coordinates(Data[nlev,:,:],
                                                  np.vstack((xsY, xsX)),
                                                  prefilter=False)#, mode='nearest')
            
        # Calculate the distance array along the cross-section
        Xdist = np.absolute((np.pi * RE / 180.) * (xslon - xslon[0]))
        Ydist = np.absolute((np.pi * RE / 180.) * (xslat - xslat[0]))
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
        Ht2D, Dist2D = np.meshgrid(self.height['data'][:], xsDist)
        
        p = ax.pcolormesh(Dist2D, Ht2D, np.ma.masked_less_equal(xs_data, -800.), 
                          vmin=vmin, vmax=vmax, cmap=cmap)
                          
        ax.set_xlabel('Distance along track (km)')
        ax.set_ylabel(' Altitude (km)')
        
        # Add title
        ax.set_title(title, fontsize=title_size)

        # Add Colorbar
        if color_bar:
            cbStr = Var['long_name'] +' ('+ Var['units'] +')'
            cb = fig.colorbar(p, orientation=cb_orient, pad=cb_pad)#,ticks=clevels)
            cb.set_label(cbStr)
            # Set the number of ticks in the colorbar based upon number of contours
            tick_locator = ticker.MaxNLocator(nbins=int(clevs/cb_tick_int))
            cb.locator = tick_locator
            cb.update_ticks()
    
        # Add title
        ax.set_title(title, fontsize=title_size)
    
    def contour_timeseries(self, field,
                           mask_procedure=None, mask_tuple=None,
                           ptype='pcolormesh',
                           cminmax=(0.,60.), clevs=25, vmin=15., vmax=60.,
                           cmap='gist_ncar',
                           color_bar=True, cb_orient='vertical',
                           cb_pad=.05, cb_tick_int=2,
                           cb_label=None, 
     
                           dForm='%H:%M',tz=None, xdate=True,
                           date_MinTicker='minute',
                           other_MajTicks=None, other_MinTicks=None,
                           other_min=None, other_max=None,
                
                           start_time=None, end_time=None,
                 
                           title=None,
                           xlab=' ', xlabFontSize=16, xpad=7,
                           ylab=' ', ylabFontSize=16, ypad=7,
                           ax=None, fig=None):
        """
        Wrapper function to produce a contoured time series plot of variable indicated
        
        Parameters::
        ----------
        field : float
            Variable to plot as time series
        mask_procedure : str
            String indicating how to apply mask via numpy, possibilities are:
            'less', 'less_equal', 'greater', 'greater_equal', 'equal', 'inside', 'outside'
        mask_tuple : (str, float[, float])
            Tuple containing the field name and value(s) below which to mask
            field prior to plotting, for example to mask all data where
        cminmax : tuple
            (min,max) values for controur levels
        clevs : integer
            Number of contour levels
        vmin : float
            Minimum contour value to display
        vmax : float
            Maximum contour value to display
        
        ptype : str
            Type of plot to make, takes 'plot', 'contour', or 'pcolormsh'
        cmap : string
            Matplotlib color map to use
        color_bar : boolean
            True to add colorbar, False does not
        cb_pad : str
            Pad to move colorbar, in the form "5%", pos is to right for righthand location
        cb_loc : str
            Location of colorbar, default is 'right', also available: 
            'bottom', 'top', 'left'
        cb_tick_int : int
            Interval to use for colorbar tick labels, higher number "thins" labels
        cb_label : string
            Label for colorbar (e.g. units 'dBZ')
            
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
        other_min : float
            Minimum value for non-date axis
        other_max : float
            Maximum value for non-date axis
    
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
        fig : Figure
            Figure to add the plot to. None will use the current figure.
        """
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
        
        # Return masked or unmasked variable
        # Subsetted if desired
        Var, tsub, Data = self._get_variable_dict_data_time_subset(field, start_time, end_time)
        if mask_procedure != None:
            Data = get_masked_data(Data, mask_procedure, mask_tuple)
            
        # Create contour level array
        clevels = np.linspace(cminmax[0], cminmax[1], clevs)
                                     
        tSub2D, Ht2D = np.meshgrid(date2num(tsub), self.height['data'][:])
        
        # Plot the time series
        ts = contour_date_ts(tSub2D, Ht2D, Data.T, 
                vmin=vmin, vmax=vmax, clevs=clevs,
                color_bar=color_bar, cb_orient=cb_orient,
                cb_pad=cb_pad, cb_tick_int=cb_tick_int,
                cb_label=cb_label,
                dForm=dForm,tz=tz, xdate=xdate, 
                date_MinTicker=date_MinTicker,
                other_MajTicks=other_MajTicks, other_MinTicks=other_MinTicks,
                other_min=other_min, other_max=other_max,
                title=title,
                xlab=xlab, xlabFontSize=xlabFontSize, xpad=xpad,
                ylab=ylab, ylabFontSize=ylabFontSize, ypad=ypad,
                ax=ax, fig=fig)
                
        return
    
    ####################
    # Get methods #
    #################### 
    
    def _get_variable_dict(self, field):
        '''Get the variable from the fields dictionary'''
        Var = self.fields[field]
       
        return Var
    
    def _get_variable_dict_data(self, field):
        '''Get the variable from the fields dictionary'''
        Var, data = self.fields[field], self.fields[field]['data'][:]
       
        return Var, data
    
    def _get_variable_dict_data_time_subset(self, field, start_time, end_time):
        '''
        Get the variable from the fields dictionary.
        Subset the time when in time series format
        '''
        Var, data = self.fields[field], self.fields[field]['data'][:]
        
        # Check to see if time is subsetted
        dt_start = _get_start_datetime(self.time, start_time)
        dt_end = _get_end_datetime(self.time, end_time)

        tsub = self.time[(self.time >= dt_start) & (self.time <= dt_end)]
        datasub = data[(self.time >= dt_start) & (self.time <= dt_end)]
       
        return Var, tsub, datasub
        
    def _get_lat_index(self, value):
        '''Calculate the exact index position within latitude array'''
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
        
    def _parse_fig(self, fig):
        """ Parse and return ax and fig parameters. """
        if fig is None:
            fig = plt.gcf()
        return fig

    