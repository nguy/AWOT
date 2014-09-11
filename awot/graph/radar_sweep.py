"""
awot.graph.radar_sweep
=========================

A group of scripts create various plots of native coordinate data 
collected by the NOAA P-3 tail Doppler radar. 

Created by Nick Guy.

"""
# HISTORY::
#  6 Mar 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
# 10 Sep 2014 - Refactored to a class sctructure to fit AWOT flow
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

from .common import plot_polar_contour
from .coord_transform import radar_coords_to_cart_track_relative, \
       radar_coords_to_cart_earth_relative, radar_coords_to_cart_aircraft_relative
       
#===============================================================
# BEGIN FUNCTIONS
#**===============================================================
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
            print "Use the radar_grid library (RadarGridPlot class)"
            return
            
        # Now initialize the RadarGridPlot Class
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