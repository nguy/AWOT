"""
awot.graph.radar_grid
=========================

A group of scripts create various plots of gridded products from 
data collected by the NOAA P-3 tail Doppler radar. 

Created by Nick Guy.

"""
# HISTORY::
# 04 Sep 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
#               Refactored from various scripts
#-------------------------------------------------------------------
# Load the needed packages
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import ticker
import numpy as np
import scipy.ndimage as scim

from .common import find_nearest_indices, get_masked_data
#import general.gplot as gp

# Define various constants that may be used for calculations
RE = 6371.  # Earth radius (average)
#===============================================================
# BEGIN FUNCTIONS
#===============================================================
class RadarGridPlot(object):
    """
    To create a RadarGridPlot instance:
    new_instance = RadarGridPlot() or new_instance = RadarGridPlot(AirborneInstance)
    
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
        elif instrument == 'ground':
            print "Connected using PyArt"
            radar_data = airborne.ground_radar_data
            return
        elif instrument == 'tdr_sweep':
            print "Use the radar_sweep library (RadarSweepPlot class)"
            return
            
        # Now initialize the RadarGridPlot Class
        self.longitude = radar_data['longitude']
        self.latitude = radar_data['latitude']
        self.height = radar_data['height']
        self.fields = radar_data['fields']
        
        if basemap != None:
            self.basemap = basemap
        try:
            self.basemap = airborne.basemap
        except:
            print "Warning: No basemap instance"

    ###########################
    # Horizontal plot modules #
    ###########################
    
    def plot_ppi(self, field, ppi_height=2., mask_procedure=None, mask_tuple=None, 
                cminmax=(0.,60.), clevs=25, vmin=15., vmax=60., clabel='dBZ',
                title=" ", title_size=20, 
                cmap='gist_ncar',
                color_bar=True, cb_pad="5%", cb_loc='right', cb_tick_int=2,
                ax=None, fig=None):
        """Produce a CAPPI (constant altitude plan position indicator) plot
        using the Tail Doppler Radar data.
    
        Parameters::
        ----------
        field : str
            3-D variable (e.g. Reflectivity [dBZ]) to use in plot
        ppi_height : float
            Height at which to plot the horizontal field
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

        Notes::
        -----
        Defaults are established during DYNAMO project analysis
        """
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
        
        # Grab the variable dictionary of interest to plot
        Var = self.fields[field]
            
        # Return masked or unmasked variable
        Var, Data = self._get_variable_dict_data(field)
        if mask_procedure != None:
            Data = get_masked_data(Data, mask_procedure, mask_tuple)
    
        # Create contour level array
        clevels = np.linspace(cminmax[0], cminmax[1], clevs)
    
        # Find the closest vertical point 
        Zind = find_nearest_indices(self.height['data'][:], ppi_height)
        
        # Convert lats/lons to 2D grid
        Lon2D, Lat2D = np.meshgrid(self.longitude['data'][:], self.latitude['data'][:])
        # Convert lats/lons to map projection coordinates
        x, y = self.basemap(Lon2D, Lat2D)
        
        #Plot the contours
#    cs = m.contourf(x,y,Z[Zind,:,:],clevels,vmin=vmin,vmax=vmax,cmap=cmap)
        cs = self.basemap.pcolormesh(x, y, Data[Zind,:,:], 
                                     vmin=vmin, vmax=vmax, cmap=cmap)
#    plt.colors.colormap.set_under('white')

        # Add Colorbar
        if color_bar:
            cbStr = Var['long_name'] + ' at ' + str("%4.1f" % ppi_height) + self.height['units']
            cb = self.basemap.colorbar(cs, location=cb_loc, pad=cb_pad)#,ticks=clevels)
            cb.set_label(cbStr)
            # Set the number of ticks in the colorbar based upon number of contours
            # and the tick interval selected
            tick_locator = ticker.MaxNLocator(nbins=int(clevs/cb_tick_int))
            cb.locator = tick_locator
            cb.update_ticks()
    
        # Add title
        ax.set_title(title, fontsize=title_size)

        return
        
    ##########
    
    def overlay_wind_vector(self, height_level=2., 
                                vtrim=4, vlw=1.3, vhw=2.5, vscale=400,
                                refVec=True, refU=10., refUposX=1.05, refUposY=1.015,
                                qcolor='k',
                                ax=None, fig=None):
        """Overlays a 2-D wind field at specified height onto map
    
        Parameters::
        ----------
        height_level : float
            Height level for winds in horizontal plot
        vtrim : float
            The number of vector arrows will be thinned by this factor
        vlw : float
            Vector arrow linewidth
        vhw : float
            Vector arrow headwidth
        vscale : int
            Vector arrow scale (smaller = longer arrow)
        refVec : boolean
            True to attach a reference vector, False to not
        refU : float
            Magnitude of reference vector [m/s]
        refU_Xpos : float
            X position of reference vector in axes-relative coordinates
        refU_Ypos : float
            Y position of reference vector in axes-relative coordinates
                
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the plot to. None will use the current figure.
  
        Notes::
        -----
        U               = Along aircraft longitudinal axis wind
        V               = Perpendicular aircraft longitudinal axis wind
        """
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
        
        # Find the closest vertical point for desired wind field
        Htind = find_nearest_indices(self.height['data'][:],height_level)
        	
        # transform to porjection grid
        U = self.fields['Uwind']['data'][Htind,:,:]
        V = self.fields['Vwind']['data'][Htind,:,:]
        lon = self.longitude['data'][:]
        lat = self.latitude['data'][:]
        uproj, vproj, xx, yy = self.basemap.transform_vector(U, V,
                                        lon, lat ,
                                        len(lon)/vtrim, len(lat)/vtrim,
                                        returnxy=True, masked=True)
        
        # Overplot the vectors
        Q = self.basemap.quiver(xx, yy, uproj, vproj, scale=vscale, 
                                headwidth=vhw, linewidths=vlw, color=qcolor)
                            
        # Make a quiver key to attach to figure.
        qkLab = str("%4.1f"%refU) + 'm/s at '+str("%4.1f"%height_level)+' km'
        qk = ax.quiverkey(Q, refUposX, refUposY, refU, qkLab)#, fontproperties={'weight': 'bold'})
    
        return
        
    ##########
    
    def plot_lf(self, field=None, mask_procedure=None, mask_tuple=None, 
                cminmax=(0.,60.), clevs=25, vmin=15., vmax=60., clabel='dBZ',
                title=" ", title_size=20, cmap='gist_ncar',
                color_bar=True, cb_pad="5%", cb_loc='right', cb_tick_int=2,
                ax=None, fig=None):
        """Produce a CAPPI (constant altitude plan position indicator) plot
        using the Tail Doppler Radar data.
    
        Parameters::
        ----------
        field : str
            3-D variable (e.g. Reflectivity [dBZ]) to use in plot
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

        Notes::
        -----
        Defaults are established during DYNAMO project analysis
        """
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
        
        # Grab the variable dictionary of interest to plot
        if field is None:
            field = 'dBZ'
            
        # Return masked or unmasked variable
        Var, Data = self._get_variable_dict_data(field)
        if mask_procedure != None:
            Data = get_masked_data(Data, mask_procedure, mask_tuple)
 
        # Create contour level array
        clevels = np.linspace(cminmax[0], cminmax[1], clevs)
                
        # Convert lats/lons to 2D grid
        Lon2D, Lat2D = np.meshgrid(self.longitude['data'][:], self.latitude['data'][:])
        # Convert lats/lons to map projection coordinates
        x, y = self.basemap(Lon2D, Lat2D)
    
        p = self.basemap.pcolormesh(x, y, Data, 
                                    vmin=vmin, vmax=vmax, cmap=cmap)

        # Add Colorbar
        if color_bar:
            cbStr = Var['long_name'] +' '+ self.height['units']
            cb = self.basemap.colorbar(p, location=cb_loc, pad=cb_pad)#,ticks=clevels)
            cb.set_label(cbStr)
            # Set the number of ticks in the colorbar based upon number of contours
            tick_locator = ticker.MaxNLocator(nbins=int(clevs/cb_tick_int))
            cb.locator = tick_locator
            cb.update_ticks()
    
        # Add title
        ax.set_title(title, fontsize=title_size)
    
    ###############
    
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
            
    def plot_line_geo(self, line_lons, line_lats,
                      line_style='r-', lw=3, alpha=0.2,
                      label0=True, label_offset = (0.01, 0.01),
                      ax=None, fig=None, **kwargs):
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
        lw : int
            LInewidth
        alpha : float
            Transparency, 0: transparent, 1: opaque
        label0 : boolean
            True to label the starting point of line, False is no label
        label_offset : [float, float]
            Offset in lon, lat degrees for the bottom left corner of the label
            text relative to the point. A value of None will use 0.01 default.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the plot to. None will use the current figure.
        """
        # parse parameters
        if ax is None:
            ax = self._parse_ax(ax)
        
        x, y = self.basemap(line_lons, line_lats)
        self.basemap.plot(x, y, line_style, lw=lw, alpha=alpha, ax=ax)
                          
        # Overplot the 0 point
        if label0:
            self.plot_point(line_lons[0], line_lats[0], 'ko',
                            label_text='0', label_offset=label_offset, markersize=0.5)
    
    #########################
    # 3-D plot methods #
    #########################
    
    def DPJgrid_3d(self, surf_field,
               surf_min=-5., surf_max=5., surf_cmap='RdBu_r',
               rstride=5, cstride=5,
               plot_contour=False, cont_field=None, ppi_height=2.,
               cminmax=(0.,60.), clevs=25, vmin=15., vmax=60., cmap='gist_ncar', alpha=0.,
               zlims=(-5.,5.), dlat=1.,dlon=1.,
               plot_track=True,
               title=" ", fig=None, ax=None):
        """Read in data from NetCDF file containing P3 flight level data created
        by NOAA AOC.  The NetCDF should be read in the main program and passed
        to this function.
        
        Parameters::
        ----------
        surf_field : string
            Name of field to use for the 3D surface plot
        ppi_height : float
            Height at which to plot the horizontal field

        surf_min : float
            Minimum surface value to display
        surf_max : float
            Maximum surface value to display
        surf_cmap : string
            Matplotlib color map to use
        rstride : int
            Row stride
        cstride : int
            Column stride
        plot_contour : boolean
            True to plot a contour along the axes
        cminmax : tuple
            (min,max) values for controur levels
        clevs : integer
            Number of contour levels
            Number of contour levels
        vmin : float
            Minimum contour value to display
        vmax : float
            Maximum contour value to display
        cmap : string
            Matplotlib color map to use
        alpha : float
            Alpha factor for opacity
        zlims : 2-tuple
            (Min, Max) tuple for Z-axis
  dlat            = Latitudinal spacing on plot
  dlon            = Longitudinal spacing on plot
  proj            = Projection to use for map
  title           = Plot title
  pName           = String name of output file
  pType           = String name of output file format
  figsize         = [x,y] size of figure to create
        plot_track : boolean
            True to overplot the aircraft track (Note must have ingested
            flight level data file as well)        
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the plot to. None will use the current figure.
            
        Notes::
        -----
        Defaults are established during DYNAMO project analysis
        """
        # parse parameters
        fig = self._parse_fig(fig)
        
        # Get variable
        Var, Data = self._get_variable_dict_data(surf_field)
    
        # Create contour level array
        clevels = np.linspace(vmin, vmax, clevs)
    
        # Set up the axes for 3D projection
        ax = fig.gca(projection='3d')
        
        # Convert lats/lons to 2D grid
        Lon2D, Lat2D = np.meshgrid(self.longitude['data'][:],self.latitude['data'][:])

        # Plot the vertical velocity as a surface plot    
        pS = ax.plot_surface(Lat2D, Lon2D, Data, 
                            vmin=surf_min, vmax=surf_max, linewidth=0, alpha=alpha,
                            rstride=rstride, cstride=cstride, cmap=surf_cmap)
                         
#    pW = ax.plot_wireframe(Lat2D,Lon2D,W,rstride=rstride,cstride=cstride,alpha=alf)
        ax.set_xlim(Lat2D.min(), Lat2D.max())
        ax.set_ylim(Lon2D.min(), Lon2D.max())
        ax.set_zlim(zlim)
        ax.set_xlabel('Latitude')
        ax.set_ylabel('Longitude')
        ax.set_zlabel(' Altitude (km)')
#    ax.view_init(20.,
        cb = fig.colorbar(pS, shrink=0.6)
        cb.set_label(Var['long_name'] + Var['units'])#r'(m s$^{-1}$)')
    
        # Plot the horizontal dBZ contour field
        if plot_contour:
            if (con_field != None):
                conVar = self._get_variable_dict(con_field)
    
                # Find the closest vertical point 
                Zind = find_nearest_indices(self.height['data'][:], ppi_height)
            
                ax.contourf(Lon2D, Lat2D, conVar['data'][Zind,:,:], clevels,
                            vmin=vmin, vmax=vmax, cmap=cmap,
                            zdir='z')#,offset=surf_min)
            else:
                print "Need to set con_field and ppi_height"
#    
                             
        if plot_track:
            ax.plot(self.latitude['data'][:], self.longitude['data'][:], 
                    self.height['data'][:]/1000., zdir='z', c='k')

        return
    
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