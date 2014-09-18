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

from .common import find_nearest_indices, get_masked_data
#import general.gplot as gp

# Define various constants that may be used for calculations
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
            print "Currently not supported.  Use PyArt to connect here"
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
                cmap='gist_ncar', color_bar=True, cb_pad="5%", cb_loc='right',
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
#    cs = m.contourf(x,y,Z[Zind,:,:],cmap=cm.s3pcpn)#clevs=clevels,
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
            tick_locator = ticker.MaxNLocator(nbins=int(clevs/2))
            cb.locator = tick_locator
            cb.update_ticks()
    
        # Add title
        ax.set_title(title, fontsize=title_size)

        return
        
    ##########
    
    def overlay_wind_vector(self, height_level=2., 
                                vtrim=4, vlw=1.3, vhw=2.5, vscale=400,
                                refVec=True, refU=10., refUposX=1.05, refUposY=1.015):
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
        U = field['Uwind']['data'][Htind,:,:]
        V = field['Vwind']['data'][Htind,:,:]
        lon = self.longitude['data'][:]
        lat = self.latitude['data'][:]
        uproj, vproj, xx, yy = self.basemap.transform_vector(U, V,
                                        lon, lat ,
                                        len(lon)/vtrim, len(lat)/vtrim,
                                        returnxy=True, masked=True)
        
        # Overplot the vectors
        Q = self.basemap.quiver(xx, yy, uproj, vproj, scale=vscale, 
                                headwidth=vhw, linewidths=vlw)
                            
        # Make a quiver key to attach to figure.
        qkLab = str("%4.1f"%refU) + 'm/s at '+str("%4.1f"%Ucoord)+' km'
        qk = ax.quiverkey(Q, refUposX, refUposY, refU, qkLab)#, fontproperties={'weight': 'bold'})
    
        return
        
    ##########
    
    def plot_lf(self, field=None, mask_procedure=None, mask_tuple=None, 
                cminmax=(0.,60.), clevs=25, vmin=15., vmax=60., clabel='dBZ',
                title=" ", title_size=20, cmap='gist_ncar',
                color_bar=True, cb_pad="5%", cb_loc='right',
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
            tick_locator = ticker.MaxNLocator(nbins=int(clevs/2))
            cb.locator = tick_locator
            cb.update_ticks()
    
        # Add title
        ax.set_title(title, fontsize=title_size)

    ##########
    
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