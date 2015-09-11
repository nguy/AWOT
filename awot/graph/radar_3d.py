"""
awot.graph.radar_3d
=========================

A group of scripts to create various 3-dimensional plots from 
data collected by the NOAA P-3 tail Doppler radar. 

Note that is experimental and not fully developed.

Author: 
    04 Sep 2014 Created by Nick Guy, NOAA/NSSL/WRDD, NRC.

"""
#-------------------------------------------------------------------
# Load the needed packages
from mpl_toolkits.mplot3d import axes3d
import numpy as np

from .common import find_nearest_indices, get_masked_data
#===============================================================
# BEGIN FUNCTIONS
#===============================================================
class Radar3DPlot(object):
    """
    To create a RadarHorizontalPlot instance:
    new_instance = RadarHorizontalPlot() or new_instance = RadarHorizontalPlot(AirborneInstance)
    
    Notable attributes
    ------------------
    
    """
    def __init__(self, airborne, instrument=None):
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
            
        # Now initialize the RadarHorizontalPlot Class
        self.longitude = radar_data['longitude']
        self.latitude = radar_data['latitude']
        self.height = radar_data['height']
        self.fields = radar_data['fields']
    
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
        """
        Read in data from NetCDF file containing P3 flight level data created
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