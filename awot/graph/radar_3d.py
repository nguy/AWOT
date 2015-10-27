"""
awot.graph.radar_3d
===================

A group of scripts to create 3-dimensional plots.

Note that is experimental and not fully developed.

"""

from __future__ import print_function
from mpl_toolkits.mplot3d import axes3d
import numpy as np

from .common import find_nearest_indices, get_masked_data

class Radar3DPlot(object):
    """
    To create a RadarHorizontalPlot instance:
    new_instance = RadarHorizontalPlot() or
    new_instance = RadarHorizontalPlot(AirborneInstance)

    Notable attributes
    ------------------

    """

    def __init__(self, radardata):
        '''Intitialize the class to create plots'''
        self.radar_data = radardata
        # Now initialize the RadarHorizontalPlot Class
        self.longitude = self.radar_data['longitude']
        self.latitude = self.radar_data['latitude']
        self.height = self.radar_data['height']
        self.fields = self.radar_data['fields']

####################
#   Plot methods   #
####################

    def DPJgrid_3d(self, surf_field,
                   surf_min=-5., surf_max=5., surf_cmap='RdBu_r',
                   rstride=5, cstride=5,
                   plot_contour=False, cont_field=None, ppi_height=2.,
                   cminmax=(0., 60.), clevs=25, vmin=15., vmax=60.,
                   cmap='gist_ncar', alpha=0.,
                   zlims=(-5., 5.), dlat=1., dlon=1.,
                   plot_track=True,
                   title=" ", fig=None, ax=None):
        """
        Read in data from AWOT radar instance, create 3D plot.

        Parameters
        ----------
        surf_field : str
            Name of field to use for the 3D surface plot.
        ppi_height : float
            Height at which to plot the horizontal field.

        surf_min : float
            Minimum surface value to display.
        surf_max : float
            Maximum surface value to display.
        surf_cmap : str
            Matplotlib color map to use.
        rstride : int
            Row stride.
        cstride : int
            Column stride.
        plot_contour : bool
            True to plot a contour along the axes.
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
        alpha : float
            Alpha factor for opacity.
        zlims : 2-tuple
            (Min, Max) tuple for Z-axis
        dlat : float
            Latitudinal spacing.
        dlon : float
            Longitudinal spacing.
        proj : str
            Projection to use.
        title : str
            Plot title.
        plot_track : bool
            True to overplot the aircraft track (Note must have ingested
            flight level data file as well).
        ax : Matplotlib axis instance
            Axis to plot on. None will use the current axis.
        fig : Matplotlib figure instance
            Figure to add the plot to. None will use the current figure.
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
        Lon2D, Lat2D = np.meshgrid(self.longitude['data'][
                                   :], self.latitude['data'][:])

        # Plot the vertical velocity as a surface plot
        pS = ax.plot_surface(Lat2D, Lon2D, Data,
                             vmin=surf_min, vmax=surf_max, linewidth=0,
                             alpha=alpha,
                             rstride=rstride, cstride=cstride, cmap=surf_cmap)

#    pW = ax.plot_wireframe(
#        Lat2D,Lon2D,W,rstride=rstride,cstride=cstride,alpha=alf)
        ax.set_xlim(Lat2D.min(), Lat2D.max())
        ax.set_ylim(Lon2D.min(), Lon2D.max())
        ax.set_zlim(zlim)
        ax.set_xlabel('Latitude')
        ax.set_ylabel('Longitude')
        ax.set_zlabel(' Altitude (km)')
#       ax.view_init(20., cb=fig.colorbar(pS, shrink=0.6)
        cb.set_label(Var['long_name'] + Var['units'])  # r'(m s$^{-1}$)')

        # Plot the horizontal dBZ contour field
        if plot_contour:
            if (con_field is not None):
                conVar = self._get_variable_dict(con_field)

                # Find the closest vertical point
                Zind = find_nearest_indices(self.height['data'][:], ppi_height)

                ax.contourf(Lon2D, Lat2D, conVar['data'][Zind, :, :], clevels,
                            vmin=vmin, vmax=vmax, cmap=cmap,
                            zdir='z')  # ,offset=surf_min)
            else:
                print("Need to set con_field and ppi_height")

        if plot_track:
            ax.plot(self.latitude['data'][:], self.longitude['data'][:],
                    self.height['data'][:] / 1000., zdir='z', c='k')

        return
