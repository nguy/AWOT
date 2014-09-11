"""
awot.display.p3tv
=========================

A group of scripts to create a display of NOAA P-3 tail Doppler radar. 

Created by Nick Guy.

UPCOMING MODIFICATIONS:
    Add the offset to the height plots

"""
# HISTORY::
#  19 Jun 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
# FUNCTIONS::
# DPJgrid_Horiz - Plot a horizontal field and optionally overlay track and winds
# plot_track_2d - Plot aircraft track on top of map
# plot_wind_2d - Plot wind vectors over map
# DPJgrid_3d - Plot a 3D field -STILL IN PROGRESS
#-------------------------------------------------------------------
# Load the needed packages
from pyart.io import read_cfradial as rcf
from pyart.graph import RadarDisplay_Airborne as getDisp
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize as cbNorm
from matplotlib.colorbar import ColorbarBase as cbBase
import numpy as np
from netCDF4 import chartostring

# Define various constants that may be used for calculations
C = 2.99792E8

# Hard code panel sizes and locations
ref_pan_axes = [0.04, 0.31, 0.51, 0.25]
vel_pan_axes = [0.04, 0.72, 0.51, 0.25]
LF_pan_axes = [0.61, 0.25, 0.36, 0.65]
    
cb_ref_axes = [0.07, 0.23, 0.25, 0.015]
cb_vel_axes = [0.07, 0.64, 0.25, 0.015]
cb_LF_axes = [0.60, 0.15, 0.25, 0.02]
#cb_ref_axes = [0.10, 0.25, 0.35, 0.27]
#cb_vel_axes = [0.10, 0.60, 0.35, 0.62]
#cb_LF_axes = [0.60, 0.10, 0.95, 0.12]
#===============================================================
class P3tv(object):
    """
    Overview
    --------
    To create a new P3tv instance:
    new_instance = P3tv() or new_instance = P3tv(filename)
    
    An object for creating plots from data in a radar object.

    Attributes
    ----------
    plots : list
        List of plots created.
    

    """

#################################################################

    def __init__(self, filename=None, verbose=False, raw=False):
    
        """
        If initialized with a filename (incl. path), will call
        cd2disp() to populate the class instance.
        If not, it simply instances the class but does not populate
        its attributes.
        verbose: Set to True for text output. Useful for debugging.
        
        """

        if isinstance(filename, str) != False:
            if raw:
                self.raw2disp(filename, verbose=verbose)
            else:
                self.cf2disp(filename, verbose=verbose)
        else:
            #Initializes class instance but leaves it to other methods to
            #populate the class attributes.
            return

#################################################################
    
    def cf2disp(self,fIn):
        """Read in the cfradial data file.
    
        INPUT PARAMETERS::
        fIn : str
            Input CFRadial format file (full path)
      
        OUTPUT::
        disp :
            PyArt radar class
        """
    # HISTORY::
    #  19 jun 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
    #--------------------------------------------------------
        rad = rcf(fIn) # Read the cfradial file
#        self.fields = rad.fields
        self.tdr_pyart_disp = getDisp(rad)
            
        # Adjust the Z components to make zero the ground
        Alt2D, rtmp = np.meshgrid(self.tdr_pyart_disp.ranges, self.tdr_pyart_disp.altitude/1000.)
        self.tdr_pyart_disp.z = self.tdr_pyart_disp.z + np.mean(self.tdr_pyart_disp.altitude)
    
        # Set a No Value where information does not exist
        nan = float('NaN')
    
        # Check for type of PRF
        if chartostring(rad.instrument_parameters['prt_mode']['data'])[0] == 'fixed':
           prf1 = 1. / rad.instrument_parameters['prt']['data'][0]
           nyq1 = rad.instrument_parameters['nyquist_velocity']['data'][0]
           prf2 = nan
           nyq2 = nan
        if chartostring(rad.instrument_parameters['prt_mode']['data'])[0] == 'dual':
           prf1 = nan#1. / rad.instrument_parameters['prt']['data'][0]
           nyq1 = nan#rad.instrument_parameters['nyquist_velocity']['data'][0]
           prf2 = nan
           nyq2 = nan
        if chartostring(rad.instrument_parameters['prt_mode']['data'])[0] == 'staggered':
           prf1 = 1. / rad.instrument_parameters['prt']['data'][0]
           nyq1 = rad.instrument_parameters['nyquist_velocity']['data'][0]
           prf2 = prf1 * rad.instrument_parameters['prt_ratio']['data'][0]
           nyq2 = prf2 * (C / rad.instrument_parameters['frequency']['data'][0]) / 4.
       

        # Grad values for display and put them in a dictionary
        self.dispVals = {'heading' : rad.heading['data'][0],
                         'drift' : rad.drift['data'][0],
                         'track' : rad.azimuth['data'][0],
                         'pitch' : rad.pitch['data'][0],
                         'roll' : rad.roll['data'][0],
                         'tilt' : rad.tilt['data'][0],
                         'AvgTilt' : np.mean(rad.tilt['data'][:]),
                         'GrndSpd' : nan,
                         'VertVel' : nan,
                         'Alt' : rad.altitude['data'][0],
                         'WndSpd' : nan,
                         'WndDir' : nan,
                         'nPulse' : rad.instrument_parameters['n_samples']['data'][0],
                         'NoiseTh' : nan,
                         'SQITh' : nan,
                         'NyqVel' : rad.instrument_parameters['nyquist_velocity']['data'][0],
                         'PRF1' : prf1,
                         'Nyq1' : nyq1,
                         'PRF2' : prf2,
                         'Nyq2' : nyq2,
                         'wavelength' : C / rad.instrument_parameters['frequency']['data'][0],
                         'platform' : rad.metadata['instrument_name'],
                         'lat' : rad.latitude['data'][0],
                         'lon' : rad.longitude['data'][0],
                         'title' : rad.metadata['title'],
                         'project' : rad.metadata['site_name']
                         }

#################################################################
    
    def lf2disp(self,fIn):
        """Read in the cfradial data file.
    
        INPUT PARAMETERS::
        fIn : str
            Input CFRadial format file (full path)
      
        OUTPUT::
        disp :
            PyArt radar class
        """
    # HISTORY::
    #  19 jun 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
    #--------------------------------------------------------
        rad = rcf(fIn) # Read the cfradial file
#        self.fields = rad.fields
        self.lf_pyart_disp = getDisp(rad)
    
#################################################################

    def add_tdr_sweep(self, field, fig=None, vmin=None, vmax=None, title=None,
                  RngMin=None, RngMax=None, HtMin=None, HtMax=None,
                  Gridlines=True, RngRings=True, mask_tuple=None, **kwargs):

        """Display a tail radar variable.
    
        INPUT PARAMETERS::
        field : str
            Name of field to display.
        
        OPTIONAL PARAMETERS::
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.
        vmin : float
            Luminance minimum value.
        vmax : float
            Luminance maximum value.
        RngMin : float
            Minimum horizontal range to display
        RngMax : float
            Maximum horizontal range to display
        HtMin : float
            Minimum altitude to display
        HtMax : float
            Maximum altitude to display
        title : str
            Title to label plot with, None to use default title generated from
            the field and sweep parameters. Parameter is ignored if title_flag
            is False.
        title_flag : bool
            True to add a title to the plot, False does not add a title.
        GridLines : Boolean
            True to turn on gridlines (default), False to turn them off
        RngRings : Boolean
            True to turn on range rings (default), False to turn off
        mask_tuple : (str, float)
            Tuple containing the field name and value below which to mask
            field prior to plotting, for example to mask all data where
            NCP < 0.5 set mask_tuple to ['NCP', 0.5]. None performs no masking.

        OUTPUT::
        p               = Plot
        """
    # HISTORY::
    #  19 Jun 2014 - Nick Guy NOAA/NSSL/WRDD, NRC 
    #-------------------------------------------------------
        if field == 'DBZ':
            ax = fig.add_axes(ref_pan_axes)
            cmap = 'gist_ncar'
        elif field == 'VEL':
            ax = fig.add_axes(vel_pan_axes)
            cmap = 'RdBu_r'
        xlim, ylim = [RngMin, RngMax], [HtMin, HtMax]
        
        # Grab the figure object
        fig = self._parse_fig(fig)
        
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
        vmin, vmax = self._parse_vmin_vmax(field, vmin, vmax, 'tdr')
                
        if title == None:
            titleTx = self.generate_title(radStr='tdr')
        else:
            titleTx = title
            
        # Plot the data
        self.tdr_pyart_disp.plot_sweep_grid(field,vmin=vmin,vmax=vmax,ax=ax,fig=fig,cmap=cmap,
                       colorbar_flag=False,mask_tuple=mask_tuple,title=titleTx)
                   
        self.tdr_pyart_disp.set_limits(xlim=xlim,ylim=ylim) # Set the axes ranges
        
        # Add gridlines and range rings if requested
        if Gridlines:
            self.tdr_pyart_disp.plot_grid_lines(ls=':',lw=2) # Enable grid lines
        if RngRings:
            self.tdr_pyart_disp.plot_range_rings([15.,30.,45.,60.],ax=ax,col='gray',ls='-',lw=1)
    
#################################################################

    def add_lf_ppi(self, field, fig=None, vmin=None, vmax=None, title=None, 
                   Gridlines=True, RngRings=True, mask_tuple=None, **kwargs):

        """Add a PPI plot of the lower fuselage radar data.
    
        INPUT PARAMETERS::
        field : str
            Name of field to display.
        
        OPTIONAL PARAMETERS::
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.
        vmin : float
            Luminance minimum value.
        vmax : float
            Luminance maximum value.
        RngMin : float
            Minimum horizontal range to display
        RngMax : float
            Maximum horizontal range to display
        """
    # HISTORY::
    #  24 Jun 2014 - Nick Guy NOAA/NSSL/WRDD, NRC 
    #--------------------------------------------------------
        if field == 'DBZ':
            cmap = 'gist_ncar'
        elif field == 'VEL':
            cmap = 'RdBu_r'
    
        # Grab the figure object
        fig = self._parse_fig(fig)
        
        # Set axis object
        ax = fig.add_axes(LF_pan_axes)
        
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
        vmin, vmax = self._parse_vmin_vmax(field, vmin, vmax, 'lf')
                
        if title == None:
            titleTx = self.generate_title(radStr='lf')
        else:
            titleTx = title
            
        # Plot the data
        self.lf_pyart_disp.plot_ppi(field, vmin=vmin,vmax=vmax,ax=ax,fig=fig,cmap=cmap,
                       colorbar_flag=False,mask_tuple=mask_tuple,title=titleTx)
        
        # Add gridlines and range rings if requested
        if Gridlines:
            self.tdr_pyart_disp.plot_grid_lines() # Enable grid lines
        if RngRings:
            rngs1 = [50.,150.,250.,350.]
            rngs2 = [100.,200.,300.,400.]
            self.tdr_pyart_disp.plot_range_rings(rngs1,ax=ax,col='gray',ls=':')
            self.tdr_pyart_disp.plot_range_rings(rngs2,ax=ax,col='gray',ls='-', lw=1)
                       
#################################################################

    def add_colorbar(self, field, fig=None, vmin=None, vmax=None, mappable=None,
                     radStr=None):

        """Add a colorbar to existing plot.
    
        INPUT PARAMETERS::
        ax : Axis
            Axis to plot on.
        """
    # HISTORY::
    #  19 Jun 2014 - Nick Guy NOAA/NSSL/WRDD, NRC 
    #--------------------------------------------------------
        # parse parameters
        fig = self._parse_fig(fig)
        
        if radStr is None:
            radStr = 'tdr'
        elif radStr == 'tdr':
            cb_orient = 'horizontal'
            if field == 'DBZ':
                ax = fig.add_axes(cb_ref_axes)
                cmap = 'gist_ncar'            
                label = 'Reflectivity (dBZ)'
                labx = (cb_ref_axes[0] + cb_ref_axes[2]) + 0.08
                laby = cb_ref_axes[1] + (cb_ref_axes[3]/2.)
            elif field == 'VEL':
                ax = fig.add_axes(cb_vel_axes)
                cmap = 'RdBu_r'
                label = r' Velocity (m s$^{-1}$)'
                labx = (cb_vel_axes[0] + cb_vel_axes[2]) + 0.08
                laby = cb_vel_axes[1] + (cb_vel_axes[3]/2.)
        elif radStr == 'lf':
            cb_orient = 'horizontal'
            ax = fig.add_axes(cb_LF_axes)
            labx = (cb_LF_axes[0] + cb_LF_axes[2]) + 0.07
            laby = cb_LF_axes[1] + (cb_LF_axes[3]/2.)
            if field == 'DBZ':
                cmap = 'gist_ncar'            
                label = 'Reflectivity (dBZ)'

        if mappable is None:
            mappable = self.tdr_pyart_disp.plots[-1]
            
        # Map the colors to match min/max values
        norm = cbNorm(vmin=vmin, vmax=vmax)
        # Create the colorbar
#        cbar = cbBase(ax,cmap=cmap,norm=norm,orientation='horizontal')
        cbar = fig.colorbar(mappable=mappable,cax=ax,cmap=cmap,norm=norm,orientation=cb_orient)
        # Add a label
        fig.text(labx,laby,label,horizontalalignment='center',
        verticalalignment='top',transform=ax.transAxes)

#################################################################

    def add_parameters(self, fig=None):
        """Place parameters on the output display
    
    
        """
    # HISTORY::
    #  23 jun 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
#--------------------------------------------------------
        # parse parameters
        fig = self._parse_fig(fig)
        ProjStr = (('Proj: %s           Platform: %s') % 
                 (self.dispVals['project'], self.dispVals['platform']))
        col1 = (('Lat: %.2f\n' + 'Lon: %.2f\n' + 'Alt: %.1f\n' + 'GrndSpd: %.1f\n' ) %
                 (self.dispVals['lat'], self.dispVals['lon'],
                  self.dispVals['Alt'], self.dispVals['GrndSpd']))
        col2 = (('Heading: %.1f\n' + 'Drift: %.1f\n' + 'Track: %.1f\n' +
                'Pitch: %.1f\n') % 
                (self.dispVals['heading'], self.dispVals['drift'], self.dispVals['track'],
                 self.dispVals['pitch']))
        col3 = (('Roll: %.1f\n' + 'Tilt: %.1f\n' + 'AveTilt: %.1f\n' +
                'VertVel: %.1f\n') % 
                (self.dispVals['roll'], self.dispVals['tilt'], self.dispVals['AvgTilt'],
                 self.dispVals['VertVel']))
        col4 = (('WndSpd: %.1f\n' + 'WndDir: %.1f\n' + '#Pulses: %i\n' +
                'Noise Th: %.1f\n' +'SQI Th: %.1f') % 
                (self.dispVals['WndSpd'], self.dispVals['WndDir'], self.dispVals['nPulse'],
                 self.dispVals['NoiseTh'], self.dispVals['SQITh']))
        col5 = (('UnambVel: %.1f\n' + 'PRF1: %.1f\n' + 'Nyq1: %.1f\n' +
                'PRF2: %.1f\n' +'Nyq2: %.1f') % 
                (self.dispVals['NyqVel'], self.dispVals['PRF1'], self.dispVals['Nyq1'],
                 self.dispVals['PRF2'], self.dispVals['Nyq2']))
                 
        fig.text(0.57,0.05,ProjStr,horizontalalignment='left',verticalalignment='top',fontsize=18)
        fig.text(0.02,0.17,col1,horizontalalignment='left',verticalalignment='top',linespacing=1.75)
        fig.text(0.12,0.17,col2,horizontalalignment='left',verticalalignment='top',linespacing=1.75)
        fig.text(0.22,0.17,col3,horizontalalignment='left',verticalalignment='top',linespacing=1.75)
        fig.text(0.32,0.17,col4,horizontalalignment='left',verticalalignment='top',linespacing=1.75)
        fig.text(0.42,0.17,col5,horizontalalignment='left',verticalalignment='top',linespacing=1.75)
#################################################################
##### PARSEING METHODS ######
#################################################################

    def _parse_ax(self, ax):
        """ Parse and return ax parameter. """
        if ax is None:
            ax = plt.gca()
        return ax

    def _parse_fig(self, fig):
        """ Parse and return ax and fig parameters. """
        if fig is None:
            fig = plt.gcf()
        return fig

    def _parse_ax_fig(self, ax, fig):
        """ Parse and return ax and fig parameters. """
        if ax is None:
            ax = plt.gca()
        if fig is None:
            fig = plt.gcf()
        return ax, fig

    def _parse_vmin_vmax(self, field, vmin, vmax, radStr):
        """ Parse and return vmin and vmax parameters. """
        if radStr == 'tdr':
            field_dict = self.tdr_pyart_disp.fields[field]
        elif radStr == 'lf':
            field_dict = self.lf_pyart_disp.fields[field]
        if vmin is None:
            if 'valid_min' in field_dict:
                vmin = field_dict['valid_min']
            else:
                vmin = -6   # default value
        if vmax is None:
            if 'valid_max' in field_dict:
                vmax = field_dict['valid_max']
            else:
                vmax = 100
        return vmin, vmax

#################################################################
##### GET METHODS ######
#################################################################

    def _get_data(self, field, mask_tuple):
        """ Retrieve and return data from a plot function. """
        data = self.fields[field]['data'][:]

        if mask_tuple is not None:  # mask data if mask_tuple provided
            mask_field, mask_value = mask_tuple
            mdata = self.fields[mask_field]['data'][:]
            data = np.ma.masked_where(mdata < mask_value, data)
        return data
        
#################################################################
##### PLOT ADJUST METHODS ######
#################################################################

    def generate_title(self,radStr=None):
        """ Generate the figure title using a default title. """
            
        if radStr is None:
            radStr = 'tdr'
        elif radStr == 'tdr':
            time_str = self.tdr_pyart_disp.time_begin.isoformat() + 'Z'
            
            if self.tdr_pyart_disp.fixed_angle <= -0.5:
                AntAng = 'Aft'
            elif self.tdr_pyart_disp.fixed_angle >= 0.5:
                AntAng = 'For'
            else:
                AntAng = 'Nadir'
            
            title = "%s %s %s " % (self.tdr_pyart_disp.radar_name, AntAng, time_str)
        elif radStr == 'lf':
            time_str = self.lf_pyart_disp.time_begin.isoformat() + 'Z'
            
            title = "%s %s " % (self.lf_pyart_disp.radar_name, time_str)
            
        return title

    def generate_filename(self, ext='png'):
        """
        Generate a filename for a plot.

        Generated filename has form:
            radar_name_field_time.ext

        Parameters
        ----------
        field : str
            Field plotted.
        ext : str
            Filename extension.

        Returns
        -------
        filename : str
            Filename suitable for saving a plot.

        """
        name_s = self.tdr_pyart_disp.radar_name.replace(' ', '_')
        time_s = self.tdr_pyart_disp.time_begin.strftime('%Y%m%d%H%M%S')
        time_s = time_s.replace('-',':')
        return '%s_%s_TDR_LF.%s' % (name_s, time_s, ext)