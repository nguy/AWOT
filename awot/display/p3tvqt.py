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
import numpy as np
from netCDF4 import chartostring

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.console
from pyqtgraph.dockarea import *

#from matplotlib.mlab import griddata as griddata2
from scipy.interpolate import griddata
from scipy.ndimage.interpolation import rotate

# Define various constants that may be used for calculations
C = 2.99792E8
#===============================================================
class P3tvQt(object):
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

    def create_display(self,winsize=None):
        
        if winsize is None:
           winX = 1000
           winY = 500
        else:
           winX = winsize[0]
           winY = winsize[1]
           
        # Create the application
        app = QtGui.QApplication([])
        win = QtGui.QMainWindow()
        self.area = DockArea()

        # Set some window parameters
        win.setCentralWidget(self.area)
        win.resize(winX,winY)
        win.setWindowTitle('Py3TV Display')

        ## Create docks, place them into the window one at a time.
        ## Note that size arguments are only a suggestion; docks will still have to
        ## fill the entire dock area and obey the limits of their internal widgets.
        self.dock1 = Dock("TDR Aft VR", size=(winX/1.3, winY*(2./5.)))
        self.dock2 = Dock("TDR Fore VR", size=(winX/1.3, winY*(2./5.)))
        self.dock3 = Dock("TDR Aft dBZ", size=(winX/1.3, winY*(2./5.)))
        self.dock4 = Dock("TDR Fore dBZ", size=(winX/1.3, winY*(2./5.)))
        self.dock5 = Dock("Lower Fuselage", size=(winX/1.7, winY*(3.6/5.)))
        self.dock6 = Dock("Parameters ", size=(winX/1.7, winY*(1.4/5.)))     ## give this dock the minimum possible size

        self.area.addDock(self.dock1, 'left')        ## place d1 at left edge of dock area (it will fill the whole space since there are no other docks yet)
        self.area.addDock(self.dock2, 'above', self.dock1)   ## stack d2 on top of d1
        self.area.addDock(self.dock3, 'bottom', self.dock1)  ## place d3 at bottom edge of d1
        self.area.addDock(self.dock4, 'above', self.dock3)   ## stack d4 on top of d3
        self.area.addDock(self.dock5, 'right')       ## place d5 at right edge of dock area
        self.area.addDock(self.dock6, 'top', self.dock5)     ## place d6 at top edge of d5
        
        return app, win

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

    def add_tdr_sweep(self, field, dock, vmin=None, vmax=None, title=None,
                  RngMin=None, RngMax=None, HtMin=None, HtMax=None,
                  Gridlines=True, RngRings=True, mask_tuple=None, **kwargs):

        """Display a tail radar variable.
    
        INPUT PARAMETERS::
        field : str
            Name of field to display.
        dock : Dock
            Which dock to add this to.
        
        OPTIONAL PARAMETERS::
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
            cmap = 'gist_ncar'
        elif field == 'VEL':
            cmap = 'RdBu_r'
                
        # Create title if desired
        if title == None:
            titleTx = self.generate_title(radStr='tdr')
        else:
            titleTx = title
        
        # Grab the plot limits
        xlim, ylim = [RngMin, RngMax], [HtMin, HtMax]
        
        # Grab the data
        data = self._get_data(field, mask_tuple, 'tdr')
        x, z = self._get_x_z(field, 'tdr')

        # parse parameters
        vmin, vmax = self._parse_vmin_vmax(field, vmin, vmax, 'tdr')

        # Grid the data for display
#        xi, zi = self._get_interp_x_z(x, z, .150, .150)
        xi, zi = self._get_interp_x_z(x, z, 150, 150)
#        datai = griddata2(x.flatten(),z.flatten(),data.flatten(),xi,zi)
        datai = self._get_interp_data(x, z, data, xi, zi, method='nearest')
        
        # Calculate scale for image
        xscale, zscale = self._get_scaling(xi, zi, datai)
#        xoffset, zoffset = self._get_offsets(xi, zi)
        print HtMin/zscale,HtMax/zscale
        # Plot the data
        pbox = pg.ViewBox()
#        pbox.setRange(xRange=[RngMin,RngMax],yRange=[HtMin,HtMax])
#        pbox.setXRange(RngMin,RngMax)
#        pbox.setYRange((HtMin/zscale)+zi.min(),(HtMax/zscale)+zi.min())#,update=False,padding=None)
        pbox.setYRange(HtMin+zi.min(),(len(zi)*zscale/HtMax)+zi.min())#,update=False,padding=None)      
        
        # Create a plot item to establish axes and gridlines, set title label
        plt = (pg.PlotItem(title=titleTx, viewBox=pbox,
              labels={'bottom': ('Cross-track distance', 'm'),'left': ('Altitude', 'm')}))
#        plt.setXRange(RngMin,RngMax)
#        plt.setYRange(HtMin,HtMax,update=False,padding=None)
        
        # Add gridlines if requested
        if Gridlines:
            plt.showGrid(x=True, y=True)
            
        # Create image item to associate with ImageView
        w = pg.ImageView(view=plt)
        w.setImage(datai.T,levels=(vmin,vmax),scale=[xscale,zscale],#[::-1,:]
                  axes={'x':1,'y':0},pos=[xi.min(),zi.min()])#[xoffset,zoffset])
        w.ui.roiBtn.hide() # Hide the ROI button on display
        w.ui.normBtn.hide() # Hide the Norm button on display
        w.view.invertY(False)
        w.view.invertX(False)
#        w.view.setYRange(HtMin/zscale, HtMax/zscale)
#        plt.setAspectLocked(True)
        
        # Add the plot to the appropriate dock
        dock.addWidget(w)

#        self.tdr_pyart_disp.set_limits(xlim=xlim,ylim=ylim) # Set the axes ranges
    
#################################################################

    def add_lf_ppi(self, field, dock, vmin=None, vmax=None, title=None, 
                   Gridlines=True, RngRings=True, mask_tuple=None, **kwargs):

        """Add a PPI plot of the lower fuselage radar data.
    
        INPUT PARAMETERS::
        field : str
            Name of field to display.
        dock : Dock
            Which dock to add this to.
        
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
                
        # Create title if desired
        if title == None:
            titleTx = self.generate_title(radStr='lf')
        else:
            titleTx = title
            
        # Set the title within plot space
#        w = pg.PlotWidget(title=titleTx)
        
        # Grab the data
        data = self._get_data(field, mask_tuple, 'lf')
        x, y = self._get_x_y(field, 'lf')
        
        # parse parameters
        vmin, vmax = self._parse_vmin_vmax(field, vmin, vmax, 'lf')
        
        # Grid the data for display
        xi, yi = self._get_interp_x_y(x, y, 1.0, 1.0)
        xi, yi = self._get_interp_x_y(x, y, 1000, 1000)
        datai = self._get_interp_data(x, y, data, xi, yi, method='nearest')
        
        # Calculate scale for image
        xscale, yscale = self._get_scaling(xi, yi, datai)
#        xoffset, yoffset = self._get_offsets(xi, yi)

        # Create a color table
        Lut = pg.GradientEditorItem()
        Lut.loadPreset('spectrum')
#        Lut.getLookupTable(256)
#        argb = pg.makeARGB(datai, lut= Lut, levels=(vmin,vmax))
        stops=np.r_[-1.0,-0.5,0.5,1.0]
        colors=np.array([[0,0,1,0.7],[0,1,0,0.2],[0,0,0,0.8],[1,0,0,1.0]])
        cm=pg.ColorMap(stops,colors)
        vals = cm.map(datai)#np.random.rand(512, 512)*2 - 1
#        print datai.shape,vals.shape
        
        # Create a plot item to establish axes and gridlines, set title label
        plt = (pg.PlotItem(title=titleTx, 
              labels={'bottom': ('E-W distance', 'm'), 'left': ('N-S Distance', 'm')}))
        
        # Add gridlines if requested
        if Gridlines:
            plt.showGrid(x=True, y=True)
        # Add range rings if requested
#        if RngRings:
#            rngs1 = [50.,150.,250.,350.]
#            rngs2 = [100.,200.,300.,400.]
#            self.plot_range_rings(plt, rngs1)#,,col='gray',ls=':')
#            self.plot_range_rings(plt, rngs2)#,col='gray',ls='-', lw=1)
                    
        # Create an image item to add color to image
        pltIm = pg.ImageItem()#datai
#        pltIm.setImage(datai.T,levels=(vmin,vmax))
#        pltIm.setImage(vals,levels=(vmin,vmax))
#        pltIm.setLookupTable(Lut)
#        pltIm.setData(vals)
                
        # Create the plot item and add the image
        # The axes call is not needed, but allows clarification of axes plotting
        w = pg.ImageView(view=plt, imageItem=pltIm)
        w.setImage(datai.T,levels=(vmin,vmax), pos=[xi.min(), yi.min()], #[::-1,:]
                  scale=[xscale,yscale],axes={'x':1,'y':0})
        w.ui.roiBtn.hide() # Hide the ROI button on display
        w.ui.normBtn.hide() # Hide the Norm button on display
        w.view.invertY(False)
        w.view.invertX(False)
#        plt.setAspectLocked(False)
        
        # Add the plot to the appropriate dock
        dock.addWidget(w)
        
        
#        if Rng
                       
#################################################################

    def updateLUT(ui,LUT=None):
        dtype = ui.dtypeCombo.currentText()
        if dtype == 'uint8':
            n = 256
        else:
            n = 4096
        LUT = ui.gradient.getLookupTable(n, alpha=ui.alphaCheck.isChecked())
        return LUT

#################################################################

    def add_parameters(self, dock):
        """Place parameters on the output display
    
    
        """
    # HISTORY::
    #  23 jun 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
#--------------------------------------------------------
        # Create formated text strings in HTML
        def html_out(Str,Value):
            html = '<div style="text-align: center"><span style="color: #FFFFFF;">'+Str+': </span><span style="color: #FF0000; font-size: 16pt;">'+Value+'</span></div>'
            return html
            
        # parse parameters
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
                 
#        fig.text(0.57,0.05,ProjStr,horizontalalignment='left',verticalalignment='top',fontsize=18)
#        fig.text(0.02,0.17,col1,horizontalalignment='left',verticalalignment='top',linespacing=1.75)
#        fig.text(0.12,0.17,col2,horizontalalignment='left',verticalalignment='top',linespacing=1.75)
#        fig.text(0.22,0.17,col3,horizontalalignment='left',verticalalignment='top',linespacing=1.75)
#        fig.text(0.32,0.17,col4,horizontalalignment='left',verticalalignment='top',linespacing=1.75)
#        fig.text(0.42,0.17,col5,horizontalalignment='left',verticalalignment='top',linespacing=1.75)
        
        ## Create text object, use HTML tags to specify color/size
        #, anchor=(-0.3,1.3))        

        plt = pg.PlotView()
#        plt.hideAxis('bottom')
#        plt.hideAxis('left')
#        plt.view.hideButtons()
        plt.ui.roiBtn.hide()
        plt.ui.normBtn.hide()
        
        tx1 = pg.TextItem(html=html_out('Proj',self.dispVals['project']))
        tx2 = pg.TextItem(html=html_out('Platform',self.dispVals['platform']))
        tx1.setPos(-30.,-15)
        tx2.setPos(-10.,-15)
        
        
        plt.addItem(tx1)
        w = pg.ImageView(view=plt)
        
        dock.addWidget(w)
        
#        TblInfo = np.array([('Lat',self.dispVals['lat'],'Heading',self.dispVals['heading'],'Roll',self)],
#                        dtype=[('1 ',object),('2 ',object),('3 ',object),('4 ',object),('5 ',object)])
#        dock.setData(data)


        
        # Add the plot to the appropriate dock
#        dock.addWidget(w)
        
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

    def _get_data(self, field, mask_tuple, radStr):
        """ Retrieve and return data from a plot function. """
        if radStr == 'tdr':
            data = self.tdr_pyart_disp.fields[field]['data'][:]
        elif radStr == 'lf':
            data = self.lf_pyart_disp.fields[field]['data'][:]

        if mask_tuple is not None:  # mask data if mask_tuple provided
            mask_field, mask_value = mask_tuple
            mdata = self.fields[mask_field]['data'][:]
            data = np.ma.masked_where(mdata < mask_value, data)
        return data
        
    def _get_x_z(self, field, radStr):
        """ Retrieve and return x and y coordinate in km. """
#        if radStr == 'tdr':
#            return self.tdr_pyart_disp.x/1000., self.tdr_pyart_disp.z/1000.
#        elif radStr == 'lf':
#            return self.lf_pyart_disp.x/1000., self.lf_pyart_disp.z/1000.
        if radStr == 'tdr':
            return self.tdr_pyart_disp.x, self.tdr_pyart_disp.z
        elif radStr == 'lf':
            return self.lf_pyart_disp.x, self.lf_pyart_disp.z
        
    def _get_x_y(self, field, radStr):
        """ Retrieve and return x and y coordinate in km. """
#        if radStr == 'tdr':
#            return self.tdr_pyart_disp.x/1000., self.tdr_pyart_disp.y/1000.
#        elif radStr == 'lf':
#            return self.lf_pyart_disp.x/1000., self.lf_pyart_disp.y/1000.
        if radStr == 'tdr':
            return self.tdr_pyart_disp.x, self.tdr_pyart_disp.y
        elif radStr == 'lf':
            return self.lf_pyart_disp.x, self.lf_pyart_disp.y
            
    def _get_interp_x_y(self, x, y, Xstep, Ystep):
        return np.arange(x.min(), x.max(), Xstep), np.arange(y.min(), y.max(), Ystep)
            
    def _get_interp_x_z(self, x, z, Xstep, Zstep):
        return np.arange(x.min(), x.max(), Xstep), np.arange(z.min(), z.max(), Zstep)
        
    def _get_interp_data(self, v1, v2, data, v1i, v2i, method=None):
        if method is None:
            method = 'nearest'
        else:
            method = method
        V1, V2 = np.meshgrid(v1i, v2i)
        return griddata((v1.flatten(), v2.flatten()), data.flatten(), (V1,V2), method=method)
        
    def _get_scaling(self, v1, v2, data):
        return (v1.max() - v1.min()) / data.shape[0], (v2.max() - v2.min()) / data.shape[1]
    
    def _get_offsets(self, v1, v2):
        return v1.min() + (v1[len(v1)/2] - np.median(v1)), v2.min() + (v2[len(v2)/2] - np.median(v2))
        
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
        
    def plot_range_rings(self, plt, range_rings, col=None, ls=None, lw=None):
        """
        Plot a series of range rings.

        Parameters
        ----------
        range_rings : list
            List of locations in km to draw range rings.
        col : str or value
            Color to use for range rings.
        ls : str
            Linestyle to use for range rings.

        """
        for range_ring_location_m in range_rings:
            self.plot_range_ring(plt, range_ring_location_m, col=col, ls=ls)

    def plot_range_ring(self, plt, range_ring_location_m, npts=100,
                        col=None, ls=None, lw=None):
        """
        Plot a single range ring.

        Parameters
        ----------
        plt : PlotItem
            PlotItem occurrence to create plot.
        range_ring_location_m : float
            Location of range ring in m.
        npts: int
            Number of points in the ring, higher for better resolution.
        col : str or value
            Color to use for range rings.
        ls : str
            Linestyle to use for range rings.

        """
        theta = np.linspace(0, 2 * np.pi, npts)
        r = np.ones([npts], dtype=np.float32) * range_ring_location_m
        x = r * np.sin(theta)
        y = r * np.cos(theta)
        if lw is None:
            lw = 2
        if ls is None:
            ls = '-'
        if col is None:
            col = 'k'
        plt.plot(x,y)
        