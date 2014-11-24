# -*- coding: utf-8 -*-
"""
Airborne.py - Class for airborne observations
"""

from .graph.common import create_basemap_instance

from awot.io.FlightDataFile import read_flight
from awot.io.RadarDataFile import read_radar
from awot.io.write_radar_netcdf import radar2nc

import inspect
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt

########################
## BEGIN MAIN CODE
########################

class AirborneData(object):

    '''
    AirborneData class to hold Airborne data dictionaries. 
    Can be returned from a reader object (e.g. FlightDataFile or RadarDataFile) 
    
    e.g.
    Data = AirborneData(reader=FlightDataFile or RadarDataFile)
    
    or 
    
    created with no data and filled through separate reader functions
    
    e.g.
    Data = AirborneData()
    Data.get_flight_data(fname=long_filename)
    
    '''

    def __init__(self, reader=None, data_type=None):

        """
        Initialize the class if a reader class is provided.
        
        Parameters:
        ------------
        reader : Either a RadarData or FlightData class
        
        If nothing is specified here then the class is instantiated.
        The user can then call the get_*_data to get pull in data.
        This is the preferred way at this time.
        """
            
        if reader is None:
           return
        else:
            if (inspect.isclass(type(reader)) or inspect.isclass(reader)):
                AirborneData = deepcopy(reader)

    ##########################################
    # Get data from files
    ##########################################
    
    def get_flight_data(self, fname=None, platform='', file_format='netcdf',
                       instrument=None):
        '''
        Return a FlightDataFile reader instance
        
        Parameters::
        ------------
        fname : str
            Long path filename
        platform : str
            Platform name (see io.FlightDataFile)
        file_format : str
            Format of file to read
        '''
        if fname is None:
            print "Must supply input file!!"
            return
        else:
            FlightData = read_flight(filename=fname, platform=platform, 
                                     file_format=file_format, instrument=instrument)
            
        # Now need to add dictionary to AirborneData class
            self.flight_data = FlightData.flight_data
    
    #####################
            
    def get_radar_data(self, fname=None, platform='', file_format='netcdf',
                       instrument=None):
        '''
        Return a RadarDataFile reader instance
        
        Parameters::
        ------------
        fname : str
            Long path filename
        platform : str
            Platform name (see io.RadarDataFile)
        file_format : str
            Format of file to read
        instrument : str
            Type (name) of instrument to process
            Currently the following arguments are valid:
                'tdr_grid' - Tail Doppler radar gridded files (e.g. dual-Doppler analysis)
                'tdr_sweep' = Tail Doppler radar (Native coordinate data)
                'lf'  - Lower Fuselage radar
                'ground' - A ground-based radar system, read in using PyArt
        '''
        if fname is None:
            print "Must supply input file!!"
            return
        else:
            RadarData = read_radar(filename=fname, platform=platform, 
                                   file_format=file_format, instrument=instrument)
            
        # Now need to add dictionary to AirborneData class
            
            if instrument is None:
               print "Need to supply instrument type"
               return
            
            # Add the dictionary to existing class
            if (platform.upper() == 'P3') or (platform.upper() == 'P-3') or \
            (platform.upper() == 'NOAA_P3') or \
            (platform.upper() == 'NOAA P3') or \
            (platform.upper() == 'ELDORA'):   
                if instrument.lower() == 'tdr_grid':
                    self.tdr_radar_data = RadarData.radar_data
                elif instrument.lower() == 'lf':
                    self.lf_radar_data = RadarData.radar_data
                elif instrument.lower() == 'tdr_sweep':
            	    self.tdr_sweep_radar_data = RadarData.radar_data
            
            elif platform.upper() == 'FALCON':
                if instrument.lower() == 'radar':
                    self.rasta_radar_data = RadarData.radar_data
                elif instrument.lower() == 'microphysics':
                    self.rasta_microphysics_data = RadarData.radar_data
                else:
                    print "Set instrument parameter to 'radar' or 'microphysics'"
                    return
            
            elif (platform.upper() == 'GROUND'):
                self.ground_radar_data = RadarData.radar_data
            else:
                print "Could not grab data, check the instrument type!"
    
    ##############################################
    # Establish a basemap instance for plotting
    ##############################################
    
    def create_basemap(self,corners=None, proj=None, resolution='l', area_thresh=1000.,
                      track_lw=1.5, alpha=1.,
                      meridians=True, parallels=True, dLon=2., dLat=2.,
                      coastlines=True, countries=True, states=False, counties=False, 
                      rivers=False, etopo=False, ax=None):
        '''Instantiate a basemap instance'''
                      
        # Create a basemap instance            
        bm = create_basemap_instance(corners=corners, proj=proj, 
                   resolution=resolution, area_thresh=area_thresh,
                   meridians=meridians, parallels=parallels, dLon=dLon, dLat=dLat,
                   coastlines=coastlines, countries=countries, states=states, 
                   counties=counties, rivers=rivers, etopo=etopo, ax=ax)
                   
        # Save the basemap instance for further plotting    
        self.basemap = bm
        
    ##########################
    # Save methods
    ##########################
    
    def save_figure(self, figName='awot_plot', figType='png', **kwargs):
        '''Save the current plot
        
        Parameters::
        ------------
        figName : str
            Figure name
        figType : str
            Figure format, default to .png
        
        '''
        plt.gca()
        plt.gcf()
        plt.savefig(figName+'.'+figType, format=figType)
        print "Saved figure: " + figName+'.'+figType
        
        # Now close the plot to make sure matplotlib is happy
        plt.close()
        
    def write_radar_netcdf(self, radar=None, Outfile=None):
        '''Save a radar instance as a NetCDF output file'''
        if radar is None:
            print "Must specify the radar instance to process"
        else:
            radar2nc(radar, Outfile=Outfile)

    ##########################################
        
    def get_dropsonde_data(self, fname=None):
        '''
        Return DropsondeDataFile reader instance
        
        Parameters::
        ------------
        fname : str
            Long path filename
            
        TODO: Get working
        '''
        if fname is None:
            print "Must supply input file!!"
            return
        else:
            pass

    ##########################################
        
    def get_lidar_data(self, fname=None):
        '''
        Return LidarDataFile reader instance
        
        Parameters::
        ------------
        fname : str
            Long path filename
            
        TODO: Get working.  Import another package?
        '''
        if fname is None:
            print "Must supply input file!!"
            return
        else:
            pass

    ##########################################
        
    def get_particle_probe(self, fname=None):
        '''
        Return ProbeDataFile reader instance
        
        Parameters::
        ------------
        fname : str
            Long path filename
            
        TODO: Get working.  Import from pyparticleprobe package?
        '''
        if fname is None:
            print "Must supply input file!!"
            return
        else:
            pass
        
    ##########################################
    


    ##########################################
    


    ##########################################
    


    ##########################################
    
    

    ##########################################
    