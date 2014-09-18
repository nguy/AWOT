# -*- coding: utf-8 -*-
"""
FlightDataFile.py - Class for reading Data from file
"""
import os

#from ..AirborneData import AirborneData
from ..io.read_p3_flight import flight_data as p3_read_flight
from ..io.read_citation_flight import flight_data as citation_read_flight

import numpy as np
########################
## BEGIN MAIN CODE
########################

def read_flight(filename=None, platform='', file_format='netcdf', \
                initialize=False):
    '''
    Takes a filename pointing to a aircraft flight level data file
     and returns a FlightLevel object.
        
    Parameters::
    ------------
    filename : str
        Long path filename of data file.
    platform : str
        Platform for processing, see FileReader.
    file_format : str
        Format of input file, see FileReader.
    
    '''
    if filename is None:
        print "Need to provide a filename!"
        return
        
    reader = FileReader(filename=filename, platform=platform, file_format=file_format)
    
    if initialize:
        FlightData = AirborneData(reader)
    
        return FlightData
    else:
        return reader
    
###################################################

class FileReader(object):

    '''
    FileReader class to process data files.  
    '''
    
    def __init__(self, filename=None, platform=None, file_format=None):

        """
        If initialized with a filename (incl. path), will call
        ***_read_flight() to populate the class instance.
        If not, it simply instances the class but does not populate
        its attributes.
        verbose: Set to True for text output. Useful for debugging.
        
        Parameters::
        ------------
        filename : str
            Filename (including path) of file to process
        platform : str
            Name of aicraft 
            Currently supported: 
                'p3' or p-3' (NOAA WP-3D)
                'citation' (Univ North Dakota Citation)
        file_format : str
            Format of input file.  Each platform currently has
            a specific file type.  This may be extended in the future.
            Currently supported:
                'netcdf'
                'ascii'
        
        """
        
        if isinstance(filename, str) != False:
        
            # Perform checks to make sure we know what platform to look up
            if platform is None :
                print "Need to specify the platform for reader"
                return
            
            if ((platform.upper() == 'P3') or (platform.upper() == 'P-3')):
                if file_format.upper() == 'NETCDF':
                    flight = p3_read_flight(filename)
            elif (platform.upper() == 'CITATION'):
                if ((file_format.upper() == 'ASCII') or \
                (file_format.upper() == 'NASA AMES') or \
                (file_format.upper() == 'NA')):
                    flight = citation_read_flight(filename)
                elif file_format.upper() == 'NETCDF':
                    print "No NetCDF reader currently exists for citation data"
                    return
                else:
                    print "Check the format call!"
                    return
                                              
            # Calculate meridional and zonal wind components
            Uwind, Vwind = self._winduv(flight)
            
            # Add to the dictionary
            flight['Uwind'] =  Uwind
            flight['Vwind'] =  Vwind
            
            # Record the data into the variable 
            self.flight_data = flight
        else:
            #Initializes class instance but leaves it to other methods to
            #populate the class attributes.
            return
            
    def _winduv(self,flight):
        """
        Calculate the horizontal windcomponents (u,v) from wind angle and speed
            U wind is positive blowing towards east
            V wind is positive blowing towards north
        """
        U = -np.cos(np.radians(flight['wind_dir'])) * flight['wind_spd']
        V = -np.sin(np.radians(flight['wind_dir'])) * flight['wind_spd']
        
        return U, V

            
        
