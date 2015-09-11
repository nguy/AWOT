# -*- coding: utf-8 -*-
"""
FlightDataFile.py - Class for reading Data from file
"""
import os
import numpy as np

from ..io.read_p3_flight import flight_data as p3_read_flight
from ..io.read_citation_flight import flight_data as citation_read_flight
from ..io.read_latmos_falcon_flight import flight_data as latmos_falcon_read_flight
from ..io.read_latmos_falcon import rasta_radar, rasta_microphysics

def read_flight(filename=None, platform='', file_format='netcdf', instrument=None,\
                initialize=False, mapping_dict=None):
    '''
    Takes a filename pointing to a aircraft flight level data file
     and returns a FlightLevel object.
        
    Parameters
    ----------
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

    reader = FileReader(filename=filename, platform=platform, file_format=file_format,
                        instrument=instrument, mapping_dict=mapping_dict)

    if initialize:
        FlightData = AirborneData(reader)

        return FlightData
    else:
        return reader


class FileReader(object):
    '''FileReader class to process data files.  '''

    def __init__(self, filename=None, platform=None, file_format=None, 
                instrument=None, mapping_dict=None):
        """
        If initialized with a filename (incl. path), will call
        ***_read_flight() to populate the class instance.
        If not, it simply instances the class but does not populate
        its attributes.
        verbose: Set to True for text output. Useful for debugging.

        Parameters
        ----------
        filename : str
            Filename (including path) of file to process
        platform : str
            Name of aicraft 
            Currently supported: 
                'p3' or p-3' (NOAA WP-3D)
                'citation' (Univ North Dakota Citation)
                'falcon' (LATMOS - SAFIRE)
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
                else:
                    print "Only netCDF format currently supported"

            elif (platform.upper() == 'CITATION'):
                if ((file_format.upper() == 'ASCII') or \
                (file_format.upper() == 'NASA AMES') or \
                (file_format.upper() == 'NA')):
                    flight = citation_read_flight(filename)
                elif file_format.upper() == 'NETCDF':
                    print "No netCDF reader currently exists for citation data"
                    return
                else:
                    print "Check the format call!"
                    return

            elif platform.upper() == 'FALCON':
                if instrument is None:
                    flight = latmos_falcon_read_flight(filename, mapping='basic')
                elif instrument.lower() == 'radar':
                    flight = rasta_radar(filename)
                elif instrument.lower() == 'microphysics':
                    flight = rasta_microphysics(filename)
                else:
                    print "Only netCDF format currently supported"

            if (instrument != 'radar') | (instrument != 'microphysics') :
                # Calculate meridional and zonal wind components
                Uwind, Vwind = self._winduv(flight)

                # Add to the dictionary
                flight['Uwind'] =  Uwind
                flight['Vwind'] =  Vwind
            else:
                flight['Uwind'] = flight['Vx_flight_level']#['Uwind']['data'][:]
                flight['Vwind'] = flight['Vy_flight_level']
                flight['Wwind'] = flight['Vz_flight_level']
                #['fields']['Vwind']['data'][:]

                # Record the data into the variable 
            self.flight_data = flight
        else:
            #Initializes class instance but leaves it to other methods to
            #populate the class attributes.
            return

    def _winduv(self,flight):
        """
        Calculate the horizontal windcomponents (u,v) 
        from wind angle and speed.
            U wind : positive blowing towards east
            V wind : positive blowing towards north
        """
        try:
            U = -np.cos(np.radians(flight['wind_dir'])) * flight['wind_spd']
        except:
            U = None
        try:
            V = -np.sin(np.radians(flight['wind_dir'])) * flight['wind_spd']
        except:
            V = None
        return U, V