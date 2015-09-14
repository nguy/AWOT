# -*- coding: utf-8 -*-
"""
RadarDataFile.py - Class for reading Data from file
"""
import os

#from ..AirborneData import AirborneData
from ..io.read_ground_radar import read_radar as read_ground_radar
from ..io.read_p3_radar import read_lf_grid, read_windsyn_tdr_netcdf 
from ..io.read_p3_radar import read_tdr_sweep, read_windsyn_binary
from ..io.read_latmos_falcon import rasta_radar, rasta_microphysics
#from ..io.read_uwka import read_wcr2


def read_radar(filename=None, platform='p3', file_format='netcdf', instrument=None,
               initialize=False ):
    '''
    Takes a filename pointing to a aircraft radar data file
     and returns a RadarData object.

    Parameters
    ----------
    filename : str
        Long path filename of data file.
    platform : str
        Platform for processing, see FileReader.
        Currenlty supported:
            'p3' - NOAA P-3 radar tail Doppler and lower fuselage radars
            'falcon' - LATMOS French falcon W-band
            'king air' - University of Wyoming King Air W-band ProSensing radar
    file_format : str
        Format of input file, see FileReader.
    instrument : str
        Type (name) of instrument to process
        Currently the following arguments are valid:
            'tdr_grid' - Tail Doppler radar gridded files (e.g. dual-Doppler analysis)
            'tdr_sweep' = Tail Doppler radar (Native coordinate data)
            'lf'  - Lower Fuselage radar
            'ground' - A ground-based radar system, read in using PyArt
            'wcr' - University of Wyoming King Air W-band ProSensing radar
    data_format : str
        Either 'grid' or 'native'.

    To Do:
    Add support for:
        Eldora X-band (historical) - should be similar to P3
        King Air/C-130 95 GHZ (Wyoming Cloud Radar; WCR)
    '''
    if filename is None:
        print "Need to provide a filename!"
        return

    reader = FileReader(filename=filename, platform=platform, 
                        file_format=file_format, instrument=instrument)

    if initialize:
        RadarData = AirborneData(reader)

        return RadarData
    else:
        return reader


class FileReader(object):
    '''FileReader class to process data files.'''
    def __init__(self, filename=None, platform=None, file_format='netcdf',
                 instrument=None):
        """
        If initialized with a filename (incl. path), will call
        ***_read_flight() to populate the class instance.
        If not, it simply instances the class but does not populate
        its attributes.
        verbose: Set to True for text output. Useful for debugging.

        Parameters
        ----------
            filename : string
                Filename (including path) of file to process
            platform : string
                Name of aicraft 
                Currently supported: 
                    'p3' or p-3' (NOAA WP-3D)
                    'eldora'
                    'citation' (Univ North Dakota Citation)
                    'kingair' or 'king air' (Univ of Wyoming King Air)
            file_format : string
                Format of input file.  Each platform currently has
                a specific file type.  This may be extended in the future.
                Currently supported:
                    'netcdf'
                    'ascii'
            instrument : str
                'tdr_grid' - Tail Doppler radar gridded files (e.g. dual-Doppler analysis)
                'tdr_sweep' = Tail Doppler radar (Native coordinate data)
                'lf'  - Lower Fuselage radar
                'ground' - A ground-based radar system, read in using PyArt
                'wcr' - Wyoming Cloud Radar, 
                        University of Wyoming King Air W-band ProSensing radar
        """
        if isinstance(filename, str) != False:
            if (platform.upper() == 'P3') or (platform.upper() == 'P-3') or \
            (platform.upper() == 'NOAA_P3') or \
            (platform.upper() == 'NOAA P3') or \
            (platform.upper() == 'ELDORA'):
                if instrument.lower() == 'tdr_grid':
                    if file_format.lower() == 'netcdf':
                        radar = read_windsyn_tdr_netcdf(filename)
                    elif file_format.lower() == 'binary':
                        radar = read_windsyn_binary(filename, 
                                                    instrument=instrument,
                                                    platform=platform)
                elif instrument.lower() == 'tdr_sweep':
                    if (file_format.lower() =='netcdf') or \
                       (file_format.lower() =='sigmet'):
                        radar = read_tdr_sweep(filename)
                    else:
                        print ("Currently limited to netCDF and Sigmet " +\
                                + "formats for sweep files")
                elif instrument.lower() == 'lf':
                    radar = read_lf_grid(filename)

            elif platform.upper() == 'FALCON':
                if instrument.lower() == 'radar':
                    radar = rasta_radar(filename)
                elif instrument.lower() == 'microphysics':
                    radar = rasta_microphysics(filename)

            elif (platform.upper() == 'C130') or \
            (platform.upper() == 'NCAR_C130') or \
            (platform.upper() == 'NCAR C130'):
                print "Sorry not supported at this time"

            elif (platform.upper() == 'KING AIR') or \
            (platform.upper() == 'KING_AIR') or \
            (platform.upper() == 'KINGAIR') or \
            (platform.upper() == 'KING-AIR') or \
            (platform.upper() == 'WCR'):
                radar = read_wcr2(filename)

            elif (platform.upper() == 'GROUND'):
                radar = read_ground_radar(filename)

            else:
                # Check to see if a ground instrument is being fed
                if instrument.lower() == 'ground':
                    radar = read_ground_radar(filename)

                elif (instrument.lower() == 'wcr'):
                    radar = read_wcr2(filename)

                else:
                    print "Check supported platform list: \
                           'p3' or 'p-3' - NOAA P-3 \
                           'eldora' - NCAR Eldor radar (same as P-3) \
                           'falcon' - LATMOS French Falcon \
                           'C130' - NCAR C130, Coming Soon \
                           'King Air' or 'King_Air' or 'WCR' - Wyoming Cloud Radar \
                           'Ground' - Any PyArt supported data format"
                    return

            # Record the data into the variable 
            self.radar_data = radar

        else:
            #Initializes class instance but leaves it to other methods to
            #populate the class attributes.
            return