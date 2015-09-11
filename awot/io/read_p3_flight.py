"""
awot.io.read_p3_flight
=========================

This is a grouping of scripts designed to process NOAA P-3 
 flight level data distributed by NOAA AOC.
"""
# Load the needed packages
from netCDF4 import Dataset, num2date
import numpy as np
import pytz

def flight_data(fname, mapping_dict=None):
    """
    Read in NetCDF data file containing P3 flight level data created
    by NOAA AOC.  
    Pull the needed variables for flight track info.
    
    Parameters
    ----------
    fname : string
        Filename [string]
    mapping_dict : dict
        Dictionary to use for mapping variable names.
        Use if provided.
        If None, use default.

    Output
    ------
        data : Dictionary of the following values
            Lat : float
                Aircraft latitude [deg]
            Lon : float
                Aircraft longitude [deg]
            Alt : float
                Aircraft altitude [m]
            PAlt : float
                Aircraft pressure altitude [m]
            Time : float
                Aircraft time array
            heading : float
                Aircraft true heading [deg]
            track : float
                Aircraft track [deg]
            vert_vel_DPJ : float
                Vertical wind via D Jorgensen calculation [m/s]
            temp : float
                Ambient Temperature [C]
            temp_dew : float
                Dewpoint Temperature [C]
            temp_virt : float
                Virtual Temperature [K]
            theta : float
                Potential Temperature [K]
            thetaE : float
                Equivalent Potential Temperature [K]
            thetaV : float
                Virtual Potential Temperature [K]
            wSpd : float
                Wind Speed [m/s]
            wDir : float
                Wind Direction [deg]
            RH : float
                Relative Humidity [%]
            SH : float
                Specific Humidity [g/kg]
            mixing_ratio : float
                Mixing ratio [g] [g/g?]
            vap_press : float
                Vapor Pressure [hPa]
            sat_vap_press : float
                Saturated Vapor Pressure [hPa]
            project : str
                Project name
            platform : str
                Platform name
            flightnum : str
                Flight number
    """
    # Read the NetCDF
    ncFile = Dataset(fname,'r')

    if mapping_dict is None:
        name_map = _get_p3_flight_namemap()
    else:
        name_map = mapping_dict

    # Cycle to through variables in file
    data = {}
    for var in name_map:
        try:
            data[var] = ncFile.variables[name_map[var]][:]
            np.ma.masked_invalid(data[var])
        except:
            data[var] = None

    # Throw out an error message if file not read
    if data['latitude'] is None:
        print "Check the variable names in file!!"

    # Get the time
    Time = _get_time(ncFile)
    data['time'] = Time

    # Pull out global attributes
    try:
        ncFile.ProjectName
        data['project'] = ncFile.ProjectName
    except:
        data['project'] = None
    try:
        ncFile.Platform
        data['platform'] = ncFile.Platform
    except:
        data['platform'] = None
    try:
        ncFile.FlightNumber
        data['flight_number'] = ncFile.FlightNumber
    except:
        data['flight_number'] = None

    ncFile.close()
    return data

def flight_level_variable(fname,Rec):
    """
    Read in NetCDF data file containing P3 flight level data created
    by NOAA AOC.  The NetCDF should be read in the main program and passed
    to this function.
    A call such as this can be used in the main program:
      FltncID=addfile(FlightFileStringName,"r")

    Parameters
    ----------
     fname : string
         Filename [string]
     Rec : string
         Variable name to be pulled out [string]

    Output
    ------
     VarOut : float
         Masked array containing variable data

    Notes
    -----
    Data file structure::
     Available variables (not full list) :
     LonGPS.3      = Novatel GPS Longitude
     LatGPS.3      = Novatel GPS Latitude
     AltGPS.3      = Novatel GPS Altitude [m]
     THdgI-GPS.1   = True heading [deg]
     TRK.1         = Track [deg]
     AltPaADDU.1   = Pressure altitude [m]
     WSZ_DPJ       = Vertical wind via D Jorgensen calculation [m/s]
     TA.1          = Ambient Temperature [C]
     TD.1          = Dewpoint Temperature [C]
     TVIRT.1       = Virtual Temperature [K]
     THETA.1       = Potential Temperature [K]
     THETAE.1      = Equivalent Potential Temperature [K]
     THETAV.1      = Virtual Potential Temperature [K]
     WS.1          = Wind Speed [m/s]
     WD.1          = Wind Direction [deg]
     HUM_REL.1     = Relative Humidity [%]
     HUM_SPEC.1    = Specific Humidity [g/kg]
     MR.1          = Mixing ratio [g] [g/g?]
     EE.1          = Vapor Pressure [hPa]
     EW.1          = Saturated Vapor Pressure [hPa]
      """
    # Read the NetCDF
    ncFile = Dataset(fname,'r')

    # Get the variable of interest
    VarOut = ncFile.variables[Rec][:]

    ncFile.close()
    return VarOut

def _get_p3_flight_namemap():
    '''Map NOAA P3 variables to AWOT structure'''
    name_map = {
               'latitude': 'LatGPS.3',
               'longitude': 'LonGPS.3',
               'altitude': 'AltGPS.3',
               'pressure_altitude': 'AltPaADDU.1',
               'true_heading': 'THdgI-GPS.1',
               'track': 'TRK.1',
               'vert_vel_DPJ': 'WSZ_DPJ',
               'temperature': 'TA.1',
               'dewpoint_temperature': 'TD.1',
               'virtual_temperature': 'TVIRT.1',
               'potential_temp': 'THETA.1',
               'equiv_potential_temp': 'THETAE.1',
               'virtual_potential_temp': 'THETAV.1',
               'wind_spd': 'WS.1',
               'wind_dir': 'WD.1',
               'relative_humidity': 'HUM_REL.1',
               'specific_humidity': 'HUM_SPEC.1',
               'mixing_ratio': 'MR.1',
               'vapor_pressure': 'EE.1',
               'sat_vapor_pressure': 'EW.1',
               }
    return name_map

def _get_time(ncFile):
    # Pull out the start time
    try:
        TimeSec = ncFile.variables['base_time'][:]
    except:
        StartTime = ncFile.StartTime
        length = len(ncFile.dimensions['Time'])
        # Create a time array 
        TimeSec = np.linspace(StartTime,StartTime + length, length)

    Time_unaware = num2date(TimeSec,'seconds since 1970-01-01 00:00:00+0:00')
    Time = Time_unaware
    return Time