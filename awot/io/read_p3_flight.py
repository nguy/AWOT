"""
awot.io.read_p3_flight
=========================

This is a grouping of scripts designed to process NOAA P-3 
 flight level data recorded during flights and put into NetCDF
 format by NOAA AOC.

Created by Nick Guy.

Original code developed in NCL between Jul 2013 - Mar 2014, 
refactored 6 Aug 2014 to python

"""
# NOTES:: This has only been tested with DYNAMO data files, versions
#         may change and another function may be needed.
# HISTORY::
#   8 Jan 2014 - Nick Guy.   NRC, NOAA/NSSL (nick.guy@noaa.gov)   
#                Converted NCL functions below to Python
# FUNCTIONS::
#  flight_level_variable - Read in a variable from flight level NetCDF
#  flight_track - Read in data to for flight track
#-------------------------------------------------------------------
# Load the needed packages
from netCDF4 import Dataset,num2date
import numpy as np
import pytz
#-------------------------------------------------------------------
# Define various constants that may be used for calculations
#===============================================================
# BEGIN FUNCTIONS
#===============================================================
def flight_data(fname):
    """Read in data from NetCDF file containing P3 flight level data created
    by NOAA AOC.  Pull out the needed variables for flight track info.
    PARAMETERS::
    ----------
        fname : string
            Filename [string]
    OUTPUT::
    ----------
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
       
    USAGE::
     data = flight_track(fname)
    """
    # Read the NetCDF
    ncFile = Dataset(fname,'r')
    
    # Pull out each variable
    Lat = ncFile.variables['LatGPS.3'][:]
    Lon = ncFile.variables['LonGPS.3'][:]
    Alt = ncFile.variables['AltGPS.3'][:]
    PAlt = ncFile.variables['AltPaADDU.1'][:]
    heading = ncFile.variables['THdgI-GPS.1'][:]
    track = ncFile.variables['TRK.1'][:]
    vert_vel_DPJ = ncFile.variables['WSZ_DPJ'][:]
    temp = ncFile.variables['TA.1'][:]
    temp_dew = ncFile.variables['TD.1'][:]
    temp_virt = ncFile.variables['TVIRT.1'][:]
    theta = ncFile.variables['THETA.1'][:]
    thetaE = ncFile.variables['THETAE.1'][:]
    thetaV = ncFile.variables['THETAV.1'][:]
    wSpd = ncFile.variables['WS.1'][:]
    wDir = ncFile.variables['WD.1'][:]
    RH = ncFile.variables['HUM_REL.1'][:]
    SH = ncFile.variables['HUM_SPEC.1'][:]
    mixing_ratio = ncFile.variables['MR.1'][:]
    vap_press = ncFile.variables['EE.1'][:]
    sat_vap_press = ncFile.variables['EW.1'][:]
    
    # Pull out the start time
    StartTime = ncFile.StartTime
    
    # Create a time array 
    TimeSec = np.linspace(StartTime,StartTime + len(Lat), len(Lat))
    
    Time_unaware = num2date(TimeSec,'seconds since 1970-01-01 00:00:00+0:00')
    Time = Time_unaware#.replace(tzinfo=pytz.UTC)
    
    # Pull out global attributes
    try:
        ncFile.ProjectName
        project = ncFile.ProjectName
    except:
        project = None
    try:
        ncFile.Platform
        platform = ncFile.Platform
    except:
        platform = None
    try:
        ncFile.FlightNumber
        flightnum = ncFile.FlightNumber
    except:
        flightnum = None
    
    # Now mask missing values
    np.ma.masked_invalid(Lat)
    np.ma.masked_invalid(Lon)
    np.ma.masked_invalid(Alt)
    np.ma.masked_invalid(PAlt)
    np.ma.masked_invalid(heading)
    np.ma.masked_invalid(track)
    np.ma.masked_invalid(vert_vel_DPJ)
    np.ma.masked_invalid(temp)
    np.ma.masked_invalid(temp_dew)
    np.ma.masked_invalid(temp_virt)
    np.ma.masked_invalid(theta)
    np.ma.masked_invalid(thetaE)
    np.ma.masked_invalid(thetaV)
    np.ma.masked_invalid(wSpd)
    np.ma.masked_invalid(wDir)
    np.ma.masked_invalid(RH)
    np.ma.masked_invalid(SH)
    np.ma.masked_invalid(mixing_ratio)
    np.ma.masked_invalid(vap_press)
    np.ma.masked_invalid(sat_vap_press)

    # Create a dictionary to transfer the data
    data = {'latitude': Lat,
            'longitude': Lon,
            'altitude': Alt,
            'pressure_altitude': PAlt,
            'time': Time,
            'true_heading': heading,
            'track': track,
            'vert_vel_DPJ': vert_vel_DPJ,
            'temperature': temp,
            'dewpoint_temperature': temp_dew,
            'virtual_temperature': temp_virt,
            'potential_temp': theta,
            'equiv_potential_temp': thetaE,
            'virtual_potential_temp': thetaV,
            'wind_spd': wSpd,
            'wind_dir': wDir,
            'relative_humidity': RH,
            'specific_humidity': SH,
            'mixing_ratio': mixing_ratio,
            'vapor_pressure': vap_press,
            'sat_vapor_pressure': sat_vap_press,
            'project' : project,
            'platform': platform,
            'flight_number': flightnum,
            }
    
    ncFile.close()
    
    return data
    
#**====================================================

def flight_track(fname):
    """Read in data from NetCDF file containing P3 flight level data created
    by NOAA AOC.  Pull out the needed variables for flight track info.
    PARAMETERS::
    ----------
     fname : string
         Filename [string]
    OUTPUT::
    ----------
     data : Dictionary of the following values
       Lat : float
           Aircraft latitude
       Lon : float
           Aircraft longitude
       Alt : float
           Aircraft altitude
       PAlt : float
           Aircraft pressure altitude
       Time : float
           Aircraft time array
    USAGE::
    ----------
     data = flight_track(fname)
    """
    # Read the NetCDF
    ncFile = Dataset(fname,'r')
    
    # Pull out each variable
    Lat = ncFile.variables['LatGPS.3'][:]
    Lon = ncFile.variables['LonGPS.3'][:]
    Alt = ncFile.variables['AltGPS.3'][:]
    PAlt = ncFile.variables['AltPaADDU.1'][:]
    
    # Pull out the start time
    StartTime = ncFile.StartTime
    
    # Create a time array 
    TimeSec = np.linspace(StartTime,StartTime + len(Lat), len(Lat))
    
    Time_unaware = num2date(TimeSec,'seconds since 1970-01-01 00:00:00+0:00')
    Time = Time_unaware#.replace(tzinfo=pytz.UTC)
    
    # Now mask missing values
    np.ma.masked_invalid(Lat)
    np.ma.masked_invalid(Lon)
    np.ma.masked_invalid(Alt)
    np.ma.masked_invalid(PAlt)

    # Create a dictionary to transfer the data
    data = {'latitude': Lat,
            'longitude': Lon,
            'altitude': Alt,
            'pressure_altitude': PAlt,
            'time': Time}
            
    ncFile.close()
    
    return data
    
#**====================================================
    
def flight_level_variable(fname,Rec):
    """Read in data from NetCDF file containing P3 flight level data created
    by NOAA AOC.  The NetCDF should be read in the main program and passed
    to this function.
    A call such as this can be used in the main program:
      FltncID=addfile(FlightFileStringName,"r")
    PARAMETERS::
    ----------
     fname : string
         Filename [string]
     Rec : string
         Variable name to be pulled out [string]
    OUTPUT::
    ----------
     VarOut : float
         Masked array containing variable data
    USAGE::
    ----------
     Lat = read_flight_level_dynamo('P3.nc','LatGPS.3')
    NOTES::
    ----------
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

    # Mask any "_FillValue" or some missing_data type attribute
    #try:
    #    VarOut = np.ma.masked_values(VarOut, VarOut.missing_value)
    #except:
    #    pass
    #try:
    #    VarOut = np.ma.masked_values(VarOut, VarOut.Missing_Value)
    #except:
    #    pass
    #try:
    #    VarOut = np.ma.masked_values(VarOut, VarOut._FillValue)
    #except:
    #    pass

    # Mask any NaN values
    #VarOut = np.ma.masked_values(VarOut, np.isnan(VarOut)
    
    ncFile.close()

    return VarOut

