"""
awot.io.read_p3_flight
=========================

This is a grouping of scripts designed to process NOAA P-3 
 flight level data distributed by NOAA AOC.
 
The data are in NetCDF format.

Author Nick Guy.  NRC, NOAA/NSSL (nick.guy@noaa.gov)
    Jul 2013  Originally developed in NCL 
    6 Aug 2014  Refactored to python for AWOT package

"""
# NOTES:: This has only been tested with DYNAMO data files, versions
#         may change and another function may be needed.
# FUNCTIONS::
#  flight_level_variable - Read in a variable from flight level NetCDF
#  flight_track - Read in data to for flight track
#-------------------------------------------------------------------
# Load the needed packages
from netCDF4 import Dataset, num2date
import numpy as np
import pytz
#-------------------------------------------------------------------
#===============================================================
# BEGIN FUNCTIONS
#===============================================================
def flight_data(fname):
    """
    Read in NetCDF data file containing P3 flight level data created
    by NOAA AOC.  
    Pull the needed variables for flight track info.
    
    PARAMETERS::
    ----------
        fname : string
            Filename [string]
    
    OUTPUT::
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
       
    USAGE::
    -----
     data = flight_track(fname)
    """
    # Read the NetCDF
    ncFile = Dataset(fname,'r')
    
    # Pull out each variable
    try:
        Lat = ncFile.variables['LatGPS.3'][:]
    except:
        Lat = None
    try:
        Lon = ncFile.variables['LonGPS.3'][:]
    except:
        Lon = None
    try:
        Alt = ncFile.variables['AltGPS.3'][:]
    except:
        Alt = None
    try:
        PAlt = ncFile.variables['AltPaADDU.1'][:]
    except:
        PAlt = None
    try:
        heading = ncFile.variables['THdgI-GPS.1'][:]
    except:
        heading = None
    try:
        track = ncFile.variables['TRK.1'][:]
    except:
        track = None
    try:
        vert_vel_DPJ = ncFile.variables['WSZ_DPJ'][:]
    except:
        vert_vel_DPJ = None
    try:
        temp = ncFile.variables['TA.1'][:]
    except:
        temp = None
    try:
        temp_dew = ncFile.variables['TD.1'][:]
    except:
        temp_dew = None
    try:
        temp_virt = ncFile.variables['TVIRT.1'][:]
    except:
        temp_virt = None
    try:
        theta = ncFile.variables['THETA.1'][:]
    except:
        theta = None
    try:
        thetaE = ncFile.variables['THETAE.1'][:]
    except:
        thetaE = None
    try:
        thetaV = ncFile.variables['THETAV.1'][:]
    except:
        thetaV = None
    try:
        wSpd = ncFile.variables['WS.1'][:]
    except:
        wSpd = None
    try:
        wDir = ncFile.variables['WD.1'][:]
    except:
        wDir = None
    try:
        RH = ncFile.variables['HUM_REL.1'][:]
    except:
        RH = None
    try:
        SH = ncFile.variables['HUM_SPEC.1'][:]
    except:
        SH = None
    try:
        mixing_ratio = ncFile.variables['MR.1'][:]
    except:
        mixing_ratio = None
    try:
        vap_press = ncFile.variables['EE.1'][:]
    except:
        vap_press = None
    try:
        sat_vap_press = ncFile.variables['EW.1'][:]
    except:
        sat_vap_press = None
        
    # Throw out an error message if file not read
    if Lat is None:
        print "Check the variable names in file!!"
    
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
    if Lat is not None:
        np.ma.masked_invalid(Lat)
    if Lon is not None:
        np.ma.masked_invalid(Lon)
    if Alt is not None:
        np.ma.masked_invalid(Alt)
    if PAlt is not None:
        np.ma.masked_invalid(PAlt)
    if heading is not None:
        np.ma.masked_invalid(heading)
    if track is not None:
        np.ma.masked_invalid(track)
    if vert_vel_DPJ is not None:
        np.ma.masked_invalid(vert_vel_DPJ)
    if temp is not None:
        np.ma.masked_invalid(temp)
    if temp_dew is not None:
        np.ma.masked_invalid(temp_dew)
    if temp_virt is not None:
        np.ma.masked_invalid(temp_virt)
    if theta is not None:
        np.ma.masked_invalid(theta)
    if thetaE is not None:
        np.ma.masked_invalid(thetaE)
    if thetaV is not None:
        np.ma.masked_invalid(thetaV)
    if wSpd is not None:
        np.ma.masked_invalid(wSpd)
    if wDir is not None:
        np.ma.masked_invalid(wDir)
    if RH is not None:
        np.ma.masked_invalid(RH)
    if SH is not None:
        np.ma.masked_invalid(SH)
    if mixing_ratio is not None:
        np.ma.masked_invalid(mixing_ratio)
    if vap_press is not None:
        np.ma.masked_invalid(vap_press)
    if sat_vap_press is not None:
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
    """
    Read in NetCDF data file containing P3 flight level data created
    by NOAA AOC.  
    Pull the needed variables for flight track info.
    
    PARAMETERS::
    ----------
     fname : string
         Filename [string]
    
    OUTPUT::
    ------
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
    if Lat is not None:
        np.ma.masked_invalid(Lat)
    if Lon is not None:
        np.ma.masked_invalid(Lon)
    if Alt is not None:
        np.ma.masked_invalid(Alt)
    if PAlt is not None:
        np.ma.masked_invalid(PAlt)
    
    # Pull out the start time
    StartTime = ncFile.StartTime
    
    # Create a time array 
    TimeSec = np.linspace(StartTime,StartTime + len(Lat), len(Lat))
    
    Time_unaware = num2date(TimeSec,'seconds since 1970-01-01 00:00:00+0:00')
    Time = Time_unaware#.replace(tzinfo=pytz.UTC)
    
    # Now mask missing values
    if Lat is not None:
        np.ma.masked_invalid(Lat)
    if Lon is not None:
        np.ma.masked_invalid(Lon)
    if Alt is not None:
        np.ma.masked_invalid(Alt)
    if PAlt is not None:
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
    """
    Read in NetCDF data file containing P3 flight level data created
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
    ------
     VarOut : float
         Masked array containing variable data
    
    USAGE::
    -----
     Lat = read_flight_level_dynamo('P3.nc','LatGPS.3')
    
    NOTES::
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

