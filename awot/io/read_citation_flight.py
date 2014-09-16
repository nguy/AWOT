"""
awot.io.read_citation_flight
=========================

This is a grouping of scripts designed to process University of
 North Dakota Citation flight level data recorded during flights 
 and NASA Ames Ascii format

Created by Nick Guy.

Original code developed in NCL between Jul 2013 - Mar 2014, 
refactored 6 Aug 2014 to python

"""
# NOTES:: This has only been tested with DYNAMO data files, versions
#         may change and another function may be needed.
# HISTORY::
#  21 Aug 2014 - Nick Guy.   NRC, NOAA/NSSL (nick.guy@noaa.gov) 
#-------------------------------------------------------------------
# Load the needed packages
import nappy
from netCDF4 import Dataset, num2date, date2num
import numpy as np
import pytz
#-------------------------------------------------------------------
# Define various constants that may be used for calculations
#===============================================================
# BEGIN FUNCTIONS
#===============================================================
def flight_data(fname):
    """Read in data from NASA Ames formatted ASCII file containing the
        University of North Dakota Citation aircraft flight level data
    INPUT::
        fname : string
            Filename [string]
    OUTPUT::
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
#---------------------------------------------------
    # Read the file
    FileIn = nappy.openNAFile(fname)
    
    # Make sure this is really a UND Citation file, always 1001
    if FileIn.getFFI() != 1001:
        print "Check if this is a UND Citation file!"
        FilenIn.close()
        return
        
    # Get the whole file as a dictionary (NA = NASA Ames)
    NA_dict = FileIn.getNADict()
    
    # Pull out global attributes
    project = NA_dict['MNAME']
    platform = NA_dict['SNAME']
    comments = NA_dict['SCOM']
    missingData = NA_dict['VMISS']
    scaleFactor = NA_dict['VSCAL']
    hdr_row_length = NA_dict['NLHEAD']

    # Now add a value to beginning of missing data scale factor arrays
    # to account for time column
    missingData.insert(0, -99.)
    scaleFactor.insert(0, 1.)

    # Get the data columns
    junk = np.genfromtxt(fname, skiprows=hdr_row_length, 
                          missing_values=missingData, filling_values=np.nan)
    
    # Get the variable names
    # This assumes that there are 4 lines of normal comments and that the 
    VarNames = NA_dict['NCOM'][2].split()
    
    # Pull out each variable data and multiply by scale factor
    Lat = np.array(junk[:,VarNames.index("POS_Lat")] \
               * scaleFactor[VarNames.index("POS_Lat")])
    Lon = np.array(junk[:,VarNames.index("POS_Lon")] \
                * scaleFactor[VarNames.index("POS_Lon")])
    Alt = np.array(junk[:,VarNames.index("POS_Alt")] \
                * scaleFactor[VarNames.index("POS_Alt")])
    PAlt = np.array(junk[:,VarNames.index("Press_Alt")] \
                * scaleFactor[VarNames.index("Press_Alt")])
    heading = np.array(junk[:,VarNames.index("POS_Head")] \
                * scaleFactor[VarNames.index("POS_Head")])
    track =np.array(junk[:,VarNames.index("POS_Trk")] \
                * scaleFactor[VarNames.index("POS_Trk")])
    vert_vel = np.array(junk[:,VarNames.index("Wind_Z")] \
                * scaleFactor[VarNames.index("Wind_Z")])
    temp = np.array(junk[:,VarNames.index("Air_Temp")] \
                * scaleFactor[VarNames.index("Air_Temp")])
    temp_dew = np.array(junk[:,VarNames.index("DEWPT")] \
                * scaleFactor[VarNames.index("DEWPT")])
    theta = np.array(junk[:,VarNames.index("Pot_Temp_T1")] \
                * scaleFactor[VarNames.index("Pot_Temp_T1")])
    wSpd = np.array(junk[:,VarNames.index("Wind_M")] \
                * scaleFactor[VarNames.index("Wind_M")])
    wDir = np.array(junk[:,VarNames.index("Wind_D")] \
                * scaleFactor[VarNames.index("Wind_D")])
    mixing_ratio = np.array(junk[:,VarNames.index("MixingRatio")] \
                * scaleFactor[VarNames.index("MixingRatio")])
    TfrostPt= np.array(junk[:,VarNames.index("FrostPoint")] \
                * scaleFactor[VarNames.index("FrostPoint")])
    Turb = np.array(junk[:,VarNames.index("TURB")] \
                * scaleFactor[VarNames.index("TURB")])
    LWC_King = np.array(junk[:,VarNames.index("King_LWC_ad")] \
                * scaleFactor[VarNames.index("King_LWC_ad")])
    TWC_Nev = np.array(junk[:,VarNames.index("Nev_TWC")] \
                * scaleFactor[VarNames.index("Nev_TWC")])
    LWC_CDP = np.array(junk[:,VarNames.index("CDP_LWC")] \
                * scaleFactor[VarNames.index("CDP_LWC")])
    Nt_2DC = np.array(junk[:,VarNames.index("2-DC_Conc")] \
                * scaleFactor[VarNames.index("2-DC_Conc")])
    Dmean_2DC = np.array(junk[:,VarNames.index("2-DC_MenD")] \
                * scaleFactor[VarNames.index("2-DC_MenD")])
    Dvol_2DC = np.array(junk[:,VarNames.index("2-DC_VolDia")] \
                * scaleFactor[VarNames.index("2-DC_VolDia")])
    Deff_2DC = np.array(junk[:,VarNames.index("2-DC_EffRad")] \
                * scaleFactor[VarNames.index("2-DC_EffRad")])
    Nt_CPC = np.array(junk[:,VarNames.index("CPCConc")] \
                * scaleFactor[VarNames.index("CPCConc")])
    
    # Create a time array 
    TimeSec = junk[:,VarNames.index("Time")]
    
    # Get the date
    Date = NA_dict['DATE']
    DateStr = str(Date[0])+'-'+str(Date[1])+'-'+str(Date[2])
    
    # Create a datetime instance from file
    Time_dt = num2date(TimeSec, 'seconds since '+DateStr+' 00:00:00+0:00')

    # Now move this back into number format
    TimeNum = date2num(Time_dt, 'seconds since 1970-1-1 00:00:00+0:00')
    
    # Finally convert this back to a standard used by this package (Epoch)
    Time_unaware = num2date(TimeNum,'seconds since 1970-1-1 00:00:00+0:00')
    Time = Time_unaware#.replace(tzinfo=pytz.UTC)

    # Now mask missing values
    Lat = np.ma.masked_values(Lat, missingData[VarNames.index("POS_Lat")])
    Lon = np.ma.masked_values(Lon, missingData[VarNames.index("POS_Lon")])
    Alt = np.ma.masked_values(Alt, missingData[VarNames.index("POS_Alt")])
    PAlt = np.ma.masked_values(PAlt, missingData[VarNames.index("Press_Alt")])
    heading = np.ma.masked_values(heading, missingData[VarNames.index("POS_Head")])
    track = np.ma.masked_values(track, missingData[VarNames.index("POS_Trk")])
    vert_vel = np.ma.masked_values(vert_vel, missingData[VarNames.index("Wind_Z")])
    temp = np.ma.masked_values(temp, missingData[VarNames.index("Air_Temp")])
    temp_dew = np.ma.masked_values(temp_dew, missingData[VarNames.index("DEWPT")])
    theta = np.ma.masked_values(theta, missingData[VarNames.index("Pot_Temp_T1")])
    wSpd = np.ma.masked_values(wSpd, missingData[VarNames.index("Wind_M")])
    wDir = np.ma.masked_values(wDir, missingData[VarNames.index("Wind_D")])
    mixing_ratio = np.ma.masked_values(mixing_ratio, missingData[VarNames.index("MixingRatio")])
    TfrostPt = np.ma.masked_values(TfrostPt, missingData[VarNames.index("FrostPoint")])
    Turb = np.ma.masked_values(Turb, missingData[VarNames.index("TURB")])
    LWC_King = np.ma.masked_values(LWC_King, missingData[VarNames.index("King_LWC_ad")])
    TWC_Nev = np.ma.masked_values(TWC_Nev, missingData[VarNames.index("Nev_TWC")])
    LWC_CDP = np.ma.masked_values(LWC_CDP, missingData[VarNames.index("CDP_LWC")])
    Nt_2DC = np.ma.masked_values(Nt_2DC, missingData[VarNames.index("2-DC_Conc")])
    Dmean_2DC = np.ma.masked_values(Dmean_2DC, missingData[VarNames.index("2-DC_MenD")])
    Dvol_2DC = np.ma.masked_values(Dvol_2DC, missingData[VarNames.index("2-DC_VolDia")])
    Deff_2DC = np.ma.masked_values(Deff_2DC, missingData[VarNames.index("2-DC_EffRad")])
    Nt_CPC = np.ma.masked_values(Nt_CPC, missingData[VarNames.index("CPCConc")])

    # Create a dictionary to transfer the data
    data = {'latitude': Lat,
            'longitude': Lon,
            'altitude': Alt,
            'pressure_altitude': PAlt,
            'time': Time,
            'true_heading': heading,
            'track': track,
            'vert_vel': vert_vel,
            'temperature': temp,
            'dewpoint_temperature': temp_dew,
            'potential_temp': theta,
            'wind_spd': wSpd,
            'wind_dir': wDir,
            'mixing_ratio': mixing_ratio,
            'frostpoint_temperature': TfrostPt,
            'turbulence_parameter': Turb,
            'LWC_king': LWC_King,
            'TWC_Nevzorov': TWC_Nev,
            'LWC_CDP': LWC_CDP,
            'total_conc_2DC': Nt_2DC,
            'drop_diam_mean_2DC': Dmean_2DC,
            'drop_diam_mean_vol_2DC': Dvol_2DC,
            'drop_eff_radius': Deff_2DC,
            'total_conc_CPC': Nt_CPC,
            'project' : project,
            'platform': platform,
            'flight_number': Date,
            }
            
    FileIn.close()

    return data 
    
#**====================================================

def flight_track(fname):
    """Read in data from NetCDF file containing P3 flight level data created
    by NOAA AOC.  Pull out the needed variables for flight track info.
    INPUT::
     fname : string
         Filename [string]
    OUTPUT::
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
     data = flight_track(fname)
    """
#---------------------------------------------------
    # Read the NetCDF
    ncFile = Dataset(fname,'r')
    
    # Pull out each variable
    Lat = np.array(junk[:,VarNames.index("POS_Lat")])
    Lon = np.array(junk[:,VarNames.index("POS_Lon")])
    Alt = np.array(junk[:,VarNames.index("POS_Alt")])
    PAlt = np.array(junk[:,VarNames.index("Press_Alt")])
    
    # Create a time array 
    TimeSec = junk[:,VarNames.index("Time")]
    
    # Get the date
    Date = NA_dict['DATE']
    DateStr = str(Date[0])+'-'+str(Date[1])+'-'+str(Date[2])
    
    # Create a datetime instance from file
    Time_dt = num2date(TimeSec, 'seconds since '+DateStr+' 00:00:00+0:00')
    
    # Now move this back into number format
    TimeNum = date2num(Time_dt, 'seconds since 1970-1-1 00:00:00+0:00')
    
    # Finally convert this back to a standard used by this package (Epoch)
    Time_unaware = num2date(TimeSec,'seconds since 1970-1-1 00:00:00+0:00')
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
    
    return data

