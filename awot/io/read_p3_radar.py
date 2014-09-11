"""
awot.io.read_p3_radar
=========================

A group of scripts to read data collected by the NOAA P-3 aircraft.
Supports both tail Doppler and lower fuselage radars. 

Created by Nick Guy.

"""
# NOTES:: This has only been tested with DYNAMO data files, versions
#         may change and another function may be needed.
# HISTORY::
#   8 Jan 2014 - Nick Guy.   NRC, NOAA/NSSL (nick.guy@noaa.gov)   
#                Converted NCL functions below to Python
#-------------------------------------------------------------------
# Load the needed packages
from netCDF4 import Dataset
import numpy as np
from Nio import open_file
import pytz
from datetime import datetime
import pyart.io as pio
#-------------------------------------------------------------------
# Define various constants that may be used for calculations
#===============================================================
# BEGIN FUNCTIONS
#===============================================================
def read_dpj_tdr_grid_netcdf(fname):
    """Read in data from NetCDF file containing P3 flight level data created
    by NOAA AOC.  The NetCDF should be read in the main program and passed
    to this function.
        
    Parameters::
    ----------
    fname : str
        Long path filename [string]
        
    Output::
    ----------
    data : Dictionary of the following values
    
    metadata : Dictionary of global attributes in file
    longitude : dict
        Longitude array [degrees]
    latitude : dict
        Latitude array [degrees]
    height : dict
        Height array [km]
    fields : Dictionary of variables in file
        reflectivity : float
            Radar Reflectivity [dBZ]
        Uwind : float
            Wind along aircraft longitudinal axis wind [m/s]
        Vwind : float
            Wind perpendicular to aircraft longitudinal axis wind [m/s]
        Wwind : float
            Vertical wind component [m/s]
        divergene : float
            Divergence x 10^3 [1/s]
        term_fall_speed : float
            Terminal fall speed [m/s]
        time_diff : float
            Time difference from nominal [s]
        datetime_start : datetime object
            HHNN (Hour Minute) module start time
        datetime_end : datetime object
            HHNN (Hour Minute) module end time
    """
#=====================================
    # Read the NetCDF
    ncFile = Dataset(fname,'r')
    ncvars = ncFile.variables
#--------------------------------------
    # Grab the metadata stored in global attributes as a dictionary
    metadata = ncFile.__dict__
    
    # Find the module start information
    # Year
    yys = ncFile.variables['YYstart'][:]
    # Month	
    mms = ncFile.variables['MMstart'][:]
    # Day
    dds = ncFile.variables['DDstart'][:]
    # Hour
    hhs = ncFile.variables['HHstart'][:]
    # Minute
    nns = ncFile.variables['MIstart'][:]
    # Second
    sss = ncFile.variables['SSstart'][:]
    #--------------------------------
    # Find the module end information
    # Year
    yye = ncFile.variables['YYend'][:]
    # Month	
    mme = ncFile.variables['MMend'][:]
    # Day
    dde = ncFile.variables['DDend'][:]
    # Hour
    hhe = ncFile.variables['HHend'][:]
    # Minute
    nne = ncFile.variables['MIend'][:]
    # Second
    sse = ncFile.variables['SSend'][:]

    # Make an output string for Date-Time start/finish    	
#    dt_start = yys+'-'+mms+'-'+dds+' '+hhs+':'+nns+':'+sss 
#    dt_end = yye+'-'+mme+'-'+dde+' '+hhe+':'+nne+':'+sse  
    # Make a Datetime instance to pass
    utc = pytz.utc
    dt_start = datetime(yys,mms,dds,hhs,nns,sss,0)#,tzinfo=utc)
    dt_end = datetime(yye,mme,dde,hhe,nne,sse,0)#,tzinfo=utc)
    #--------------------------------------
    # Pull in the Lat / Lon information
    Lon = _ncvar_to_dict(ncvars['Lon'])
    Lat = _ncvar_to_dict(ncvars['Lat'])
    
    # Pull in the Height 
    Ht = _ncvar_to_dict(ncvars['Height'])
    #---------------------------------------
    # Add fields to their own dictionary
    #---------------------------------------
    fields = {}
    # Read the Reflectivity
    fields['dBZ'] = _ncvar_to_dict(ncvars['dBZ'])
    
    # Read the U wind (Along aircraft longitudinal axis)
    fields['Uwind'] = _ncvar_to_dict(ncvars['U'])
    
    # Read the V wind (Perpendicular to aircraft longitudinal axis)
    fields['Vwind'] = _ncvar_to_dict(ncvars['V'])
    
    # Read the vertical wind
    fields['Wwind'] = _ncvar_to_dict(ncvars['W'])
    
    # Try reading some variables that are not always included
    # Read the divergence
    try:
        fields['divergence'] = _ncvar_to_dict(ncvars['Div'])
    except:
        fields['divergence'] = None
    
    # Read the terminal velocity
    try:
        fields['term_fall_speed'] = _ncvar_to_dict(ncvars['Vt'])
    except:
        fields['term_fall_speed'] = None
    
    # Read the time difference
    try:
        fields['time_diff'] = _ncvar_to_dict(ncvars['Tdiff'])
    except:
        fields['time_diff'] = None
    
    # Mask bad data values
#    badval = float(ncFile.MissingValue)
#    np.ma.masked_values(fields['dBZ']['data'], badval)
    # Pull out specific information on platform/instrument   
    try:
        ncFile.Platform
        platform = ncFile.title
    except:
        platform = 'p3'    

    try:
        ncFile.Instrument
        instrument = ncFile.instrument
    except:
        instrument = 'tdr_grid'
    
    # Create a dictionary to transfer the data
    data = {'metadata' : metadata,
            'longitude': Lon,
            'latitude': Lat,
            'height': Ht,
            'fields': fields,
            'datetime_start' : dt_start,
            'datetime_end' : dt_end,
            'platform' : platform,
            'instrument' : instrument
            }
    
    ncFile.close()

    return data

###############################################

def tdr_grid_variable(fname,Rec):
    """Read in a variable from a gridded NetCDF file 
    containing NOAA P-3 tail Doppler radar fields.
    
    Parameters::
    ----------
    fname : str
        Long path filename 
    VarName : str
        Variable name to access
        
    Usage::
    ----------
        VarOut = tdr_grid_variable(fname, Rec)
    """
    # Read the NetCDF
    ncFile = Dataset(fname,'r')

    # Get the variable of interest
    VarOut = ncFile.variables[VarName][:]
    
    ncFile.close()
    return VarOut
    
#######################################################
    
def read_lf_grid(fname):
    """Read in data from NetCDF file containing P3 flight level data created
    by NOAA AOC.  The NetCDF should be read in the main program and passed
    to this function.
        
    Parameters::
    ----------
    fname : str
        Long path filename [string]
        
    Output::
    ----------
    data : Dictionary of the following values
    
    metadata: Dictionary of metadata in file
    longitude : float
        Longitude array [degrees]
    latitude : float
        Latitude array [degrees]
    height : float
        Height array [km]
    fields: Dictionary of variables
        reflectivity : float
            Radar Reflectivity [dBZ]
    datetime_start : datetime object
        HHNN (Hour Minute) module start time
    datetime_end : datetime object
        HHNN (Hour Minute) module end time
    """
#=====================================
    # Read the NetCDF
    ncFile = Dataset(fname,'r')
    ncvars = ncFile.variables
#-------------------------------------
    # Grab the metadata stored in global attributes as a dictionary
    metadata = ncFile.__dict__
    
    # Find the module start information
    # Year
    yys = ncFile.variables['YYstart'][:]
    # Month	
    mms = ncFile.variables['MMstart'][:]
    # Day
    dds = ncFile.variables['DDstart'][:]
    # Hour
    hhs = ncFile.variables['HHstart'][:]
    # Minute
    nns = ncFile.variables['MIstart'][:]
    # Second
    sss = ncFile.variables['SSstart'][:]
    #--------------------------------
    # Find the module end information
    # Year
    yye = ncFile.variables['YYend'][:]
    # Month	
    mme = ncFile.variables['MMend'][:]
    # Day
    dde = ncFile.variables['DDend'][:]
    # Hour
    hhe = ncFile.variables['HHend'][:]
    # Minute
    nne = ncFile.variables['MIend'][:]
    # Second
    sse = ncFile.variables['SSend'][:]

    # Make an output string for Date-Time start/finish    	
#    dt_start = yys+'-'+mms+'-'+dds+' '+hhs+':'+nns+':'+sss 
#    dt_end = yye+'-'+mme+'-'+dde+' '+hhe+':'+nne+':'+sse  
    # Make a Datetime instance to pass
    utc = pytz.utc
    dt_start = datetime(yys,mms,dds,hhs,nns,sss,0)#,tzinfo=utc)
    dt_end = datetime(yye,mme,dde,hhe,nne,sse,0)#,tzinfo=utc)
    #---------------------------------------------------
    # Pull in the Lat / Lon information
    Lon = _ncvar_to_dict(ncvars['Lon'])
    Lat = _ncvar_to_dict(ncvars['Lat'])
    
    # Pull in the Height of aircraft
    Ht_level = _ncvar_to_dict(ncvars['Zlev'])
    
    # Pull in the Lat / Lon information
#    Lon = ncFile.variables['Lon'][:]
#    Lat = ncFile.variables['Lat'][:]
    
    # Pull in the Height of aircraft
#    Ht_level = ncFile.variables['Zlev'][:]
    #---------------------------------------------------
    fields = {}
    # Read the Reflectivity
    fields['dBZ'] = _ncvar_to_dict(ncvars['dBZ'])
    
#    badval = float(ncFile.MissingValue)
    
#    np.ma.masked_values(fields['dBZ']['data'], badval)
    
    # Pull out specific information on platform/instrument   
    try:
        ncFile.Platform
        platform = ncFile.title
    except:
        platform = 'p3'    

    try:
        ncFile.Instrument
        instrument = ncFile.instrument
    except:
        instrument = 'lf'
    
    # Create a dictionary to transfer the data
    data = {'metadata' : metadata,
            'longitude': Lon,
            'latitude': Lat,
            'height': Ht_level,
            'fields' : fields,
            'datetime_start' : dt_start,
            'datetime_end' : dt_end,
            'platform' : platform,
            'instrument' : instrument
            }
    
    ncFile.close()

    return data
    
#######################################################

def read_tdr_sweep(fname):
    """Read in data from NetCDF file containing P3 flight level data created
    by NOAA AOC.  The NetCDF should be read in the main program and passed
    to this function.
        
    Parameters::
    ----------
    fname : str
        Long path filename [string]
        
    Output::
    ----------
    data : Dictionary of the following values
    """
    data = pio.read_cfradial(fname)
    
    return data
    
####################
# Variable methods #
####################
def _ncvar_to_dict(ncvar):
    """ Convert a NetCDF Dataset variable to a dictionary. 
    Appropriated from PyArt package
    """
    d = dict((k, getattr(ncvar, k)) for k in ncvar.ncattrs())
    d['data'] = ncvar[:]
    if np.isscalar(d['data']):
        # netCDF4 1.1.0+ returns a scalar for 0-dim array, we always want
        # 1-dim+ arrays with a valid shape.
        d['data'] = np.array(d['data'])
        d['data'].shape = (1, )
    return d    
