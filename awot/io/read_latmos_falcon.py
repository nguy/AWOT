"""
awot.io.read_latmos_falcon
=========================

This is a grouping of scripts designed to process (French) 
 Falcon distributed by LATMOS.

Note: This has only been tested with DYNAMO data files, versions
may change and another function may be needed.
 
"""
# Load the needed packages
from netCDF4 import Dataset, num2date, date2num
import datetime
import numpy as np
import pytz

def flight_data(fname):
    '''A wrapper to call the rasta_dynamic.'''
    data = rasta_dynamic(fname)
    return data

def rasta_wind(fname):
    '''A wrapper to call the rasta_dynamic module.'''
    data = rasta_dynamic(fname)
    return data

def rasta_radar(fname):
    '''A wrapper to call the rasta_dynamic module.'''
    data = rasta_dynamic(fname)
    return data

def rasta_radonvar(fname):
    '''A wrapper to call the rasta_microphysics module.'''
    data = rasta_microphysics(fname)
    return data

def rasta_dynamic(fname):
    '''
    Read in NetCDF data file containing Falcon dynamic properties 
    retrieved from radar measurements.

    Parameters
    ----------
    fname : string
        Filename [string]

    Output
    ------
    data : Dictionary of the following values

	    latitude : float
    	    Aircraft latitude [decimal degrees]
	    longitude : float
    	    Aircraft longitude [decimal degrees]
	    altitude : float
    	    Aircraft altitude via GPS [km]
	    height : float
    	    Height [km]
	    temperature : float
    	    Temperature at aircraft altitude [degrees C]
	    pressure: float
    	    Pressure at aircraft altitude [hPa]
	    Vx_flight_level : float
    	    Vx at flight level from in situ measurements [m/s]
	    Vy_flight_level : float
    	    Vy at flight level from in situ measurements [m/s]
	    Vz_flight_level : float
    	    Vz at flight level from in situ measurements [m/s]
	    zonal_wind : float
    	    Zonal wind component [m/s]
	    meridional_wind : float
    	    Meriodional wind component [m/s]
	    fields : Dictionary of variables in file
    	    reflectivity : float
        	    Radar Reflectivity [dBZ]
	        Uwind : float
    	        Wind along aircraft longitudinal axis wind [m/s]
	        Vwind : float
    	        Wind perpendicular to aircraft longitudinal axis wind [m/s]
        	Wwind : float
	            Vertical wind component [m/s]
    	    term_fall_speed : float
        	    Terminal fall speed [m/s]
	        term_fall_speed_weighted : float
    	        Terminal fall speed weighted by Vt-Z [m/s]
        metadata : Dictionary of global attributes in file
        project : str
            Project Name
        platform : str
            Platform name or identifier
        flight_number : str
            Flight number or identifer
    '''
    # Read the NetCDF
    ncFile = Dataset(fname,'r')
    ncvars = ncFile.variables

    # Create dictionary to hold output data
    data = {}
    
    # Grab the metadata stored in global attributes as a dictionary
    try:
        data['metadata'] = ncFile.__dict__
    except:
        data['metadata'] = None

    # Find the indices of not missing points
    Good = np.where(~np.isnan(ncFile.variables['time'][:]))

    # Create a dictionary to hold data
    data['time'] = _get_latmos_time(fname, ncFile, Good)
    
    # Grab a name map for RASTA dynamic file data
    name_map_data = _get_dynamic_flight_name_map()

    # Loop through the variables and pull data
    for varname in name_map_data:
        if varname in ncvars:
            data[varname] = _nc_var_masked(ncFile, name_map_data[varname], Good)
            if varname is 'altitude':
                data[varname] = data[varname] * 1000.
        else:
            data[varname] = None

    try:
        Ht = _nc_height_var_to_dict(ncvars['height'])
    except:
        Ht = None
    data['height'] = Ht

    try:
        mask_hydro = _nc_radar_var_to_dict(ncvars['Mask'], Good)
    except:
        mask_hydro = None
    data['mask_hydro'] = mask_hydro

    # Add fields to their own dictionary
    fields = {}

    # Grab a name map for RASTA dynamic field data
    name_map_fields = _get_dynamic_field_name_map()

    # Loop through the variables and pull data
    for varname in name_map_fields:
        try:
            fields[varname] = _ncvar_radar_to_dict(ncvars[name_map_fields[varname]], Good)
        except:
            fields[varname] = None

    # Save to output dictionary
    data['fields'] = fields

    # Pull out global attributes
    try:
        ncFile.ProjectName
        data['project'] = ncFile.ProjectName
    except:
        data['project'] = fname.split("_")[0]
    try:
        ncFile.FlightNumber
        data['flight_number'] = ncFile.FlightNumber
    except:
        data['flight_number'] = fname.split("_")[2]

    # Create a dictionary to transfer the data
    data['platform'] = 'falcon'
    data['data_format'] = 'falcon_vertical'

    ncFile.close()
    return data


def rasta_microphysics(fname):
    '''
    Read in NetCDF data file containing Falcon microphysical properties.

    Parameters
    ----------
    fname : string
        Filename [string]

    Output
    ------
    data : Dictionary of the following values
        latitude : float
            Aircraft latitude [decimal degrees]
        longitude : float
            Aircraft longitude [decimal degrees]
        height : float
            Height [km]
        fields : Dictionary of variables in file
            extinction : float
                Visible extinction [1/m]
            time : float
                Aircraft time array
            n0start : float
                Normalized concentration parameter [m^-4]
            iwc : float
                Ice water content [kg m^-3]
            effective_radius : float
                Effective radius [micron]
            Dm : float
                Equivalen volume weighted diameter [micron]
            Nt : float
                Number_concentration [m^-3]
            dBZ : float
                Radar reflectivity [dBZ]
            Z_fwd : float
                Forward-modelled radar reflectivity [dBZ]
            term_velocity : float
                Terminal fall velocity [m s^-1]
            term_velocity_fwd : float
                Forward-modelled terminal fall velocity [m s^-1]
            temperature : float
                Temperature [degrees C]
            aM : float
                Retrieved aM (M (D)=aD^b) [unitless]
            metadata : Dictionary of global attributes in file
            project : str
                Project Name
            platform : str
                Platform name or identifier
            flight_number : str
                Flight number or identifer
    '''

    # Read the NetCDF
    ncFile = Dataset(fname,'r')
    ncvars = ncFile.variables

    # Create dictionary to hold output data
    data = {}
    
    # Grab the metadata stored in global attributes as a dictionary
    try:
        data['metadata'] = ncFile.__dict__
    except:
        data['metadata'] = None

    # Find the indices of not missing points
    Good = np.where(~np.isnan(ncFile.variables['time'][:]))

    # Create a dictionary to hold data
    data['time'] = _get_latmos_time(fname, ncFile, Good)
    
    # Pull out each variable
    try:
        Lat = ncFile.variables['latitude'][Good]
        data['latitude'] = Lat
    except:
        data['latitude'] = None
    try:
        Lon = ncFile.variables['longitude'][Good]
        data['longitude'] = Lon
    except:
        data['longitude'] = None
    try:
        Ht = _nc_height_var_to_dict(ncvars['height'])
        data['height'] = Ht
    except:
        data['height'] = None

    # Add fields to their own dictionary
    fields = {}

    # Grab a name map for RASTA microphysics field data
    name_map = _get_microphysics_field_name_map()

    # Loop through the variables and pull data
    for varname in name_map:
        try:
            fields[varname] = _ncvar_radar_to_dict(ncvars[name_map[varname]], Good)
        except:
            fields[varname] = None

    if fields['Dm'] is not None:
        fields['Dm']['data'][:] = fields['Dm']['data'][:] * 1000.
        fields['Dm']['units'] = 'mm'
    data['fields'] = fields

    # Pull out global attributes
    try:
        ncFile.ProjectName
        data['project'] = ncFile.ProjectName
    except:
        data['project'] = fname.split("_")[0]
    try:
        ncFile.FlightNumber
        data['flight_number'] = ncFile.FlightNumber
    except:
        data['flight_number'] = fname.split("_")[2]

    # Now mask missing values
    if Lat is not None:
        np.ma.masked_invalid(Lat)
    if Lon is not None:
        np.ma.masked_invalid(Lon)

    # Create a dictionary to transfer the data
    data['platform'] = 'falcon'
    data['data_format'] = 'falcon_vertical'

    ncFile.close()
    return data

###############################
##  Create Variable methods  ##
###############################

def _get_dynamic_flight_name_map():
    '''Map RASTA dynamic variables to AWOT structure.'''
    name_map = {
               'latitude': 'latitude',
               'longitude': 'longitude',
               'height': 'height',
               'altitude': 'altitude',
               'temperature': 'temperature',
               'pressure': 'pressure',
               'Vx_flight_level': 'Vx_insitu',
               'Vy_flight_level': 'Vy_insitu',
               'Vz_flight_level': 'Vz_insitu',
               'zonal_wind':  'VE',
               'meridional_wind': 'VN',
               }
    return name_map

def _get_microphysics_field_name_map():
    '''Map RASTA microphysics variables to AWOT structure..'''
    name_map = {
               'extinction': 'extinction',
               'n0star': 'n0star',
               'iwc': 'iwc',
               'effective_radius': 'effective_radius',
               'Dm': 'Dm',
               'Nt': 'Nt',
               'reflectivity': 'Z',
               'Z_fwd': 'Z_fwd',
               'term_velocity': 'vt',
               'term_velocity_fwd': 'vt_fwd',
               'temperature': 'T',
               'aM': 'aM',
               }
    return name_map

def _get_dynamic_field_name_map():
    '''
    Map RASTA dynamic variables to AWOT structure.
        'Z'          : Radar reflectivity
        'Vx'         : U wind (Along aircraft longitudinal axis)
        'Vy'         : V wind (Perpendicular to aircraft longitudinal axis)
        'Vz'         : Vertical wind
        'Vt'         : Terminal velocity
        'Vt_weighted : Weighted Terminal velocity
    '''
    name_map = {
               'reflectivity': 'Z',
               'Uwind': 'Vx',
               'Vwind': 'Vy',
               'Wwind': 'Vz',
               'term_fall_speed': 'Vt',
               'term_fall_speed_weighted': 'Vt_weighted'
               }
    return name_map
    
def _get_latmos_time(fname, ncFile, Good_Indices):
    """Pull the time from RASTA file and convert to AWOT useable."""
    # Pull out the date, convert the date to a datetime friendly string
    # Adds dashes between year, month, and day
    # This assumes that the date is the 2nd instance in the filename!!
    yyyymmdd = fname.split("_")[1]
    StartDate = yyyymmdd[0:4] + '-' + yyyymmdd[4:6] + '-' + yyyymmdd[6:8]

    # Create the time array
    # First find the good indices, there can be missing values 
    # in Falcon Data (seriously come on!)
    # So we have remove these indices from all fields in the future

    # Now convert the time array into a datetime instance
    dtHrs = num2date(ncFile.variables['time'][Good_Indices], 'hours since ' + StartDate + '00:00:00+0:00')
    # Now convert this datetime instance into a number of seconds since Epoch
    TimeSec = date2num(dtHrs, 'seconds since 1970-01-01 00:00:00+0:00')
    # Now once again convert this data into a datetime instance
    Time_unaware = num2date(TimeSec, 'seconds since 1970-01-01 00:00:00+0:00')
    Time = Time_unaware#.replace(tzinfo=pytz.UTC)
    return Time

def _nc_var_masked(ndFile, ncvar, Good_Indices):
    """Convert a NetCDF variable into a masked variable."""
    d = ncFile.variables[ncvar][Good_Indices]
    np.ma.masked_invalid(d)
    return d

def _nc_radar_var_to_dict(ncvar, Good_Indices):
    """
    Convert a NetCDF Dataset variable to a dictionary. 
    Appropriated from PyArt package.
    """
    d = dict((k, getattr(ncvar, k)) for k in ncvar.ncattrs())
    d['data'] = ncvar[:]
    if np.isscalar(d['data']):
        # netCDF4 1.1.0+ returns a scalar for 0-dim array, we always want
        # 1-dim+ arrays with a valid shape.
        d['data'] = np.array(d['data'][Good_Indices, :])
        d['data'].shape = (1, )
    return d

def _nc_height_var_to_dict(ncvar):
    """
    Convert a NetCDF Dataset variable to a dictionary. 
    Appropriated from PyArt package.
    """
    d = dict((k, getattr(ncvar, k)) for k in ncvar.ncattrs())
    d['data'] = ncvar[:]
    if np.isscalar(d['data']):
        # netCDF4 1.1.0+ returns a scalar for 0-dim array, we always want
        # 1-dim+ arrays with a valid shape.
        d['data'] = np.array(d['data'][:])
        d['data'].shape = (1, )

    d['data'][:] = Ht['data'][:] * 1000.
    d['units'] = 'meters'

    return d