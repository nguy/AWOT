"""
awot.io.read_nasa_hiwrap
=========================

Grouping of scripts designed to NASA HIWRAP data files.
http://har.gsfc.nasa.gov/index.php?section=13

"""
# Load the needed packages
from netCDF4 import Dataset, num2date, date2num
import datetime
import numpy as np

from ..io.common import (_ncvar_subset_to_dict, _ncvar_subset_masked,
                         _ncvar_to_dict, _var_not_found,
                         _get_epoch_units)


def read_hiwrap_netcdf(fname, mapping_dict=None, field_mapping=None):
    '''
    Read in NetCDF data file containing HIWRAP radar and flight data.

    Parameters
    ----------
    fname : string
        Filename [string].
    mapping_dict : dict
        Mapping dictionary to use for flight variable data.
        If None, then the default mapping is used.
    field_mapping : dict
        Mapping dictionary to use for field variable data.
        If None, then the default mapping is used.

    Output
    ------
    data : dict
        AWOT dictionary instance.
        latitude : float
            Aircraft latitude [decimal degrees].
        longitude : float
            Aircraft longitude [decimal degrees].
        altitude : float
            Aircraft altitude via GPS [km].
        height : float
            Height [km].
        temperature : float
            Temperature at aircraft altitude [degrees C].
        pressure: float
            Pressure at aircraft altitude [hPa].
        Vx_flight_level : float
            Vx at flight level from in situ measurements [m/s].
        Vy_flight_level : float
            Vy at flight level from in situ measurements [m/s].
        Vz_flight_level : float
            Vz at flight level from in situ measurements [m/s].
        zonal_wind : float
            Zonal wind component [m/s].
        meridional_wind : float
            Meriodional wind component [m/s].
        fields : dict
            reflectivity : float
                Radar Reflectivity [dBZ].
            Uwind : float
                Wind along aircraft longitudinal axis wind [m/s].
            Vwind : float
                Wind perpendicular to aircraft longitudinal axis wind [m/s].
            Wwind : float
                Vertical wind component [m/s].
            term_fall_speed : float
                Terminal fall speed [m/s].
            term_fall_speed_weighted : float
                Terminal fall speed weighted by Vt-Z [m/s].
        metadata : dict
            Dictionary of global attributes in file.
        project : str
            Project Name.
        platform : str
            Platform name or identifier.
        flight_number : str
            Flight number or identifier.
    '''
    # Read the NetCDF
    ncFile = Dataset(fname, 'r')
    ncvars = ncFile.variables

    # Create dictionary to hold output data
    data = {}

    # Grab the metadata stored in global attributes as a dictionary
    try:
        data['metadata'] = ncFile.__dict__
    except:
        data['metadata'] = None

    # Get missing data value
    baddata = ncFile.variables['missing'][:]
    # Find the indices of not missing points
    Good = np.where(~np.isnan(ncFile.variables['time'][:]))

    # Create a dictionary to hold data
    data['time'] = _get_hiwrap_time(fname, ncFile, Good)

    # Grab a name map for HIWRAP flight data
    if mapping_dict is None:
        name_map_flightdata = _get_hiwrap_flight_name_map()
    else:
        name_map_flightdata = mapping_dict

    # Loop through the variables and pull data
    # Adjust altitude to meters from kilometers
    for varname in name_map_flightdata:
        if name_map_flightdata[varname] in ncvars:
            data[varname] = _ncvar_to_dict(
                ncvars[name_map_flightdata[varname]])
            if varname is 'altitude':
                data[varname]['data'] = data[varname]['data'] * 1000.
                data[varname]['units'] = 'meters'
        else:
            data[varname] = None
            _var_not_found(varname)

    # Add fields to their own dictionary
    fields = {}

    # Grab a name map for HIWRAP field data
    if field_mapping is None:
        name_map_fields = _get_hiwrap_field_name_map()
    else:
        name_map_fields = field_mapping

    # Loop through the variables and pull data
    for varname in name_map_fields:
        try:
            fields[varname] = _ncvar_to_dict(
                ncvars[name_map_fields[varname]])
            # Apply mask to any points with missing value
            # indicated by file
            fields[varname]['data'] = np.ma.masked_equal(fields[varname]['data'], baddata)
        except:
            fields[varname] = None
            _var_not_found(varname)

    # Save to output dictionary
    data['fields'] = fields

    # Pull out global attributes
    try:
        data['project'] = ncFile.experiment
    except:
        data['project'] = fname.split("_")[0]
    try:
        data['flight_number'] = ncFile.FlightNumber
    except:
        data['flight_number'] = None

    # Include HIWRAP specific metadata if it exists
    try:
        data['metadata']['year'] = _ncvar_to_dict(ncvars['year'])
    except:
        data['metadata']['year'] = None
    try:
        data['metadata']['radar_frequency'] = _ncvar_to_dict(ncvars['freq'])
    except:
        data['metadata']['radar_frequency'] = None
    try:
        data['metadata']['incidence_angle'] = _ncvar_to_dict(ncvars['incid'])
    except:
        data['metadata']['incidence_angle'] = None
    try:
        data['metadata']['range_gate_spacing'] = _ncvar_to_dict(ncvars['gatesp'])
    except:
        data['metadata']['range_gate_spacing'] = None

    # Create a dictionary to transfer the data
    data['platform'] = 'hiwrap'
    data['data_format'] = 'hiwrap_vertical'

    ncFile.close()
    return data

###############################
#   Create Variable methods   #
###############################


def _get_hiwrap_flight_name_map():
    '''Map HIWRAP NetCDF flight variables to AWOT structure.'''
    name_map = {
               'latitude': 'lat',
               'longitude': 'lon',
               'altitude': 'alt',
               'roll': 'roll',
               'pitch': 'pitch',
               'tilt': 'tilt',
               'heading': 'head',
               'rotation': 'rot',
               'track': 'track',
               'range': 'range',
               'east_ground_speed': 'evel',
               'north_ground_speed': 'nvel',
               'vert_speed': 'wvel',
               'vert_aircraft_velocity': 'vacft'
               }
    return name_map


def _get_hiwrap_field_name_map():
    '''Map HIWRAP NetCDF field variables to AWOT structure.'''
    name_map = {
               'surface_gate': 'sgate',
               'power': 'pwr',
               'reflectivity': 'ref',
               'velocity': 'dopcorr'
               }
    return name_map

def _get_hiwrap_time(fname, ncFile, Good_Indices):
    """Pull the time from HIWRAP file and convert to AWOT useable."""
    # The time structure is odd here (to me) and is in
    # seconds since last Sunday - wtf

    # Pull out the date, convert the date to a datetime friendly string
    # Adds dashes between year, month, and day
    # This assumes that the date is the 4th instance in the filename!!
    yyyymmdd = fname.split("_")[3]

    # Find the date for Sunday previous and check this (should be = 6 for Sunday)
    startday = int(yyyymmdd[6:8]) - int(divmod(ncFile.variables['time'][0], 24 * 3600)[0])
    if datetime.date(ncFile.variables['year'][:], int(yyyymmdd[4:6]), startday).weekday() != 6:
        print("Time could be incorrect, check file to see if time units "
              "are 'computer time (sec from last Sunday at 12 am)'")

    StartDate = yyyymmdd[0:4] + '-' + yyyymmdd[4:6] + '-' + str(startday)

    # Create the time array
    # First find the good indices, there can be missing values
    # in Falcon Data (seriously come on!)
    # So we have remove these indices from all fields in the future

    # Now convert the time array into a datetime instance
    dtHrs = num2date(ncFile.variables['time'][
        Good_Indices], 'seconds since ' + StartDate + '00:00:00+0:00')
    # Now convert this datetime instance into a number of seconds since Epoch
    TimeSec = date2num(dtHrs, _get_epoch_units())
    # Now once again convert this data into a datetime instance
    Time_unaware = num2date(TimeSec, _get_epoch_units())
    Time = {'data': Time_unaware, 'units': _get_epoch_units(),
            'title': 'Time', 'full_name': 'Time (UTC)'}

    return Time
