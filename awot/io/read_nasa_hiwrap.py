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

try:
    import h5py
    _H5PY_AVAILABLE = True
except ImportError:
    _H5PY_AVAILABLE = False

from ..io.common import (_ncvar_subset_to_dict, _ncvar_subset_masked,
                         _ncvar_to_dict, _h5var_to_dict,
                         _var_not_found,_get_epoch_units)


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
    data : AWOT dictionary instance.
        latitude : float
            Aircraft latitude [decimal degrees].
        longitude : float
            Aircraft longitude [decimal degrees].
        altitude : float
            Aircraft altitude via GPS [km].
        roll : float
            Aircraft roll angle [degrees].
        pitch : float
            Aircraft pitch angle [degrees].
        tilt : float
            Aircraft tilt angle [degrees].
        heading : float
            Aircraft heading angle [degrees (from north)].
        rotation : float
            Aircraft rotation angle [degrees].
        track : float
            Aircraft track angle [degrees].
        range: float
            Radar range from aircraft [meters].
        east_ground_speed : float
            Aircraft eastward velocity [m/s].
        north_ground_speed : float
            Aircraft northward velocity [m/s].
        vert_aircraft_velocity : float
            Vertical aircraft speed [m/s].
        vert_speed : float
            Vertical wind compenent [m/s].
        surface_gate : float
            Index of surface gate in profile [unitless].
        height : float
            Height of radar observations [m].
        topo : float
            Height of surface topography [m].
        fields : dict
            reflectivity : float
                Radar Reflectivity [dBZ].
            velocity : float
                Doppler velocity, corrected [m/s].
            power : float
                Radar return power[dB].
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
    data['time'] = _get_old_hiwrap_time(fname, ncFile, Good)

    # Grab a name map for HIWRAP flight data
    if mapping_dict is None:
        name_map_flightdata = _get_old_hiwrap_flight_name_map()
    else:
        name_map_flightdata = mapping_dict

    # Loop through the variables and pull data
    # Adjust altitude to meters from kilometers
    for varname in name_map_flightdata:
        if name_map_flightdata[varname] in ncvars:
            data[varname] = _ncvar_to_dict(
                ncvars[name_map_flightdata[varname]])
## NG TEST FILE INCORRECT UNITS - MAY NEED TO REINSTATE
#            if data[varname]['units'] == 'km':
#                print(varname)
#                data[varname]['data'] = data[varname]['data'] * 1000.
#                data[varname]['units'] = 'meters'
        else:
            data[varname] = None
            _var_not_found(varname)

    # Replace negative range gates - used for calibration purposes
    # This likely has no affect given data tested from GRIP campaign
    gate_mask = np.ma.less(ncvars['range'][:], 0.)
    data['range']['data'][gate_mask] = np.nan

    # Calculate the height field by subtracting the range from altitude
    r2D, alt2D = np.meshgrid(data['range']['data'], data['altitude']['data'])
    data['height'] = {'name': "Height",
                      'long_name': "Height above surface",
                      'data': (alt2D - r2D),
                      'units': data['altitude']['units']}

    # Calculate the topo height
    ## NG NEED TO CHECK TO MAKE SURE CORRECT
    sfcgateheight = np.array(
       [data['range']['data'][int(s)] for s in ncvars['sgate'][:]])
    topo = data['altitude']['data'][:] - sfcgateheight[:]
    data['topo'] = {'name' : "topo",
                    'long_name' : "Height of Topography",
                    'data' : topo,
                    'units' : 'meters'
                    }

    # Add fields to their own dictionary
    fields = {}

    # Grab a name map for HIWRAP field data
    if field_mapping is None:
        name_map_fields = _get_old_hiwrap_field_name_map()
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
            fields[varname]['data'][:, gate_mask] = np.nan
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

def read_hiwrap_h5(fname, mapping_dict=None, field_mapping=None):
    '''
    Read in HDF5 data file containing HIWRAP radar and flight data.

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
               'profile_start_range': 'rng0',
               'surface_doppler_velocity': 'surfvel',
               'antenna_elevation_angle': 'antElevation',
    data : AWOT dictionary instance.
        latitude : float
            Aircraft latitude [decimal degrees].
        longitude : float
            Aircraft longitude [decimal degrees].
        altitude : float
            Aircraft altitude via GPS [km].
        roll : float
            Aircraft roll angle [degrees].
        pitch : float
            Aircraft pitch angle [degrees].
        tilt : float
            Aircraft tilt angle [degrees].
        heading : float
            Aircraft heading angle [degrees (from north)].
        track : float
            Aircraft track angle [degrees].
        range: float
            Radar range from aircraft [meters].
        east_ground_speed : float
            Aircraft eastward velocity [m/s].
        north_ground_speed : float
            Aircraft northward velocity [m/s].
        vert_aircraft_velocity : float
            Vertical aircraft speed [m/s].
        vert_speed : float
            Vertical wind compenent [m/s].
        surface_gate : float
            Index of surface gate in profile [unitless].
        profile_start_range : Index
            Index of of the profile start
        height : float
            Height of radar observations [m].
        topo : float
            Height of surface topography [m].
        fields : dict
            reflectivity : float
                Radar Reflectivity [dBZ].
            velocity : float
                Doppler velocity, corrected [m/s].
            power : float
                Radar return power[dB].
        metadata : dict
            Dictionary of global attributes in file.
        project : str
            Project Name.
        platform : str
            Platform name or identifier.
        flight_number : str
            Flight number or identifier.
    '''

    # check that h5py is available
    if not _H5PY_AVAILABLE:
        raise MissingOptionalDependency(
            "h5py is required to use read_odim_h5 but is not installed")
    # Read the NetCDF
    h5File = h5py.File(fname, 'r')

    # Create dictionary to hold output data
    data = {}

    # Grab the metadata stored in global attributes as a dictionary
    try:
        data['metadata'] = h5File.attrs()
    except:
        data['metadata'] = {'global': None}

    # Get missing data value
    baddata = np.array(h5File['missing'])
    # Find the indices of not missing points
#    Good = np.where(~np.isnan(h5File['time'][:]))

    # Create a dictionary to hold data
    data['time'] = _get_hiwrap_time(h5File)

    # Grab a name map for HIWRAP flight data
    if mapping_dict is None:
        name_map_flightdata = _get_hiwrap_flight_name_map()
    else:
        name_map_flightdata = mapping_dict

    # Loop through the variables and pull data
    # Adjust altitude to meters from kilometers
    for varname in name_map_flightdata:
        if name_map_flightdata[varname] in h5File.keys():
            data[varname] = _h5var_to_dict(
                h5File[name_map_flightdata[varname]])
        else:
            data[varname] = None
            _var_not_found(varname)
    # Replace negative range gates - used for calibration purposes
    gate_mask = np.ma.less(h5File['rangevec'][:], 0.)
    data['range']['data'][gate_mask] = np.nan

    # Add the height field, need to calculate
    r2D, alt2D = np.meshgrid(data['range']['data'], data['altitude']['data'])
    data['height'] = {'name': "Height",
                      'long_name': "Height above surface",
                      'data': (alt2D - r2D),
                      'units': data['altitude']['units']}

    # Calculate the topo height
    ## NG NEED TO CHECK TO MAKE SURE CORRECT
    sfcgaterange = np.array(
       [data['range']['data'][int(s)] for s in h5File['sgate'][:]])
    topo = data['altitude']['data'][:] - sfcgaterange[:]
    data['topo'] = {'name' : "topo",
                    'long_name' : "Height of Topography",
                    'data' : topo,
                    'units' : 'meters'
                    }

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
            fields[varname] = _h5var_to_dict(
                h5File[name_map_fields[varname]])
            # Apply mask to any points with missing value
            # indicated by file
            fields[varname]['data'] = fields[varname]['data'].T
            fields[varname]['data'] = np.ma.masked_equal(fields[varname]['data'], baddata)
            fields[varname]['data'][:, gate_mask] = np.nan
        except:
            fields[varname] = None
            _var_not_found(varname)
    # Save to output dictionary
    data['fields'] = fields

    # Pull out global attributes
    try:
        data['project'] = np.array(h5File['ExperName'])[0]
    except:
        data['project'] = fname.split("_")[0]
    try:
        data['flight_number'] = np.array(h5File['FlightNumber'])[0]
    except:
        data['flight_number'] = None
    try:
        data['radar_name'] = np.array(h5File['radarName'])[0]
    except:
        data['radar_name'] = None

    # Include HIWRAP specific metadata if it exists
    try:
        data['metadata']['year'] = _h5var_to_dict(h5File['utcYear'])
    except:
        data['metadata']['year'] = None
    try:
        data['metadata']['radar_frequency'] = _h5var_to_dict(h5File['Frequency'],
                                                 units="GHz", standard_name="Radar Frequency")
    except:
        data['metadata']['radar_frequency'] = None
    try:
        data['metadata']['radar_wavelength'] = _h5var_to_dict(h5File['Wavelength'],
                                                 units="meters", standard_name="Radar Wavelength")
    except:
        data['metadata']['radar_wavelength'] = None
    try:
        data['metadata']['incidence_angle'] = _h5var_to_dict(h5File['incid'],
                                                 units="degrees", standard_name="Incidence Angle")
    except:
        data['metadata']['incidence_angle'] = None
    try:
        data['metadata']['range_gate_spacing'] = _h5var_to_dict(h5File['gatesp'],
                                                 units="meters", standard_name="Gate spacing")
    except:
        data['metadata']['range_gate_spacing'] = None
    try:
        data['metadata']['antenna_elevation_angle'] = _h5var_to_dict(h5File['antElevation'],
                                                 units="degrees", standard_name="Antenna Elevation")
    except:
        data['metadata']['antenna_elevation_angle'] = None
    try:
        data['metadata']['range_gate_spacing'] = _h5var_to_dict(h5File['nbeams'],
                                                 standard_name="Number of beams")
    except:
        data['metadata']['range_gate_spacing'] = None
    try:
        data['metadata']['dual_nyquist_velocity'] = _h5var_to_dict(h5File['vnyqDual'],
                                                 units="m/s", standard_name="Nyquist Velocity",
                                                 long_name="Dual-frequency Nyquist Velocity")
    except:
        data['metadata']['dual_nyquist_velocity'] = None
#    try:
#        data['metadata']['creation_date'] = np.array(h5File['creationDate'])[0]
#    except:
#        data['metadata']['creation_date'] = None

    data['platform'] = 'hiwrap'
    data['data_format'] = 'hiwrap_vertical'

    h5File.close()
    return data
###############################
#   Create Variable methods   #
###############################

def _get_old_hiwrap_flight_name_map():
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
               'vert_aircraft_velocity': 'vacft',
               'vert_speed': 'wvel',
               'surface_gate': 'sgate',
               }
    return name_map


def _get_old_hiwrap_field_name_map():
    '''Map HIWRAP NetCDF field variables to AWOT structure.'''
    name_map = {
               'power': 'pwr',
               'reflectivity': 'ref',
               'velocity': 'dopcorr'
               }
    return name_map

def _get_hiwrap_flight_name_map():
    '''Map HIWRAP HDF flight variables to AWOT structure.'''
    name_map = {
               'latitude': 'lat',
               'longitude': 'lon',
               'altitude': 'height',
               'roll': 'roll',
               'pitch': 'pitch',
#               'tilt': 'tilt',
               'heading': 'head',
#               'rotation': 'rot',
               'track': 'track',
               'range': 'rangevec',
               'profile_start_range': 'rng0',
               'east_ground_speed': 'evel',
               'north_ground_speed': 'nvel',
               'vert_aircraft_velocity': 'vacft',
               'vert_speed': 'wvel',
               'surface_doppler_velocity': 'surfvel',
               'surface_gate': 'sgate',
               'antenna_elevation_angle': 'antElevation',
               }
    return name_map


def _get_hiwrap_field_name_map():
    '''Map HIWRAP HDF field variables to AWOT structure.'''
    name_map = {
               'power': 'stitchedPower',
               'reflectivity': 'stitchedReflectivity',
               'velocity': 'stitchedVelocity',
               'surface_power': 'surfpwr_int',
               'surface_velocity': 'surfvel',
               }
    return name_map

def _get_old_hiwrap_time(fname, ncFile, Good_Indices):
    """
    Pull the time from HIWRAP file and convert to AWOT useable.
    The time structure is odd here (to me) and is in
    seconds since last Sunday.
    The assumption that the data is the 4th 'field' in the filename
    is required to make this work.
    """
    # Pull out the date, convert the date to a datetime friendly string
    # Adds dashes between year, month, and day
    yyyymmdd = fname.split("_")[3]

    # Find the date for Sunday previous and check this (should be = 6 for Sunday)
    startday = int(yyyymmdd[6:8]) - int(divmod(ncFile.variables['time'][0], 24 * 3600)[0])
    if datetime.date(ncFile.variables['year'][:], int(yyyymmdd[4:6]), startday).weekday() != 6:
        print("Time could be incorrect, check file to see if time units "
              "are 'computer time (sec from last Sunday at 12 am)'")

    StartDate = yyyymmdd[0:4] + '-' + yyyymmdd[4:6] + '-' + str(startday)

    # Create the time array
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

def _get_hiwrap_time(h5File):
    """
    Pull the time from the HIWRAP file and convert to AWOT useable.
    This method assumes that the first values for time variables
    are indeed correct.
    """
    year = np.array(h5File['utcYear'])
    month = np.array(h5File['utcMonth'])
    day = np.array(h5File['utcDay'])
    seconds = np.array(h5File['timeUTC']) * 3600.
    t_initial = ("%s-%s-%s 00:00:0.00Z"%(str(year[0]), str(month[0]), str(day[0])))

    dtHrs = num2date(seconds, 'seconds since %s'%t_initial)
    # Now convert this datetime instance into a number of seconds since Epoch
    TimeSec = date2num(dtHrs, _get_epoch_units())
    # Now once again convert this data into a datetime instance
    Time_unaware = num2date(TimeSec, _get_epoch_units())
    Time = {'data': Time_unaware, 'units': _get_epoch_units(),
            'title': 'Time', 'full_name': 'Time (UTC)'}
    return Time
