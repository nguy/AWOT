"""
awot.io.read_uwka
=================

Scripts to read Wyoming Cloud Radar NetCDF data files.
 http://flights.uwyo.edu/wcr/

"""
# NOTES:: Testing was done on test data from 2014 flights.
#-------------------------------------------------------------------
# Load the needed packages
from netCDF4 import Dataset, num2date, date2num
import datetime
import numpy as np

from ..io.common import (_ncvar_subset_to_dict, _ncvar_subset_masked,
                         _ncvar_to_dict, _var_found, _var_not_found,
                         _get_epoch_units)

def read_wcr2(fname, field_mapping=None, file_mapping=None):
    '''
    Read Wyoming Cloud Radar 2 (WCR2) NetCDF level 2 data file.

    Parameters
    ----------
    fname : str
        Filename
    field_mapping : dict
        Mapping dictionary to use for field variable data.
        If None, then the default mapping is used.
    file_mapping : dict
            Mapping dictionary to use for pulling in data.
            If None, then the default mapping is used.

    Output
    ------
    data : dict
        AWOT dictionary instance.
        latitude : float
            Aircraft latitude [decimal degrees].
        longitude : float
            Aircraft longitude [decimal degrees].
        height : float
            Height of center of radar range gate [km].
	    altitude : float
    	    Aircraft altitude via GPS [km].
        tas : float
            Platform true airspeed [m/s].
        ground_speed : float
            Platform ground speed [m/s].
        dBZ_minimum : float
            Minimum detectable reflectivity at 1 km.
        aspect : float
            WCR range gate / (WCR time integer * TAS).
        beam_vector : float
            (East, North, Up)  beam vector unit.
        aircraft_wind : float
            In situ wind component at platform altitude along WCR beam.
            Positive is away from radar.
	    fields : Dictionary of variables in file
    	    dBZ : float
        	    Radar Equivalent Reflectivity Factor [dBZ].
	        velocity : float
    	        Mean Doppler radial velocity [m/s].
	        mask : int
    	        Target mask, see variable for notes.
        metadata : dict
            Dictionary of global attributes in file.
        project : str
            Project name.
        platform : str
            Platform name or identifier.
        flight_number : str
            Flight number or identifier.
    '''
    # Create a dictionary to hold data
    data = {}

    # Read the NetCDF
    ncFile = Dataset(fname,'r')
    ncvars = ncFile.variables

    # Grab the metadata stored in global attributes as a dictionary
    try:
        data['metadata'] = ncFile.__dict__
    except:
        data['metadata'] = None

    try:
        data['project'] = ncFile.ProjectName
    except:
        data['project'] = fname.split("_")[0]
    try:
        data['flight_number'] = ncFile.DataDate
    except:
        data['flight_number'] = fname.split("_")[2]

    # Set the platform
    try:
        data['platform'] = ncFile.Platform
    except:
        data['platform'] = 'n2uw'

    # Find the indices of not missing points
    Good = np.where(~np.isnan(ncFile.variables['time'][:]))

    data['time'] = _get_time(fname, ncFile, Good)

    # Grab the name map
    name_map_data = _get_wcr_data_map_level2()

    # Loop through the variables and pull data
    for varname in name_map_data:
        if name_map_data[varname] in ncvars:
##            data[varname] = _nc_var_masked(ncFile, name_map_data[varname], Good)
            data[varname] = _ncvar_to_dict(ncvars[name_map_data[varname]])
        else:
            data[varname] = None
            _var_not_found(varname)
    # Add fields to their own dictionary
    fields = {}

    # Grab a name map for WCR field data
    name_map_fields = _get_wcr_field_name_map()

    # Loop through the variables and pull data
    for varname in name_map_fields:
        if name_map_fields[varname] in ncvars:
            fields[varname] = _ncvar_subset_to_dict(ncvars[name_map_fields[varname]], Good)
        else:
            fields[varname] = None
            _var_not_found(varname)

    # Find the surface variable
    # See http://flights.uwyo.edu/uwka/wcr/projects/owles13/PROCESSED_DATA/WCR_L2_OWLES13.20131015.cdl
    # for details of mask properties
    surface = np.empty_like(data['latitude']['data'])
    condition = np.equal(fields['mask']['data'], 32)
    for nn in range(len(surface)):
         if np.any(condition[nn, :]):
             surface[nn] = data['height']['data'][np.where(condition[nn, :])[0][0]]
    data['surface'] = {'name' : "surface",
                       'long_name' : "Height of Surface",
                       'data' : surface,
                       'units' : 'meters'
                       }

    # Save to output dictionary
    data['fields'] = fields

    # Set the data format
    try:
        data['data_format'] = ncFile.WCR_BeamName
    except:
        data['data_format'] = 'wcr2'

    ncFile.close()

    return data

####################
##  _get methods  ##
####################

def _get_wcr_data_map_level2():
    '''Map WCR variable names to AWOT dictionary.'''
    name_map = {
               'latitude': 'LAT',
               'longitude': 'LON',
               'altitude': 'ALT',
               'tas': 'TAS',
               'ground_speed': 'GS',
               'aircraft_wind':  'acwcbeam',
               'aspect': 'wcraspect',
               'height': 'altrange',
               'height_ldr': 'ldrrange',
               'height_zdr': 'zdrrange',
               'reflectivity_minimum': 'min_reflectivity',
               'beam_vector': 'wcrbeamvector',
               }
    return name_map

def _get_wcr_field_name_map():
    '''Map WCR file names to AWOT dictionary.'''
    name_map = {
               'reflectivity': 'reflectivity',
               'velocity': 'velocity',
               'mask': 'wcrmask',
               'ldr': 'ldr',
               'zdr': 'zdr',
               }
    return name_map

def _get_time(fname, ncFile, Good_Indices):
    """Pull the time from WCR NetCDF file and convert to AWOT useable."""

    # Pull out the date, convert the date to a datetime friendly string

    # Now convert the time array into a datetime instance
    Time_unaware = num2date(ncFile.variables['time'][Good_Indices], _get_epoch_units())
    Time = {'data': Time_unaware, 'units': _get_epoch_units(),
            'title': 'Time', 'full_name': 'Time (UTC)'}
    return Time
