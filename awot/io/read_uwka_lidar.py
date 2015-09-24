"""
awot.io.read_uwka_lidar
=========================

Scripts for reading NetCDF data files from the 
Wyoming Cloud Lidar, part of the University of Wyoming
King Air Research Aircraft facility.
 http://flights.uwyo.edu/wcl/
"""
# Load the needed packages
from __future__ import print_function
from netCDF4 import Dataset, num2date, date2num
import datetime
import numpy as np

from ..io.common import (_ncvar_subset_to_dict, _ncvar_subset_masked,
                         _var_found, _var_not_found)

def read_wcl(fname):
    '''
    Read in NetCDF data file containing Wyoming Cloud Lidar data.

    Parameters
    ----------
    fname : str
        Filename [string]

    Output
    ------
    data : dict
        AWOT dictionary instance.
        fields : dict
            Dictionary containing fields from WCL file.
        metadata : dict
            Dictionary of global attributes in file.
        project : str
            Project name.
        platform : str
            Platform name or identifier.
        flight_number : str
            Flight number or identifer.
        data_format: str
            AWOT identifying string.
    '''
    # Create a dictionary to hold data
    data = {}

    # Read the NetCDF
    ncFile = Dataset(fname, 'r')
    ncvars = ncFile.variables

    # Grab the metadata stored in global attributes as a dictionary
    # Pull out other global attribute information
    data['metadata'] = ncFile.__dict__
    try:
        data['project'] = ncFile.ProjectName
    except:
        project = fname.split("_")[0]
    try:
        data['flight_number'] = ncFile.FlightNumber
    except:
        flightnum = fname.split("_")[2]

    # Set the platform and data format keys
    data['platform'] = ncFile.Platform
    data['data_format'] = 'wcl'

    # Find the indices of not missing points
    Good = np.where(~np.isnan(ncvars['time'][:]))

    Time = _get_time(fname, ncFile, Good)
    # Put time into output dictionary
    data['time'] = Time

    # Grab a name map for WCL data
    name_map_data = _get_wcl_name_map()

    # Loop through the variables and pull data
    for varname in name_map_data:
        if name_map_data[varname] in ncvars:
            data[varname] = _ncvar_subset_masked(
                ncFile, name_map_data[varname], Good)
        else:
            data[varname] = None
            _var_not_found(name_map_data[varname])

    # Add fields to their own dictionary
    fields = {}

    # Grab a name map for WCR field data
    name_map_fields = _get_wcl_field_name_map()

    # Loop through the variables and pull data
    for varname in name_map_fields:
        if name_map_fields[varname] in ncvars:
            fields[varname] = _ncvar_subset_to_dict(
                ncvars[name_map_fields[varname]], Good)
            _var_found(name_map_fields[varname])
        else:
            fields[varname] = None
            _var_not_found(name_map_fields[varname])

    # Save to output dictionary
    data['fields'] = fields

    ncFile.close()
    return data

###########################
# Create Variable methods #
###########################


def _get_wcl_name_map():
    '''Map out WCL variable names to AWOT dictionary.'''
    name_map = {
               'Time': 'Time',
               'time': 'time',
               'latitude': 'LAT',
               'longitude': 'LON',
               'altimeter_height': 'Ralt',
               'height': 'ALT',
               'roll': 'Roll',
               'pitch': 'Pitch',
               'temperature': 'trf',
               'pressure': 'pmb',
               'qc_flag': 'Prof_qc_flag',
               'range': 'Range',
               }
    return name_map

def _get_wcl_field_name_map():
    '''Map out WCL variable names to AWOT dictionary.'''
    name_map = {
               'copol_power': 'CopolPower',
               'cross_power': 'CrossPower',
               'copol_overlap': 'CopolOverlap',
               'cross_overlap': 'CrossOverlap',
               'copol_background': 'CopolBG',
               'cross_background': 'CrossBG',
               'copol_background_std': 'CopolBGSTD',
               'cross_background_std': 'CrossBGSTD',
               'copol_saturation': 'CopolSatur',
               'cross_saturation': 'CrossSatur',
               'zenith_angle': 'Zenith',
               'attitude_vectors': 'BeamVector',
               }
    return name_map

def _get_time(fname, ncFile, Good_Indices):
    """Pull the time from WCL file and convert to AWOT useable."""
    # Pull out the date, convert the date to a datetime friendly string
    # Now convert the time array into a datetime instance
    Time = num2date(ncFile.variables['time'][
                    Good_Indices], 'seconds since 1970-01-01 00:00:00+0:00')
    return Time

# def _nc_var_masked(ncFile, ncvar, Good_Indices):
#     """Convert a NetCDF variable into a masked variable."""
#     d = ncFile.variables[ncvar][Good_Indices]
#     np.ma.masked_invalid(d)
#     return d

# def _nc_var_to_dict(ncvar, Good_Indices):
#     """ Convert a NetCDF Dataset variable to a dictionary.
#     Appropriated from PyArt package.
#     """
#     d = dict((k, getattr(ncvar, k)) for k in ncvar.ncattrs())
#     d['data'] = ncvar[:]
#     if np.isscalar(d['data']):
#         # netCDF4 1.1.0+ returns a scalar for 0-dim array, we always want
#         # 1-dim+ arrays with a valid shape.
#         d['data'] = np.array(d['data'][Good_Indices, :])
#         d['data'].shape = (1, )
#     return d
