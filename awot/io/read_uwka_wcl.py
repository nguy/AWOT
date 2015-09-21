"""
awot.io.read_uwka_lidar
=========================

Scripts for processing data from the Wyoming Cloud Lidar.
 http://flights.uwyo.edu/wcl/

Tested 17 Aug 2015, may not be fully functional.
"""
# NOTES:: Testing was done on test data from 2014 flights.

# Load the needed packages
from __future__ import print_function
from netCDF4 import Dataset, num2date, date2num
import datetime
import numpy as np
import pytz


def read_wcl(fname):
    '''
    Read in NetCDF data file containing Wyoming cloud radar data.

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
            Height of center of radar range gate [km]
            altitude : float
            Aircraft altitude via GPS [km]
        tas : float
            Platform true airspeed [m/s]
        ground_speed : float
            Platform ground speed [m/s]
        dBZ_minimum : float
            Minimum detectable reflectivity at 1 km
        aspect : float
            WCR range gate / (WCR time integer * TAS)
        beam_vector : float
            (East, North, Up)  beam vector unit
        aircraft_wind : float
            In situ wind component at platform altitude along WCR beam
            Positive is away from radar

            fields : Dictionary of variables in file
            dBZ : float
                    Radar Equivalent Reflectivity Factor [dBZ]
                velocity : float
                Mean Doppler radial velocity [m/s]
                mask : int
                Target mask, see variable for notes
        metadata : Dictionary of global attributes in file
        project : str
            Project Name
        platform : str
            Platform name or identifier
        flight_number : str
            Flight number or identifer
    '''
    # Create a dictionary to hold data
    data = {}

    # Read the NetCDF
    ncFile = Dataset(fname, 'r')
    ncvars = ncFile.variables

    # Grab the metadata stored in global attributes as a dictionary
    metadata = ncFile.__dict__

    # Find the indices of not missing points
    Good = np.where(~np.isnan(ncFile.variables['time'][:]))

    Time = _get_time(fname, ncFile, Good)

    # Put time into output dictionary
    data['time'] = Time

    # Grab a name map for RASTA dynamic file data
    name_map_data = _get_wcl_name_map()

    # Loop through the variables and pull data
    for varname in name_map_data:
        if name_map_data[varname] in ncvars:
            data[varname] = _nc_var_masked(
                ncFile, name_map_data[varname], Good)
        else:
            data[varname] = None
            print(name_map_data[varname] + " does not exist in file...")

    if 'altrange' in ncvars:
        data['height'] = _nc_radar_var_to_dict(ncvars['altrange'], Good)
    else:
        print("No height variable in the file")

    # Add fields to their own dictionary
    fields = {}

    # Grab a name map for WCR field data
    name_map_fields = _get_wcr_field_name_map()

    # Loop through the variables and pull data
    for varname in name_map_fields:
        if name_map_fields[varname] in ncvars:
            fields[varname] = _nc_radar_var_to_dict(
                ncvars[name_map_fields[varname]], Good)
            print("Found " + name_map_fields[varname])
        else:
            fields[varname] = None
            print(name_map_fields[varname] + " does not exist in file...")

    # Save to output dictionary
    data['fields'] = fields

    # Pull out global attributes
    try:
        ncFile.ProjectName
        project = ncFile.ProjectName
    except:
        project = fname.split("_")[0]
    try:
        ncFile.FlightNumber
        flightnum = ncFile.DataDate
    except:
        flightnum = fname.split("_")[2]

    # Set the platform
    platform = ncFile.Platform

    # Create a dictionary to transfer the data
    data['metadata'] = metadata
    data['project'] = project
    data['platform'] = platform
    data['flight_number'] = flightnum
    data['data_format'] = 'wcr_vertical'

    ncFile.close()
    return data

###########################
# Create Variable methods #
###########################


def _get_wcl_name_map():
    '''Map out names used in RASTA microphysics NetCDF files to AWOT'''
    name_map = {
        'latitude': 'LAT',
        'longitude': 'LON',
        'tas': 'TAS',
               'ground_speed': 'GS',
               'altitude': 'ALT',
               'range': 'Range',
               'ralt': 'Ralt',
               'roll': 'Roll',
               'pitch': 'Pitch',
               'temperature': 'trf',
               'pressure': 'pmb',
               'bad_flag': 'Prof_qc_flag',

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
               'zenith': 'Zenith',
               'attitude': 'BeamVector'
    }

    # get the id's of the variables to be read
    Time_id = ncdf_varid(cdfid, 'Time')
    timeSec_id = ncdf_varid(cdfid, 'time')

    return name_map


def _get_wcr_field_name_map():
    '''
    Map out names used in RASTA dynamic NetCDF files to AWOT.
        'reflectivity' : Radar reflectivity
        'velocity'     : Mean radial Doppler velocity
        'wcrmask'      : Target mask
    '''
    name_map = {
        'reflectivity': 'reflectivity',
        'velocity': 'velocity',
        'mask': 'wcrmask',
    }
    return name_map


def _get_time(fname, ncFile, Good_Indices):
    """Pull the time from RASTA file and convert to AWOT useable."""
    # Pull out the date, convert the date to a datetime friendly string
    # Now convert the time array into a datetime instance
    Time = num2date(ncFile.variables['time'][
                    Good_Indices], 'seconds since 1970-01-01 00:00:00+0:00')
    return Time


def _nc_var_masked(ncFile, ncvar, Good_Indices):
    """Convert a NetCDF variable into a masked variable."""
    d = ncFile.variables[ncvar][Good_Indices]
    np.ma.masked_invalid(d)
    return d


def _nc_radar_var_to_dict(ncvar, Good_Indices):
    """ Convert a NetCDF Dataset variable to a dictionary.
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
