"""
    awot.io.read_hiaper_radar
    =================
    
    Scripts to read NCAR HIAPER Cloud Radar NetCDF data files.
    
    https://www.eol.ucar.edu/observing_facilities/hiaper-gulfstream-gv
    
    Testing was done from data file provided by Pei Tsai.
    """

from netCDF4 import Dataset, num2date, date2num
import datetime
import numpy as np
from . import common


def read_hcr(fname, field_mapping=None, file_mapping=None):
    '''
    Read NCAR HIAPER Cloud Radar HCR NetCDF data file.
    
    ATTENTION!
    There are extra variables that are ignored by the reader at this time.
    
    Parameters
    ----------
    fname : str
        Filename
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
    height : float
        Height_agl of center of radar range gate [km].
    roll: float
        Aircraft roll attitude [decimal degrees].
    pitch: float
        Aircraft pitch attitude[decimal degrees].
    drift: float
        Aircraft drift.
    roatation: float
        Radar rotation rae for CF radial components [decimal degrees].
    tilt: float 
        Radar tilt angle [decimal degrees].
    height: float
        Beam range height[km].
    project : str
        Project name.
    platform : str
        Platform name or identifier.
    intrument_name: str
        Name of the instrument HCR
    '''
    # Create a dictionary to hold data
    data = {}
    
    # Read the NetCDF
    ncFile = Dataset(fname, 'r')
    ncvars = ncFile.variables
    
    # Grab the metadata stored in global attributes as a dictionary
    try:
        data['metadata'] = ncFile.__dict__
    except:
        data['metadata'] = None
    
    try:
        data['project'] = ncFile.institution
    except:
        data['project'] = None

    try:
        data['flight_number'] = ncFile.start_datetime
    except:
        data['flight_number'] = None
    
    # Set the platform
    try:
        data['platform'] = ncFile.site_name
    except:
        data['platform'] = 'GV'

    try:
        data['intrument_name'] = ncFile.intrument_name
    except:
        data['platform'] = 'HCR'


    # Find the indices of not missing points
    Good = np.where(~np.isnan(ncFile.variables['time'][:]))

    data['time'] = _get_time(fname, ncFile, Good)
    
    # Grab the name map
    name_map_data = _get_hcr_data_map()
    
    # Loop through the variables and pull data
    for varname in name_map_data:
        if name_map_data[varname] in ncvars:
            data[varname] = common._ncvar_to_dict(
                ncvars[name_map_data[varname]])
        else:
            data[varname] = None
            common._var_not_found(varname)
    # Add fields to their own dictionary
    fields = {}

    # Grab a name map for HCR field data
    name_map_fields = _get_hcr_field_name_map()
    
    # Loop through the variables and pull data
    for varname in name_map_fields:
        if name_map_fields[varname] in ncvars:
            fields[varname] = common._ncvar_subset_to_dict(
                ncvars[name_map_fields[varname]], Good)
        else:
            fields[varname] = None
            common._var_not_found(varname)

    #cant remember do we need this SURFACE DATA?

    # Find the surface variable
    # See PROCESSED_DATA/WCR_L2_OWLES13.20131015.cdl
    # for details of mask properties
    ##surface = np.empty_like(data['latitude']['data'])
    #    condition = np.equal(fields['mask']['data'], 32)
    #    for nn in range(len(surface)):
    #        if np.any(condition[nn, :]):
    #            surface[nn] = data['height']['data'][np.where(
    #                                                          condition[nn, :])[0][0]]
    #    data['surface'] = {'name': "surface",
    #        'long_name': "Height of Surface",
    #            'data': surface,
    #                'units': 'meters'
    #                }

    # Save to output dictionary
    data['fields'] = fields
        
        # Set the data format
    try:
        data['data_format'] = ncFile.version
    except:
        data['data_format'] = 'hcr_vertical'
        
        ncFile.close()
        
    return data

##################
#  _get methods  #
##################


def _get_hcr_data_map():
    '''Map HCR variable names to AWOT dictionary.'''
    name_map = {
    'latitude': 'latitude',
    'longitude': 'longitude',
    'altitude': 'altitude_agl',
    'heading':'heading',
    'roll':'roll',
    'pitch':'pitch',
    'drift':'drift',
    'rotation':'rotation',
    'tilt': 'tilt',
    'height': 'range',
    }
    return name_map


def _get_hcr_field_name_map():
    '''Map HCR radar variables to AWOT dictionary.'''
    name_map = {
    'reflectivity': 'DBZ',
    'reflectivity_hc': 'DBZHC',
    'reflectivity_vc': 'DBZVC',
    'velocity': 'VEL',
    'raw_velocity': 'VEL_RAW',
    'spectrum_width': 'WIDTH',
    'snr': 'SNR',
    'ncp': 'NCP',
    'ldr': 'LDR',
    'ldrh': 'LDRH',
    'ldrv': 'LDRV',
    }
    return name_map


def _get_time(fname, ncFile, Good_Indices):
    """Pull the time from HIAPER NetCDF file and convert to AWOT useable."""
    
    # Pull out the date, convert the date to a datetime friendly string
    
    # Now convert the time array into a datetime instance
    
    
    #bad indent
    Time_unaware = num2date(ncFile.variables['time'][Good_Indices],
                            ncFile.variables['time'].units)
                            
    Time = {'data': Time_unaware, 'units': common.EPOCH_UNITS,
            'title': 'Time', 'full_name': 'Time (UTC)'}
    return Time
