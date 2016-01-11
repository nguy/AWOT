"""
awot.util.convert
=================

These scripts are various convenience utilities.

"""
from netCDF4 import num2date

from ..io.common import _build_dict

def pyart_radar_to_awot(radar, instrument=None, platform=None):
    """
    This function creates an AWOT radar instance from a Py-ART
    radar instance.

    Parameters
    ----------
    radar : Py-ART radar instance
        A Py-ART radar instance created using the Py_ART read
        functions.
    instrument : str
        If set this supersedes the instrument key in AWOT dictionary.
    platform : str
        If set this supersedes the platform key in AWOT dictionary.
    """
    awotradar ={'metadata': radar.metadata.copy(),
                 'longitude': radar.longitude.copy(),
                 'latitude': radar.latitude.copy(),
                 'height': radar.altitude.copy(),
                 'fields': radar.fields.copy(),
                 'platform': radar.metadata['instrument_name'],
                 'instrument': radar.metadata['instrument_name'],
                 'time': radar.time.copy(),
                 'data_format': 'ground'
                 }
    Time = num2date(radar.time['data'], radar.time['units'])
    awotradar['time']['data'] = Time
    if instrument is not None:
        awotradar['instrument'] = instrument
    if platform is not None:
        awotradar['platform'] = platform
    return awotradar


def to_awot_flight(lon_dict=None, lat_dict=None, alt_dict=None,
                   time_dict=None, other_dict=None,
                   lon_array=None, lon_unit=None, lon_longname=None, lon_stdname=None,
                   lat_array=None, lat_unit=None, lat_longname=None, lat_stdname=None,
                   alt_array=None, alt_unit=None, alt_longname=None, alt_stdname=None,
                   time_array=None,time_unit=None, time_longname=None, time_stdname=None):
    """
    A convenience function to create an AWOT object from data loaded in some
    other fashion. For example if a reader does not exist in AWOT.

    Dictionaries may be provided or individual data can be passed and an
    attempt to build a minimal object will be attempted.

    Parameters
    ----------
    lon_dict : dict
        Dictionary to use for longitude.
        Should include 'data', 'units', and 'standard_name' keys minimally.
    lat_dict : dict
        Dictionary to use for latitude.
        Should include 'data', 'units', and 'standard_name' keys minimally.
    alt_dict : dict
        Dictionary to use for altitude.
        Should include 'data', 'units', and 'standard_name' keys minimally.
    time_dict : dict
        Dictionary to use for time.
        Should include 'data', 'units', and 'standard_name' keys minimally.
        Data should be datetime instances.
    other_dict : list
        List of other dictionaries to be included.
        If fields are passed, combine the dictionaries into a
        'fields' dictionary.

    lon_array : flt array
        Array of longitude values. AWOT uses masked Numpy arrays.
    lon_unit : str
        String indicating the units of longitude array.
    lon_longname : str
        The string to assign to the 'name' key of longitude dictionary.
    lon_stdname : str
        The string to assign to the 'standard_name' key of longitude dictionary.
    lat_array : flt array
        Array of latitude values. AWOT uses masked Numpy arrays.
    lat_unit : str
        String indicating the units of latitude array.
    lat_longname : str
        The string to assign to the 'name' key of latitude dictionary.
    lat_stdname : str
        The string to assign to the 'standard_name' key of latitude dictionary.
    alt_array : flt array
        Array of longitudinal values. AWOT uses masked Numpy arrays.
    alt_unit : str
        String indicating the units of altitude array.
    alt_longname : str
        The string to assign to the 'name' key of altitude dictionary.
    alt_stdname : str
        The string to assign to the 'standard_name' key of altitude dictionary.
    time_array : flt array
        Array of time values. AWOT uses datetime arrays.
    time_unit : str
        String indicating the units of time array.
    time_longname : str
        The string to assign to the 'name' key of time dictionary.
    time_stdname : str
        The string to assign to the 'standard_name' key of time dictionary.
    """
    awot_flight = {}

    # Build AWOT object if dictionaries are provided
    if lon_dict is not None:
        awot_flight['longitude'] = lon_dict
    if lon_dict is not None:
        awot_flight['latitude'] = lat_dict
    if lon_dict is not None:
        awot_flight['altitude'] = alt_dict
    if lon_dict is not None:
        awot_flight['time'] = time_dict
    if lon_dict is not None:
        awot_flight['longitude'] = lon_dict

    for d in other_dict:
        awot_flight[d] = other_dict[d]

    # Attempt to build AWOT dictionaries if arrays are provided
    if lon_array is not None:
        awot_flight = _build_dict(awot_flight, 'longitude',
                        lon_array, lon_unit, lon_longname, lon_stdname)
    if lat_array is not None:
        awot_flight = _build_dict(awot_flight, 'latitude',
                        lat_array, lat_unit, lat_longname, lat_stdname)
    if alt_array is not None:
        awot_flight = _build_dict(awot_flight, 'altitude',
                        alt_array, alt_unit, alt_longname, alt_stdname)
    if time_array is not None:
        awot_flight = _build_dict(awot_flight, 'time',
                        time_array, time_unit, time_longname, time_stdname)

    return awot_flight

def _build_dict(flight, key, data, units, longname, stdname):
    flight[key] = {'data': data,
                   'units' : units,
                   'long_name' : longname,
                   'standard_name' : stdname}
    return flight