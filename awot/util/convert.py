"""
awot.util.convert
=================

These scripts are various convenience utilities.

"""
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
    awotradar ={'metadata': radar.metadata,
                 'longitude': radar.longitude,
                 'latitude': radar.latitude,
                 'height': radar.altitude,
                 'fields': radar.fields,
                 'platform': radar.metadata['instrument_name'],
                 'instrument': radar.metadata['instrument_name'],
                 'time': radar.time,
                 'data_format': 'ground'
                 }
    if instrument is not None:
        awotradar['instrument'] = instrument
    if platform is not None:
        awotradar['platform'] = platform
    return awotradar


def to_awot_flight(lon_dict=None, lat_dict=None, alt_dict=None,
                   time_dict=None, other_dict=None,
                   lon_array=None, lon_unit=None, lon_name=None, lon_title=None,
                   lat_array=None, lat_unit=None, lat_name=None, lat_title=None,
                   alt_array=None, alt_unit=None, alt_name=None, alt_title=None,
                   time_array=None,time_unit=None, time_name=None, time_title=None):
    """
    A convenience function to create an AWOT object from data loaded in some
    other fashion. For example if a reader does not exist in AWOT.

    Dictionaries may be provided or individual data can be passed and an
    attempt to build a minimal object will be attempted.

    Parameters
    ----------
    lon_dict : dict
        Dictionary to use for longitude.
        Should include 'data', 'units', and 'title' keys minimally.
    lat_dict : dict
        Dictionary to use for latitude.
        Should include 'data', 'units', and 'title' keys minimally.
    alt_dict : dict
        Dictionary to use for altitude.
        Should include 'data', 'units', and 'title' keys minimally.
    time_dict : dict
        Dictionary to use for time.
        Should include 'data', 'units', and 'title' keys minimally.
        Data should be datetime instances.
    other_dict : list
        List of other dictionaries to be included.
        If fields are passed, combine the dictionaries into a
        'fields' dictionary.

    lon_array : flt array
        Array of longitude values. AWOT uses masked Numpy arrays.
    lon_unit : str
        String indicating the units of longitude array.
    lon_name : str
        The string to assign to the 'name' key of longitude dictionary.
    lon_title : str
        The string to assign to the 'title' key of longitude dictionary.
    lat_array : flt array
        Array of latitude values. AWOT uses masked Numpy arrays.
    lat_unit : str
        String indicating the units of latitude array.
    lat_name : str
        The string to assign to the 'name' key of latitude dictionary.
    lat_title : str
        The string to assign to the 'title' key of latitude dictionary.
    alt_array : flt array
        Array of longitudinal values. AWOT uses masked Numpy arrays.
    alt_unit : str
        String indicating the units of altitude array.
    alt_name : str
        The string to assign to the 'name' key of altitude dictionary.
    alt_title : str
        The string to assign to the 'title' key of altitude dictionary.
    time_array : flt array
        Array of time values. AWOT uses datetime arrays.
    time_unit : str
        String indicating the units of time array.
    time_name : str
        The string to assign to the 'name' key of time dictionary.
    time_title : str
        The string to assign to the 'title' key of time dictionary.
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
        awot_flight = _build_flight_dict(awot_flight, 'longitude',
                        lon_array, lon_unit, lon_name, lon_title)
    if lat_array is not None:
        awot_flight = _build_flight_dict(awot_flight, 'latitude',
                        lat_array, lat_unit, lat_name, lat_title)
    if alt_array is not None:
        awot_flight = _build_flight_dict(awot_flight, 'altitude',
                        alt_array, alt_unit, alt_name, alt_title)
    if time_array is not None:
        awot_flight = _build_flight_dict(awot_flight, 'time',
                        time_array, time_unit, time_name, time_title)

    return awot_flight

def _build_flight_dict(flight, key, data, units, name, title):
    flight[key] = {'data': data,
                   'units' : units,
                   'full_name' : name,
                   'title' : title}
    return flight