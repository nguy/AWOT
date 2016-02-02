"""
awot.io.read_ground_radar
=========================

These scripts are a wrapper using the pyart interface to access a large
number of file formats.

"""
# Load the needed packages
import pyart
from netCDF4 import num2date, date2num
from . import common


def read_ground_radar(fname, map_to_awot=True,
                      instrument=None, platform=None):
    """
    A wrapper using the Py-ART read interface.

    In this structure the 'fields' attribute will be
    structured as in the PyArt package.

    Parameters
    ----------
    fname : string
        Filename.
    map_to_awot : bool
        If True it maps the input to the AWOT radar structure [default].
        If False the Py-ART radar instance is maintained.
    instrument : str
        If set this supersedes the instrument key in AWOT dictionary.
    platform : str
        If set this supersedes the platform key in AWOT dictionary.
    """
    rad = pyart.io.read(fname)

    if map_to_awot:
        # build the fields dictionary
        fields = {}
        for fldName in rad.fields:
            fields[fldName] = rad.fields[fldName]

        Time = _convert_pyart_time(rad)
        # Create a dictionary to transfer the data
        radar = {'metadata': rad.metadata,
                 'longitude': rad.gate_longitude,
                 'latitude': rad.gate_latitude,
                 'height': rad.gate_z,
                 'fields': fields,
                 'platform': rad.metadata['instrument_name'],
                 'instrument': rad.metadata['instrument_name'],
                 'time': Time,
                 'data_format': 'ground'
                 }
        if instrument is not None:
            radar['instrument'] = instrument
        if platform is not None:
            radar['platform'] = platform
    else:
        radar = rad
    return radar


def _convert_pyart_time(radar):
    """Pull the time from WCR NetCDF file and convert to AWOT useable."""

    # Pull out the date, convert the date to a datetime friendly string

    # Now convert the time array into a datetime instance
    timedate = pyart.util.datetime_utils.datetimes_from_radar(radar)
    Timesec = date2num(timedate, common.EPOCH_UNITS)
    Time_unaware = num2date(Timesec, common.EPOCH_UNITS)
    Time = {'data': Time_unaware, 'units': common.EPOCH_UNITS,
            'title': 'Time', 'full_name': 'Time (UTC)'}
    return Time
