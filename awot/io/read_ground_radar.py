"""
awot.io.read_ground_radar
=========================

These scripts are a wrapper using the pyart interface to access a large
number of file formats.
.

"""
# Load the needed packages
import pyart


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

        # Create a dictionary to transfer the data
        radar = {'metadata': rad.metadata,
                 'longitude': rad.longitude,
                 'latitude': rad.latitude,
                 'height': rad.altitude,
                 'fields': fields,
                 'platform': rad.metadata['instrument_name'],
                 'instrument': rad.metadata['instrument_name'],
                 'data_format': 'ground'
                 }
        if instrument is not None:
            radar['instrument'] = instrument
        if platform is not None:
            radar['platform'] = platform
    else:
        radar = rad
    return radar
