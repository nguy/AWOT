"""
awot.io.read_radar_sweep
=========================

These scripts are a wrapper using the pyart interface to access a large
number of file formats.

"""
# Load the needed packages
try:
    import pyart
    _PYART_AVAILABLE = True
except ImportError:
    _PYART_AVAILABLE = False


def read_tdr_sweep(fname, map_to_awot=True,
                   instrument=None, platform=None):
    """
    A wrapper using the Py-ART read interface.

    In this structure the 'fields' attribute will be
    structured as in the Py-ART package.

    Parameters
    ----------
    fname : str
        Long path filename.
    map_to_awot : bool
        If True it maps the input to the AWOT radar structure [default].
        If False the Py-ART radar instance is maintained.
    instrument : str
        If set this supersedes the instrument key in AWOT dictionary.
    platform : str
        If set this supersedes the platform key in AWOT dictionary.
    """
    if not _PYART_AVAILABLE:
        raise MissingOptionalDependency(
            "pyart is required to use read_tdr_sweep but is not installed")
    rad = pyart.io.read(fname)

    if map_to_awot:
        # build the fields dictionary
        fields = {}
        for fldName in rad.fields.keys():
            fields[fldName] = rad.fields[fldName]

        # Create a dictionary to transfer the data
        radar = {'metadata': rad.metadata,
                 'longitude': rad.longitude,
                 'latitude': rad.latitude,
                 'altitude': rad.altitude,
                 'fields': fields,
                 'range': rad.range,
                 'rotation': rad.rotation,
                 'drift': rad.drift,
                 'heading': rad.heading,
                 'pitch': rad.pitch,
                 'roll': rad.roll,
                 'tilt': rad.tilt,
                 'platform': rad.metadata['instrument_name'],
                 'instrument': rad.metadata['instrument_name'],
                 'data_format': 'tdr_sweep'
                 }
        if instrument is not None:
            radar['instrument'] = instrument
        if platform is not None:
            radar['platform'] = platform
    else:
        radar = rad
    return radar


def read_lf_sweep(fname, map_to_awot=True,
                  instrument=None, platform=None):
    """
    A wrapper using the Py-ART read interface to read in
    lower fuselage radar (e.g. NOAA P-3) data files.
    Both cfradial and raw Sigmet formats are accepted.


    In this structure the 'fields' attribute will be
    structured as in the PyArt package.

    Parameters
    ----------
    fname : str
        Long path filename.
    map_to_awot : bool
        If True it maps the input to the AWOT radar structure [default].
        If False the Py-ART radar instance is maintained.
    instrument : str
        If set this supersedes the instrument key in AWOT dictionary.
    platform : str
        If set this supersedes the platform key in AWOT dictionary.
    """
    if not _PYART_AVAILABLE:
        raise MissingOptionalDependency(
            "pyart is required to use read_lf_sweep but is not installed")
    rad = pyart.io.read(fname)

    if map_to_awot:
        # build the fields dictionary
        fields = {}
        for fldName in rad.fields.keys():
            fields[fldName] = rad.fields[fldName]

        # Create a dictionary to transfer the data
        radar = {'metadata': rad.metadata,
                 'longitude': rad.longitude,
                 'latitude': rad.latitude,
                 'altitude': rad.altitude,
                 'fields': fields,
                 'range': rad.range,
                 'rotation': rad.rotation,
                 'drift': rad.drift,
                 'heading': rad.heading,
                 'pitch': rad.pitch,
                 'roll': rad.roll,
                 'tilt': rad.tilt,
                 'platform': rad.metadata['instrument_name'],
                 'instrument': rad.metadata['instrument_name'],
                 'data_format': 'tdr_sweep'
                 }
        if instrument is not None:
            radar['instrument'] = instrument
        if platform is not None:
            radar['platform'] = platform
    else:
        radar = rad
    return rad
