"""
awot.util.helper
=================

These scripts are various convenience utilities.

"""
from netCDF4 import num2date

from ..graph.common import _get_start_datetime, _get_end_datetime
from ..io.common import _build_dict

def time_subset_awot_dict(time, data, start_time, end_time):
    '''
    Get the variable from the fields dictionary.
    Subset the time when in time series format.

    Parameters
     ----------
    time : dict
        AWOT time dictionary
    data : dict
        AWOT data dictionary.
    start_time : str
        UTC time to use as start time for subsetting in datetime format.
        (e.g. 2014-08-20 12:30:00)
    end_time : str
        UTC time to use as an end time for subsetting in datetime format.
        (e.g. 2014-08-20 16:30:00)
    '''
    # Check to see if time is subsetted
    dt_start = _get_start_datetime(time, start_time)
    dt_end = _get_end_datetime(time, end_time)
    datasub = data.copy()
    datasub['data'] = data['data'][(time['data'] >= dt_start) &
                   (time['data'] <= dt_end)]
    return datasub

def add_dict_to_awot(awot, keyname, newdict=None, data=None, units=None,
                     longname=None, stdname=None):
    '''
    Add a dictionary to an AWOT data instance.

    Parameters
    ----------
    awot : dict
        AWOT data dictionary instance.
    keyname : str
        The key to be used when adding to the AWOT dictionary.
    newdict : dict
        Optional AWOT dictionary if built by user. Can be passed and
        added to main AWOT instance. If None, a dictionary will be built
        from the variables passed below.
    data : array
        The data associated with keyname dictionary. Optional.
    units : str
        The units associated with keyname dictionary. Optional.
    longname : str
        The long name associated with keyname dictionary. Optional.
    stdname : str
        The standard name associated with keyname dictionary. Optional.
    '''
    if newdict is None:
       newdict = _build_dict(data, units, longname, stdname)
    awot[keyname] = newdict
    return

def add_dict_to_awot_fields(awot, keyname, newdict=None, data=None, units=None,
                     longname=None, stdname=None):
    '''
    Add a dictionary to the fields dictionary in an AWOT data instance.

    Parameters
    ----------
    awot : dict
        AWOT data dictionary instance.
    keyname : str
        The key to be used when adding to the AWOT dictionary.
    newdict : dict
        Optional AWOT dictionary if built by user. Can be passed and
        added to main AWOT instance. If None, a dictionary will be built
        from the variables passed below.
    data : array
        The data associated with keyname dictionary. Optional.
    units : str
        The units associated with keyname dictionary. Optional.
    longname : str
        The long name associated with keyname dictionary. Optional.
    stdname : str
        The standard name associated with keyname dictionary. Optional.
    '''
    if newdict is None:
       newdict = _build_dict(data, units, longname, stdname)
    awot['fields'][keyname] = newdict
    return