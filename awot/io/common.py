"""
awot.io.common
==============

Common IO routines.

"""

from __future__ import print_function
import numpy as np
from netCDF4 import num2date, date2num

#################################
#  variable/dictionary methods  #
#################################
EPOCH_UNITS = 'seconds since 1970-1-1T00:00:00+0:00'
EPOCH_UNITS2 = 'seconds since 1970-1-1 00:00:00+0:00'


def _build_dict(data, units, longname, stdname):
    newdict = {'data': data,
               'units': units,
               'long_name': longname,
               'standard_name': stdname}
    return newdict


def _ncvar_subset_masked(ncFile, ncvar, Good_Indices):
    """
    Convert a NetCDF variable into a masked variable.
    Assumes a 1D variable
    """
    d = ncFile.variables[ncvar][Good_Indices]
    np.ma.masked_invalid(d)
    return d


def _ncvar_subset_to_dict(ncvar, Good_Indices):
    """
    Convert a NetCDF Dataset variable to a dictionary.
    Appropriated from Py-ART package.
    Assumes subsetting in first column.
    """
    d = dict((k, getattr(ncvar, k)) for k in ncvar.ncattrs())
    d['data'] = ncvar[:]
    if np.isscalar(d['data']):
        # netCDF4 1.1.0+ returns a scalar for 0-dim array, we always want
        # 1-dim+ arrays with a valid shape.
        d['data'] = np.array(d['data'][Good_Indices, :])
        d['data'].shape = (1, )
    return d


def _ncvar_to_dict(ncvar):
    """
    Convert a NetCDF Dataset variable to a dictionary.
    Appropriated from PyArt package.
    """
    d = dict((k, getattr(ncvar, k)) for k in ncvar.ncattrs())
    d['data'] = ncvar[:]
    if np.isscalar(d['data']):
        # netCDF4 1.1.0+ returns a scalar for 0-dim array, we always want
        # 1-dim+ arrays with a valid shape.
        d['data'] = np.array(d['data'][:])
        d['data'].shape = (1, )
    return d


def _ncvar_to_dict_masked(ncvar, Good_Indices):
    """
    Convert a NetCDF Dataset variable to a dictionary.
    Appropriated from PyArt package.
    """
    d = dict((k, getattr(ncvar, k)) for k in ncvar.ncattrs())
    d['data'] = ncvar[Good_Indices]
    if np.isscalar(d['data']):
        # netCDF4 1.1.0+ returns a scalar for 0-dim array, we always want
        # 1-dim+ arrays with a valid shape.
        d['data'] = np.array(d['data'][:])
        d['data'].shape = (1, )
    return d


def _nasa_ames_var_to_dict(var, standard_name, long_name):
    d = {}
    d['standard_name'] = standard_name
    d['long_name'] = long_name
    d['units'] = " "
    d['data'] = var
    return d


def _h5var_to_dict(dataset, units=None, long_name=None, standard_name=None):
    """ Convert an HDF5 Dataset to a dictionary."""
    d = {}
    if dataset.dtype.char == "S":
        d['data'] = np.array(dataset)[0]
    else:
        d['data'] = np.array(dataset)
    if len(dataset.attrs) > 0:
        for attrname in list(dataset.attrs):
            d[attrname] = dataset.attrs.get(attrname)
    else:
        d['standard_name'] = standard_name
        d['long_name'] = long_name
        d['units'] = units
    return d


def _var_found(var):
    '''Print variable found message.'''
    print("Found %s" % var)


def _var_not_found(var):
    '''Print variable not found message.'''
    print("%s does not exist in file..." % var)

##################
#  time methods  #
##################


def _get_epoch_dict(timesec, time_units):
    '''Output AWOT standard Epoch time dictionary.'''
    # Convert the time array into a datetime instance
    dt = num2date(timesec, time_units)
    # Now convert this datetime instance into a number of seconds since Epoch
    eptime = date2num(dt, EPOCH_UNITS)
    # Now once again convert this data into a datetime instance
    dt_ep = num2date(eptime, EPOCH_UNITS)
    epdict = {'data': dt_ep, 'units': EPOCH_UNITS,
              'standard_name': 'Time', 'long_name': 'Time (UTC)'}
    return epdict


def convert_to_epoch_dict(datetime_dict):
    '''Output AWOT standard Epoch time dictionary.'''
    # Check to see if already using AWOT epoch units
    if ((datetime_dict['units'] == EPOCH_UNITS) or
        (datetime_dict['units'] == EPOCH_UNITS2)):
        dt_ep = datetime_dict['data']
    else:
        # Convert the datetime instance into a number array
        try:
            eptime = date2num(datetime_dict['data'], EPOCH_UNITS)
        except:
            dt = num2date(datetime_dict['data'], datetime_dict['units'])
            eptime = date2num(dt, EPOCH_UNITS)
        # Convert number array into AWOT epoch datetime instance
        dt_ep = num2date(eptime, EPOCH_UNITS)
    # Build dictionary to return
    epdict = {'data': dt_ep, 'units': EPOCH_UNITS,
              'standard_name': 'Time', 'long_name': 'Time in %s' % EPOCH_UNITS}
    return epdict


def check_epoch_units(datetime_dict):
    '''Check to see if dictionary has AWOT epoch units.'''
    if ((datetime_dict['units'] == EPOCH_UNITS) or
        (datetime_dict['units'] == EPOCH_UNITS2)):
        return True
    else:
        return False

########################
#  image save methods  #
########################


def save_figure(self, figName='awot_plot', figType='png', **kwargs):
    '''Save the current plot

    Parameters
    ----------
    figName : str
        Figure name
    figType : str
        Figure format, default to .png
    '''
    plt.gca()
    plt.gcf()
    plt.savefig(figName+'.'+figType, format=figType)
    print("Saved figure: %s.%s" % (figName, figType))

    # Now close the plot to make sure matplotlib is happy
    plt.close()
