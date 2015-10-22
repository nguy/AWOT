"""
common.py
"""

from __future__ import print_function
import numpy as np

######################
#  variable methods  #
######################


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

def _get_epoch_units():
    """Set common time units for AWOT. Using Epoch."""
    return 'seconds since 1970-1-1 00:00:00+0:00'

def _var_found(var):
    '''Print variable found message.'''
    print("Found %s" % var)

def _var_not_found(var):
    '''Print variable not found message.'''
    print("%s does not exist in file..." % var)

##############################
#  image save methods  #
##############################

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
