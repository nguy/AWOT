"""
common.py
"""

from __future__ import print_function
import numpy as np

################################
#   variable convert methods  ##
################################


def _nc_var_masked(ncFile, ncvar, Good_Indices):
    """Convert a NetCDF variable into a masked variable."""
    d = ncFile.variables[ncvar][Good_Indices]
    np.ma.masked_invalid(d)
    return d


def _nc_radar_var_to_dict(ncvar, Good_Indices):
    """
    Convert a NetCDF Dataset variable to a dictionary.
    Appropriated from Py-ART package.
    """
    d = dict((k, getattr(ncvar, k)) for k in ncvar.ncattrs())
    d['data'] = ncvar[:]
    if np.isscalar(d['data']):
        # netCDF4 1.1.0+ returns a scalar for 0-dim array, we always want
        # 1-dim+ arrays with a valid shape.
        d['data'] = np.array(d['data'][Good_Indices, :])
        d['data'].shape = (1, )
    return d


def _get_time_units():
    """Set common time units for AWOT. Using Epoch."""
    return 'seconds since 1970-1-1 00:00:00+0:00'


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
    print("Saved figure: " + figName+'.'+figType)

    # Now close the plot to make sure matplotlib is happy
    plt.close()
