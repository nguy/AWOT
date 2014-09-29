"""
awot.io.read_ground_radar
=========================

These scripts are a wrapper using the pyart interface to access a large
number of file formats.  
In this structure the 'fields' attribute will be structured as in 
the PyArt package.  

Created by Nick Guy.
04 Sep 2014.

"""
#-------------------------------------------------------------------
# Load the needed packages
import pyart.io as pio
#===============================================================
# BEGIN FUNCTIONS
#===============================================================
def read_radar(fname, instrument=None, platform=''):
    """Read in data using the pyart reader.
    PARAMETERS::
    ----------
        fname : string
            Filename [string]
    """
    if instrument is None:
        instrument = 'ground'
    
    rad = pio.read(fname)
    
    # build the fields dictionary
    fields = {}
    for fldName in rad.fields:
        fields[fldName] = rad.fields[fldName]

    # Create a dictionary to transfer the data
    data = {'metadata': rad.metadata,
            'longitude': rad.longitude,
            'latitude': rad.latitude,
            'height': rad.altitude,
            'fields': fields,
            'platform': rad.metadata['instrument_name'],
            'instrument': instrument
            }
            
    return data