"""
awot.io.read_ground_radar
=========================

These scripts are a wrapper using the pyart interface to access a large
number of file formats.  
.  

"""
# Load the needed packages
import pyart

def read_radar(fname, instrument=None, platform=''):
    """
    A wrapper using the pyart read interface.
    
    In this structure the 'fields' attribute will be 
    structured as in the PyArt package
    
    Parameters
    ----------
        fname : string
            Filename [string]
    """
    rad = pyart.io.read(fname)

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
    return radar