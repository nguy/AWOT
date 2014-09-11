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
def read_radar(fname):
    """Read in data using the pyart reader.
    PARAMETERS::
    ----------
        fname : string
            Filename [string]
    """
    rad = pio.read(fname)
    
        # Create a dictionary to transfer the data
    data = {'latitude': rad.latitude,
            'longitude': rad.longitude,
            'height': rad.altitude,
            'fields': rad.fields
            }
            
    return data