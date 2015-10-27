# -*- coding: utf-8 -*-
"""
awot.io.SondeDataFile
=====================

Class for reading Data from file.

"""

from __future__ import absolute_import, print_function
import os
from ..io.read_sounding_data import get_sounding, get_dropsonde

########################
#  BEGIN MAIN CODE
########################


def read_sonde(filename=None, instrument='dropsonde',
               initialize=False):
    '''
    Takes a filename pointing to a sounding data file
    and returns a SondeData object.

    Parameters::
    ------------
    filename : str
        Long path filename of data file.
    instrument : str
        Type (name) of instrument to process
        Currently the following arguments are valid:
        'dropsonde' - Dropsonde from aircraft
        'sounding' - Upsonde (rawinsonde) from ground

    To Do:
    Add support for:

    '''
    if filename is None:
        print("Need to provide a filename!")
        return

    reader = FileReader(filename=filename, instrument=instrument)

    if initialize:
        SondeData = AirborneData(reader)

        return SondeData
    else:
        return reader
# ##################################################


class FileReader(object):

    '''
    FileReader class to process data files.
    '''

    def __init__(self, filename=None, instrument=None):

        """
        Read in a file

        """

        if isinstance(filename, str) != False:
            if instrument.lower() == 'dropsonde':
                sonde = get_dropsonde(filename)
            elif instrument.lower() == 'sounding':
                sonde = get_sounding(filename)

            # Record the data into the variable
            self.sonde_data = sonde

        else:
            # Initializes class instance but leaves it to other methods to
            # populate the class attributes.
            return
