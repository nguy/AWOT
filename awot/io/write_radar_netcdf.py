"""
awot.io.write_radar_netcdf
=========================

A group of scripts to write radar data to NetCDF.
Especially that collected by the NOAA P-3 aircraft.
Supports both tail Doppler and lower fuselage radars. 

Created by Nick Guy.

"""
# NOTES:: This has only been tested with DYNAMO data files, versions
#         may change and another function may be needed.
# HISTORY::
#   8 Jan 2014 - Nick Guy.   NRC, NOAA/NSSL (nick.guy@noaa.gov)   
#                Converted NCL functions below to Python
#-------------------------------------------------------------------
# Load the needed packages
import netCDF4 as nc4
import numpy as np
#-------------------------------------------------------------------
# Begin methods
######################
# TDR file methods #
######################
def radar2nc(radar, Outfile=None):
    """Write a NetCDF data file with data in a radar dictionary.
        
    Parameters::
    ----------
    radar : dict
        Dictionary of data retrieved from an input reader 
        (e.g. read_p3_radar.py)
    Outfile : str
        String name for output netcdf file
    """
    if Outfile is None:
        Outfile = radar['metadata']['Flight_ID'].replace(" ", "") + '_windsyn'
    nc_fid = nc4.Dataset(Outfile + '.nc', 'w', format='NETCDF4')
    nc_fid.description = "Airborne radar data NetCDF"
    
    # Define dimensions
    Imax = len(radar['longitude']['data'][:])
    Jmax = len(radar['latitude']['data'][:])
    Kmax = len(radar['height']['data'][:])
    
    pid = nc_fid.createDimension('property', 1)
    xid = nc_fid.createDimension('lon', Imax)
    yid = nc_fid.createDimension('lat', Jmax)
    zid = nc_fid.createDimension('height', Kmax)
    
    # Set global attributes
    nc_fid.source = radar['metadata']['source_file']
    nc_fid.creation_date = radar['metadata']['creation_date'].isoformat()
    nc_fid.number_radars = radar['metadata']['number_radars']
    nc_fid.Nmosm = radar['metadata']['Nmosm']
    nc_fid.dBZ_calculation = radar['metadata']['dBZ_calculation']
    nc_fid.Integration_direction = radar['metadata']['Integration_direction']
    nc_fid.dual_Doppler_Analysis = radar['metadata']['dual_Doppler_Analysis']
    nc_fid.number_smooth = radar['metadata']['number_smooth']
    nc_fid.Obrien_div_correction = radar['metadata']['Obrien_div_correction']
    nc_fid.Flight_ID = radar['metadata']['Flight_ID']
    nc_fid.title = "Windsyn computatons of P-3 tail radar data"
    nc_fid.Start_datetime = radar['datetime_start'].isoformat() + 'Z'
    nc_fid.End_datetime = radar['datetime_end'].isoformat() + 'Z'
    nc_fid.platform = radar['platform']
    nc_fid.instrument = radar['instrument']
    
    # Create Output variables
    lonid = nc_fid.createVariable('Lon', np.float32, ('lon',))
    lonid.units = radar['longitude']['units']
    lonid.long_name = radar['longitude']['long_name']
    
    latid = nc_fid.createVariable('Lat', np.float32, ('lat',))
    latid.units = radar['latitude']['units']
    latid.long_name = radar['latitude']['long_name']
    
    htid = nc_fid.createVariable('Height', np.float32, ('height',))
    htid.units = radar['height']['units']
    htid.long_name = radar['height']['long_name']
    
    # Loop through the fields to create variables for each
    for variable in radar['fields'].keys():
        if radar['fields'][variable] is not None:
            vid = nc_fid.createVariable(variable, np.float32, ('height', 'lat', 'lon'))
            vid.units = radar['fields'][variable]['units']
            vid.long_name = radar['fields'][variable]['long_name']
            vid.fill_value = radar['fields'][variable]['_FillValue']
            vid[:] = radar['fields'][variable]['data'][:]
        
    # Close the NetCDF file
    nc_fid.close()