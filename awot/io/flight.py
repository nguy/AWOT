"""
flight.py - Routines for reading flight data from file
"""
from __future__ import absolute_import, print_function
import os
import numpy as np

from datetime import datetime
from netCDF4 import Dataset, num2date, date2num
from ..io.common import (_get_epoch_units,
                         _ncvar_to_dict, _var_not_found)
from ..io.name_maps_flight import _get_name_map

#########################
#   NetCDF Methods      #
#########################


def read_netcdf(fname, time_var=None, mapping_dict=None, platform=None):
    """
    Read in NetCDF formatted flight data.
    Output variable names are controlled by mapping_dict and platform
    keywords. If nothing is chosen a direct copy of names in file will
    be output.

    This reader checks to see if NCAR Research Aircraft Facility
    (RAF) Nimbus conventions are followed. This convention packs 
    data with a higher rate than 1 Hz into a 2-dimensional array.
    See the following website for a description of this data
    convention: http://www.eol.ucar.edu/raf/Software/netCDF.html

    Parameters
    ----------
    fname : str
        Filename.
    time_var : str
        Name of time variable. Choosing this overrides any mapping 
        dictionary that is used.
    mapping_dict : dict
        Dictionary to use for mapping variable names. Use if provided.
        If None, use default.
    platform: str
        Platform name. If no mapping dictionary is provided, then
        this can provide an internal default mapping dictionary.
    """
    # Read the NetCDF
    ncFile = Dataset(fname, 'r')

    # If this is a T-28 file, use a separate read function
##    t28_names = ['t28', 't-28', 'sdsmt', 'sdsmtt28', 'sdsm&t',
##                 'sdsm&tt-28', 'sdsmtt-28', 'sdsm&tt28']
##    if platform.lower().replace(" ", "") in t28_names:
##        from .read_t28 import read_t28_netcdf
##        return read_t28_netcdf(ncFile)

    # Check to see if this file follows RAF Nimbus conventions
    try:
        jnk = ncFile.Conventions
        isRAF, RAFrate, RAFdim = _determine_RAF(ncFile)
    except:
        isRAF, RAFrate, RAFdim = False, None, None

    # Grab a name map for data
    if mapping_dict is not None:
        name_map = mapping_dict
    else:
        # Check to see if platform is specified
        if platform is not None:
            name_map = _get_name_map(platform)
        else:
            name_map = _make_name_map_from_varlist(ncFile.variables.keys())

    # Cycle to through variables in file
    data = _make_data_dictionary(ncFile, name_map, isRAF, RAFrate=RAFrate)

    # Calculate U,V wind if not present
    if 'Uwind' not in name_map:
        Uwind, Vwind = _winduv(data)
        # Add to the dictionary
        data['Uwind'] = Uwind
        data['Vwind'] = Vwind

    # Time can be a fickle little beast, so even if it is in name
    # map, we need to massage it into a format we can work with

    if time_var is not None:
        data['time'] = _get_time(
            ncFile, isRAF, RAFrate=RAFrate, timevar=time_var)
    else:
        if 'time' in data.keys():
            data['time'] = _get_time(
                ncFile, isRAF, RAFrate=RAFrate, timevar=name_map['time'])
        else:
            data['time'] = _get_time(
                ncFile, isRAF, RAFrate=RAFrate, timevar=None)

    # Pull out global attributes
    try:
        ncFile.ProjectName
        data['project'] = ncFile.ProjectName
    except:
        data['project'] = None
    try:
        ncFile.Platform
        data['platform'] = ncFile.Platform
    except:
        data['platform'] = None
    try:
        ncFile.FlightNumber
        data['flight_number'] = ncFile.FlightNumber
    except:
        data['flight_number'] = None

    ncFile.close()
    return data


def read_netcdf_variable(fname, Rec):
    """
    Read a single variable from a NetCDF data file.
    Parameters
    ----------
     fname : string
         Filename.
     Rec : string
         Variable name to be pulled out [string].
    Output
    ------
     VarOut : float
         Masked array containing variable data.
      """
    # Read the NetCDF, grab the variable, close file
    ncFile = Dataset(fname, 'r')
    VarOut = _ncvar_to_dict(ncFile.variables[Rec])
    ncFile.close()
    return VarOut


def _determine_RAF(ncFile):
    if 'nimbus' in ncFile.Conventions.lower():
        coord0 = ncFile.coordinates.split(' ')[0]
        try:
            RAFrate = ncFile.variables[coord0].shape[1]
        except:
            RAFrate = 1
        RAFdim = ('sps' + str(RAFrate))
        return True, RAFrate, RAFdim
    else:
        return False, None, None


def _get_time(ncFile, isRAF, RAFrate=None, timevar=None):
    # Pull out the start time
    if timevar is not None:
        varname = timevar
    elif (timevar is None) & ('base_time' in ncFile.variables.keys()):
        varname = 'base_time'
    elif (timevar is None) & ('time' in ncFile.variables.keys()):
        varname = 'time'
    elif (timevar is None) & ('Time' in ncFile.variables.keys()):
        varname = 'Time'

    # Sanity check on the time variable in case a dummy time variable
    # is used - yes this really happens...Ugh
    if varname is not None:
        if np.array(ncFile.variables[varname][:]).min() < 0:
            varname = None

    if varname is not None:
        print("Using '%s' to make AWOT time variable" % varname)
        TimeSec = np.array(ncFile.variables[varname][:]).ravel()

        # Check if it is a high rate file and 2D - yep instances of this
        # out there as well...
        if (isRAF) & (RAFrate > 1) & (RAFrate not in ncFile.variables[varname].shape):
            Timehirate = np.linspace(
            TimeSec[0], TimeSec[-1], len(TimeSec) * RAFrate)
            TimeSec = Timehirate
    else:
        print("No time variable found, using StarTime to make AWOT time variable")
        StartTime = ncFile.StartTime
        length = len(ncFile.dimensions['Time'])
        # Create a time array
        TimeSec = np.linspace(StartTime, StartTime + length, length)
        
    try:
        time_units = ncFile.variables[varname].units
    except:
        time_units = _get_epoch_units()        

    # Now convert the time array into a datetime instance
    dtHrs = num2date(TimeSec, time_units)
    # Now convert this datetime instance into a number of seconds since Epoch
    TimeSec = date2num(dtHrs, _get_epoch_units())
    # Now once again convert this data into a datetime instance
    Time_unaware = num2date(TimeSec, _get_epoch_units())
    Time = {'data': Time_unaware, 'units': _get_epoch_units(),
            'title': 'Time', 'full_name': 'Time (UTC)'}
    return Time


def _make_data_dictionary(ncFile, name_map, isRAF, RAFrate=None):
    data = {}

    for var in name_map:
        if name_map[var] in ncFile.variables.keys():
            data[var] = _ncvar_to_dict(ncFile.variables[name_map[var]])
            if (isRAF) & (RAFrate in ncFile.variables[name_map[var]].shape):
                data[var]['data'] = np.array(ncFile.variables[name_map[var]][:]).ravel()
            try:
                mask = data[var]['data'].mask
            except:
                data[var]['data'] = np.ma.masked_array(data[var]['data'], mask=False)
        else:
            data[var] = None
    return data

#########################
#   NASA AMES Methods   #
#########################


def read_nasa_ames(filename, mapping_dict=None, platform=None):
    '''
    Read NASA AMES FFI 1001 formatted data files.
    The header tells all about the file. Find format here:
    https://espoarchive.nasa.gov/content/Ames_Format_Specification_v20#tth_sEc5.1
    Output variable names are controlled by mapping_dict and platform
    keywords. If nothing is chosen a direct copy of names in file will
    be output.
    Parameters
    ----------
    filename : str
        Name of file (long path okay).
    mapping_dict : dict
        Dictionary to use for mapping variable names. Use if provided.
        If None, use skip the dictionary and use names in file.
        Supersedes the platform keyword argument.
    platform : str
        If value is set and valid, use the default mapping
        dictionary provided for that platform.
    '''
    f = open(filename, 'r')

    # Retrieve header information
    hdr = _get_ames_header(f)

    # Read in the data from file
    junk = np.genfromtxt(filename, skiprows=int(hdr['NLHEAD']),
                         missing_values=hdr['VMISS'], filling_values=np.nan)

    # Grab a name map for data
    if mapping_dict is not None:
        name_map = mapping_dict
    else:
        # Check to see if platform is specified
        if platform is not None:
            name_map = _get_name_map(platform)
        else:
            name_map = _make_name_map_from_varlist(hdr['VNAME'])

    # Loop through the variables and pull data
    readfile = {}

    print(len(hdr['VMISS']))
    for jj, name in enumerate(hdr['VNAME']):
        readfile[name] = np.array(junk[:, jj] * hdr['VSCAL'][jj])
        readfile[name] = np.ma.masked_values(readfile[name], hdr['VMISS'][jj])

    # Now work out the time into a datetime instance
    # Get the date
    DateStr = str(hdr['DATE'][0]) + '-' + \
        str(hdr['DATE'][1]) + '-' + str(hdr['DATE'][2])
    StartTime = datetime(hdr['DATE'][0], hdr['DATE'][
                         1], hdr['DATE'][2], 0, 0, 0)
    TimeSec = np.array(readfile['time'][:]) + date2num(
        StartTime, units=_get_epoch_units())

    # Finally convert this back to a standard used by this package (Epoch)
    Time_unaware = num2date(TimeSec[:], units=_get_epoch_units())
    Time = Time_unaware  # .replace(tzinfo=pytz.UTC)

    del readfile['time']
    readfile['time'] = Time

    data = {}
    for varname in name_map:
        try:
            data[varname] = readfile[name_map[varname]]
        except:
            data[varname] = None

    # Calculate U,V wind if not present
    if 'Uwind' not in name_map:
        Uwind, Vwind = _winduv(data)
        # Add to the dictionary
        data['Uwind'] = Uwind
        data['Vwind'] = Vwind

    # Add some values to the data dictionary
    data['project'] = hdr['MNAME']
    data['platform'] = hdr['SNAME']
    data['flight_number'] = hdr['NCOM'][0].split(':')[1]

    return data


def _get_ames_header(f):
    '''
    NLHEAD : Number of header lines
    FFI : NASA AMES FFI format number
    ONAME : Originator/PI Name
    ORG : Name of organization
    SNAME : Instrument/platform name
    MNAME : Project/mission name
    IVOL : Current volume number (almost always 1)
    NVOL : Number of volumes for data (almost always 1)
    DATE : YYYY MM DD UTC begin date
    RDATE : Reduction/revision UTC date
    DX : Interval between successive values (data rate)
    XNAME : Name/Description of DX variable above
    NV : Number of primary variables in file
    VSCL : Scaling factor for each variable column
    VMISS : Missing value for each variable column
    VNAME : Name of first variable
    NSCOML : Number of special comment lines within header
    SCOM : Special comments about file/data, etc.
    NNCOML : Number of normal comment lines within header
    NCOM : Normal comments
    '''
    hdr = {}
    hdr['NLHEAD'], hdr['FFI'] = f.readline().split()

    # Check that the file is indeed NASA AMES 1001
    if hdr['FFI'] != '1001':
        print("Check file type, looks like it's not FFI 1001")
        return

    hdr['ONAME'] = f.readline().rstrip('\n')
    hdr['ORG'] = f.readline().rstrip('\n')
    hdr['SNAME'] = f.readline().rstrip('\n')
    hdr['MNAME'] = f.readline().rstrip('\n')
    hdr['IVOL'], hdr['NVOL'] = f.readline().split()
    yy1, mm1, dd1, yy2, mm2, dd2 = f.readline().split()
    hdr['DATE'] = (int(yy1), int(mm1), int(dd1))
    hdr['RDATE'] = (int(yy2), int(mm2), int(dd2))
    hdr['DX'] = f.readline().rstrip('\n')
    hdr['XNAME'] = f.readline().rstrip('\n')
    hdr['NV'] = int(f.readline().rstrip('\n'))
    vscl = f.readline().split()
    hdr['VSCAL'] = [int(x) for x in vscl]
    vmiss = f.readline().split()
    hdr['VMISS'] = [float(x) for x in vmiss]
    hdr['VNAME'] = ['time']
    for nvar in range(hdr['NV']):
        hdr['VNAME'].append(f.readline().rstrip('\n'))
    hdr['NSCOML'] = int(f.readline().rstrip('\n'))
    hdr['SCOM'] = []
    for nscom in range(hdr['NSCOML']):
        hdr['SCOM'].append(f.readline().rstrip('\n'))
    hdr['NNCOML'] = int(f.readline().rstrip('\n'))
    hdr['NCOM'] = []
    for nncom in range(hdr['NNCOML']):
        hdr['NCOM'].append(f.readline().rstrip('\n'))
    # Insert elements to account for time column
    hdr['VSCAL'].insert(0, 1)
    hdr['VMISS'].insert(0, np.nan)
    f.close()
    return hdr


def _latmos_name_map():
    '''Map out names used in SAFIRE/LATMOS Falcon data to AWOT'''
    name_map = {
        'time': 'time',
        'latitude': 'latitude : from GPS (degree)',
        'longitude': 'longitude : from GPS (degree)',
        'altitude': 'altitude : from GPS (meter)',
        'true_heading': 'platform_orientation : from INS (degree)',
        'temperature':
        'air_temperature : from deiced Rosemount sensor (Celsius)',
        'dewpoint_temperature': 'dew_point_temperature : from 1011B top dew-point hygrometer (Celsius)',
        'wind_spd': 'wind_speed : Attitudes and speed wrt ground from INS, air angles from radome, air speed from pitot (m/s)',
        'wind_dir': 'wind_from_direction : Attitudes and speed wrt ground from INS, air angles from radome, air speed from pitot (degree)',
        'relative_humidity': 'relative_humidity : from Aerodata sensor (%)',
        'mixing_ratio': 'humidity_mixing_ratio : from Aerodata sensor (gram/kg)',
        'pressure': 'air_pressure : from front sensor, corrected for the so-called static defect (hPa)',
        'roll': 'platform_roll_angle : from INS (degree)',
        'pitch': 'platform_pitch_angle : from INS (degree)',
        'aircraft_air_speed': 'platform_speed_wrt_air : from pitot (m/s)',
        'platform_ground_speed': 'platform_speed_wrt_ground : from GPS (m/s)',
        'platform_ground_speed2': 'platform_speed_wrt_ground : from INS (kt)',
        'aircraft_vert_accel': 'platform_acceleration_along_vertical_axis : from INS (meter second-2)',
        'altitude2': 'altitude : from INS (meter)',
        'mixing_ratio2': 'humidity_mixing_ratio : from top dew-point hygrometer (GE 1011B) (gram/kg)',
        'platform_course': 'platform_course : from INS (degree)',
        'platform_upward_ground_speed': 'upward_platform_speed_wrt_ground : from INS (m/s)',
        'platform_upward_ground_speed': 'upward_platform_speed_wrt_ground : from GPS (m/s)',
        'attack_angle': 'angle_of_attack : from sensor on the boom (degree)',
        'sideslip_angle': 'angle_of_sideslip : from sensor on the boom (degree)',
        'Uwind': 'eastward_wind : Attitudes and speed wrt ground from INS, air angles from radome, air speed from pitot (m/s)',
        'Vwind': 'northward_wind : Attitudes and speed wrt ground from INS, air angles from radome, air speed from pitot (m/s)',
        'air_vertical_velocity': 'upward_air_velocity : Attitudes and speed wrt ground from INS, air angles from radome, air speed from pitot (m/s)',
        'wind_dir2': 'wind_from_direction : Attitudes and speed wrt ground from INS, air angles from radome, air speed from pitot (degree)',
        'wind_spd2': 'wind_speed : Attitudes and speed wrt ground from INS, air angles from radome, air speed from pitot (m/s)',
    }
    return name_map


def _und_citation_name_map():
    '''Map out names used in UND Citation data to AWOT'''
    name_map = {
        'time': 'time',
        'latitude': 'POS_Lat',
        'longitude': 'POS_Lon',
        'altitude': 'POS_Alt',
        'pressure_altitude': 'Press_Alt',
        'true_heading': 'POS_Head',
        'track': 'POS_Trk',
        'temperature': 'Air_Temp',
        'dewpoint_temperature': 'DEWPT',
        'theta': 'Pot_Temp_T1',
        'wind_spd': 'Wind_M',
        'wind_dir': 'Wind_D',
        'mixing_ratio': 'MixingRatio',
        'frost_point_temp': 'FrostPoint',
        'lwc': 'King_LWC_ad',
        'twc': 'Nev_TWC',
        'Conc_2DC': '2-DC_Conc',
        'Dmean_2DC': '2-DC_MenD',
        'Dvol_2DC': '2-DC_VolDia',
        'Deff_2DC': '2-DC_EffRad',
        'Conc_CPC': 'CPCConc',
        'air_vertical_velocity': 'Wind_Z',
        'turb': 'TURB',
    }
    return name_map

######################
#   Shared methods   #
######################


def _make_name_map_from_varlist(varlist):
    '''
    Use the passed list of variable names to construct a list.
    Duplicates all variable names.
    '''
    name_map = {}
    for var in varlist:
        name_map[var] = var
    return name_map


def _winduv(data):
    """
    Calculate the horizontal windcomponents (u, v)
    from wind angle and speed.
        U wind : positive blowing towards east
        V wind : positive blowing towards north
    """
    try:
        U = {}
        U['data'] = -np.cos(np.radians(data['wind_dir']['data'][:])) * data['wind_spd']['data'][:]
        U['units'] = data['wind_spd']['units']
        U['title'] = "U Wind"
        U['full_name'] = "Zonal Wind, positive blowing towards east"
    except:
        U = None
    try:
        V = {}
        V ['data']= -np.sin(np.radians(data['wind_dir']['data'][:])) * data['wind_spd']['data'][:]
        V['units'] = data['wind_spd']['units']
        V['title'] = "V Wind"
        V['full_name'] = "Meridional Wind, positive blowing towards east"
    except:
        V = None
    return U, V
