"""
awot.io.read_p3_radar
=========================

A group of scripts to read data collected by the NOAA P-3 aircraft.
Supports both tail Doppler and lower fuselage radars. 

Author:
    08 Jan 2014 - Created by Nick Guy, NOAA/NSSL/WRDD, NRC.
                  Converted NCL functions to Python

"""
# NOTES:: This has only been tested with DYNAMO data files, versions
#         may change and another function may be needed.  
#-------------------------------------------------------------------
# Load the needed packages
from netCDF4 import Dataset
import numpy as np
import pytz
from datetime import datetime
import pyart
#-------------------------------------------------------------------
# Begin methods
#####################################################################
# TDR file methods #
######################
def read_windsyn_tdr_netcdf(fname, mapping=None):
    """Read in a NetCDF data file containing a pseudo-dual-Doppler
    analysis of NOAA P-3 tail radar data from the windsyn program.  
    
    See the src directory for more details on the windsyn program.
        
    Parameters::
    ----------
    fname : str
        Long path filename [string]
    mapping : str
            Mapping dictionary to use for pulling data in
            If not set, equivalent to None, then the basic variables
            are input.  
            However, if 'extended' is selected then 
            additional variables are mapped.
        
    Output::
    ------
    data : Dictionary of the following values
    
    metadata : Dictionary of global attributes in file
    longitude : dict
        Longitude array [degrees]
    latitude : dict
        Latitude array [degrees]
    height : dict
        Height array [km]
    fields : Dictionary of variables in file
        reflectivity : float
            Radar Reflectivity [dBZ]
        Uwind : float
            Wind along aircraft longitudinal axis wind [m/s]
        Vwind : float
            Wind perpendicular to aircraft longitudinal axis wind [m/s]
        Wwind : float
            Vertical wind component [m/s]
        divergene : float
            Divergence x 10^3 [1/s]
        term_fall_speed : float
            Terminal fall speed [m/s]
        time_diff : float
            Time difference from nominal [s]
        datetime_start : datetime object
            HHNN (Hour Minute) module start time
        datetime_end : datetime object
            HHNN (Hour Minute) module end time
    """
     
    # Read in NetCDF
    ncFile = Dataset(fname,'r')
    ncvars = ncFile.variables
    
    # Grab the metadata stored in global attributes as a dictionary
    metadata = ncFile.__dict__
    
    # Find the module start information
    # Year
    yys = ncFile.variables['YYstart'][:]
    # Month	
    mms = ncFile.variables['MMstart'][:]
    # Day
    dds = ncFile.variables['DDstart'][:]
    # Hour
    hhs = ncFile.variables['HHstart'][:]
    # Minute
    nns = ncFile.variables['MIstart'][:]
    # Second
    sss = ncFile.variables['SSstart'][:]
    #--------------------------------
    # Find the module end information
    # Year
    yye = ncFile.variables['YYend'][:]
    # Month	
    mme = ncFile.variables['MMend'][:]
    # Day
    dde = ncFile.variables['DDend'][:]
    # Hour
    hhe = ncFile.variables['HHend'][:]
    # Minute
    nne = ncFile.variables['MIend'][:]
    # Second
    sse = ncFile.variables['SSend'][:]
     
    # Make a Datetime instance to pass
    utc = pytz.utc
    dt_start = datetime(yys,mms,dds,hhs,nns,sss,0)#,tzinfo=utc)
    dt_end = datetime(yye,mme,dde,hhe,nne,sse,0)#,tzinfo=utc)
    #--------------------------------------
    # Pull in the Lat / Lon information
    Lon = _ncvar_to_dict(ncvars['Lon'])
    Lat = _ncvar_to_dict(ncvars['Lat'])
    
    # Pull in the Height 
    Ht = _ncvar_to_dict(ncvars['Height'])
    #---------------------------------------
    
    # Add fields to their own dictionary
    fields = {}
    
    # Grab a name map for NOAA P-3 TDR data
    name_map = _get_tdr_name_map(mapping=mapping)
    
    # Loop through the variables and pull data
    for varname in name_map:
        fields[varname] = _ncvar_to_dict(ncvars[name_map[varname]])
    
    # Mask bad data values
    #badval = float(ncFile.MissingValue)
    #np.ma.masked_values()
    
    # Pull out specific information on platform/instrument   
    try:
        ncFile.Platform
        platform = ncFile.title
    except:
        platform = 'p3'    

    try:
        ncFile.Instrument
        instrument = ncFile.instrument
    except:
        instrument = 'tdr_grid'
    
    # Create a dictionary to transfer the data
    data = {'metadata' : metadata,
            'longitude': Lon,
            'latitude': Lat,
            'height': Ht,
            'fields': fields,
            'datetime_start' : dt_start,
            'datetime_end' : dt_end,
            'platform' : platform,
            'instrument' : instrument
            }
    
    ncFile.close()

    return data

###############################################

def tdr_grid_variable(fname,Rec):
    """Read in a variable from a gridded NetCDF file 
    containing NOAA P-3 tail Doppler radar fields.
    
    Parameters::
    ----------
    fname : str
        Long path filename 
    VarName : str
        Variable name to access
        
    Usage::
    -----
        VarOut = tdr_grid_variable(fname, Rec)
    """
    # Read the NetCDF
    ncFile = Dataset(fname,'r')

    # Get the variable of interest
    VarOut = ncFile.variables[VarName][:]
    
    ncFile.close()
    return VarOut
    
#######################################################

def read_windsyn_binary(fname, platform=None, instrument=None, radar_num=None):
    """Read in unformatted Fortran binary data files created by the windsyn program.
        
    Parameters::
    ----------
    fname : str
        Long path filename [string]
    platform : str
        Name of platform that instrument is aboard
    instrument : str
        Name of instrument taking observations
    radar_num : int
        Number of radar from which to pull datetime start/finish information
        
    Output::
    ------
    data : Dictionary of the following values
    
    Notes::
    -----
    The windsyn program stores the pseudo-dual(or multi)-Doppler analysis 
    data in two files with header information (.hdr) and data (.dpw).
    Both files are written using Big Endian order.
    
    The header file is written with the following Fortran call:
    Open (3,Err=999,File=File,Iostat=Ierr,Form='UNFORMATTED', Convert='Big_Endian')
    
    The data solution file is written with the following Fortran call:
    Open (3, Err=999, File=File, Iostat=Ierr, Access='Direct',
     #      Recl=Imax * 4, Convert='Big_Endian')
     
     Older versions of the .dpw files may have only 6 variables.  In this case, 
     it may be desirable to re-run the windsyn analysis.
    """
    
    # Store the platform information or set with default guess 
    if platform is None:
        print "Guessing that platform is P-3, this can be set with platform keyword"
        platform = 'p3'  

    if instrument is None:
        print "Guessing that instrument is TDR, this can be set with platform keyword"
        instrument = 'tdr_grid'
        
    # Read in the .hdr file
    hdr_dict = _construct_windsyn_hdr_dict(fname)
    
    # Set dimension
    Imax = hdr_dict['Imax']
    Jmax = hdr_dict['Jmax']
    Kmax = hdr_dict['Kmax']
    ##########
    
    # Retrieve data from the .dpw file
    fi = open(fname+'.dpw', 'rb')
    
    U = np.fromfile(fi, dtype='>f', 
        count=(Imax * Jmax * Kmax)).reshape(Kmax, Jmax, Imax)
    V = np.fromfile(fi, dtype='>f', 
        count=(Imax * Jmax * Kmax)).reshape(Kmax, Jmax, Imax)
    W = np.fromfile(fi, dtype='>f', 
        count=(Imax * Jmax * Kmax)).reshape(Kmax, Jmax, Imax)
    Div = np.fromfile(fi, dtype='>f', 
        count=(Imax * Jmax * Kmax)).reshape(Kmax, Jmax, Imax)
    dBZ = np.fromfile(fi, dtype='>f', 
        count=(Imax * Jmax * Kmax)).reshape(Kmax, Jmax, Imax)
    Vt = np.fromfile(fi, dtype='>f', 
        count=(Imax * Jmax * Kmax)).reshape(Kmax, Jmax, Imax)
    Tdif = np.fromfile(fi ,dtype='>f', 
        count=(Imax * Jmax * Kmax)).reshape(Kmax, Jmax, Imax)
    Tave = np.fromfile(fi ,dtype='>f', 
        count=(Imax * Jmax * Kmax)).reshape(Kmax, Jmax, Imax)
    
    fi.close()
    
    # Create dictionary definitions
    Lon, Lat = _get_lon_lat_from_header(hdr_dict)
    
    # Create the fields
    fields = {}
    # Read the Reflectivity
    fields['dBZ'] = _windsyn_var_to_dict(dBZ, hdr_dict, 
                    units='dBZ', 
                    long_name='Reflectivity')
    
    # Read the U wind (Along aircraft longitudinal axis)
    fields['Uwind'] = _windsyn_var_to_dict(U, hdr_dict, 
                      units=r'm s$^{-1}$', 
                      long_name='East-West Relative Wind Velocity')
    
    # Read the V wind (Perpendicular to aircraft longitudinal axis)
    fields['Vwind'] = _windsyn_var_to_dict(V, hdr_dict, 
                      units=r'm s$^{-1}$', 
                      long_name='North-South Relative Wind Velocity')
    
    # Read the vertical wind
    fields['Wwind'] = _windsyn_var_to_dict(W, hdr_dict, 
                      units=r'm s$^{-1}$', 
                      long_name='Vertical Velocity')
    
    # Try reading some variables that are not always included
    # Read the divergence
    try:
        fields['divergence'] = _windsyn_var_to_dict(Div*1000., hdr_dict, 
                              units=r'X10$^{3}$ s$^{-1}$',
                              long_name='Divergence')
    except:
        fields['divergence'] = None
    
    # Read the terminal velocity
    try:
        fields['term_fall_speed'] = _windsyn_var_to_dict(Vt, hdr_dict, 
                                    units=r'm s$^{-1}$', 
                                    long_name='Terminal Fall Velocity')
    except:
        fields['term_fall_speed'] = None
    
    # Read the time difference
    try:
        fields['time_diff'] = _windsyn_var_to_dict(Tdiff, hdr_dict, 
                              units=r's$', 
                              long_name='Time Difference from nominal')
    except:
        fields['time_diff'] = None
    
    # Create a dictionary to transfer the data
    data = {'metadata' : _get_metadata_from_header(hdr_dict),
            'longitude': Lon,
            'latitude': Lat,
            'height': _get_height_from_header(hdr_dict),
            'fields': fields,
            'datetime_start': _get_datetime_start_from_header(hdr_dict),
            'datetime_end': _get_datetime_start_from_header(hdr_dict),
            'platform': platform,
            'instrument': instrument
            }

    return data
    
#######################################################
# Lower fuselage gridded methods #
##################################
    
def read_lf_grid(fname):
    """Read in a NetCDF data file containing gridded NOAA P-3 
    lower fuselage radar data.
        
    Parameters::
    ----------
    fname : str
        Long path filename [string]
        
    Output::
    ----------
    data : Dictionary of the following values
    
    metadata: Dictionary of metadata in file
    longitude : float
        Longitude array [degrees]
    latitude : float
        Latitude array [degrees]
    height : float
        Height array [km]
    fields: Dictionary of variables
        reflectivity : float
            Radar Reflectivity [dBZ]
    datetime_start : datetime object
        HHNN (Hour Minute) module start time
    datetime_end : datetime object
        HHNN (Hour Minute) module end time
    """
    # Read in NetCDF
    ncFile = Dataset(fname,'r')
    ncvars = ncFile.variables
    
    # Grab the metadata stored in global attributes as a dictionary
    metadata = ncFile.__dict__
    
    # Find the module start information
    # Year
    yys = ncFile.variables['YYstart'][:]
    # Month	
    mms = ncFile.variables['MMstart'][:]
    # Day
    dds = ncFile.variables['DDstart'][:]
    # Hour
    hhs = ncFile.variables['HHstart'][:]
    # Minute
    nns = ncFile.variables['MIstart'][:]
    # Second
    sss = ncFile.variables['SSstart'][:]
    #--------------------------------
    # Find the module end information
    # Year
    yye = ncFile.variables['YYend'][:]
    # Month	
    mme = ncFile.variables['MMend'][:]
    # Day
    dde = ncFile.variables['DDend'][:]
    # Hour
    hhe = ncFile.variables['HHend'][:]
    # Minute
    nne = ncFile.variables['MIend'][:]
    # Second
    sse = ncFile.variables['SSend'][:]

    # Make a Datetime instance to pass
    utc = pytz.utc
    dt_start = datetime(yys,mms,dds,hhs,nns,sss,0)#,tzinfo=utc)
    dt_end = datetime(yye,mme,dde,hhe,nne,sse,0)#,tzinfo=utc)
    #---------------------------------------------------
    # Pull in the Lat / Lon information
    Lon = _ncvar_to_dict(ncvars['Lon'])
    Lat = _ncvar_to_dict(ncvars['Lat'])
    
    # Pull in the Height of aircraft
    Ht_level = _ncvar_to_dict(ncvars['Zlev'])
    #---------------------------------------------------
    
    fields = {}
    # Read the Reflectivity
    fields['dBZ'] = _ncvar_to_dict(ncvars['dBZ'])
    
    #badval = float(ncFile.MissingValue)
    #np.ma.masked_values(fields['dBZ']['data'], badval)
    
    # Pull out specific information on platform/instrument   
    try:
        ncFile.Platform
        platform = ncFile.title
    except:
        platform = 'p3'    

    try:
        ncFile.Instrument
        instrument = ncFile.instrument
    except:
        instrument = 'lf'
    
    # Create a dictionary to transfer the data
    data = {'metadata' : metadata,
            'longitude': Lon,
            'latitude': Lat,
            'height': Ht_level,
            'fields' : fields,
            'datetime_start' : dt_start,
            'datetime_end' : dt_end,
            'platform' : platform,
            'instrument' : instrument
            }
    
    ncFile.close()

    return data
    
#######################################################
# Sweep file methods #
######################

def read_tdr_sweep(fname):
    """Read in a tail radar sweep file using the PyArt interface.
    Both cfradial and raw Sigmet formats are accepted.
        
    Parameters::
    ----------
    fname : str
        Long path filename [string]
        
    Output::
    ----------
    data : Dictionary of the following values
    """
    data = pyart.io.read(fname)
    
    return data

def read_lf_sweep(fname):
    """Read in a lower fuselage radar sweep file using the PyArt interface.
    Both cfradial and raw Sigmet formats are accepted.
        
    Parameters::
    ----------
    fname : str
        Long path filename [string]
        
    Output::
    ----------
    data : Dictionary of the following values
    """
    data = pyart.io.read(fname)
    
    return data
    
###################################################
# Get methods using the header dictionary #
###################################################

def _get_height_from_header(hdr):
    """Calculated the height array given windsyn binary header"""
#    Height = {'data': np.arange(hdr['Z0'], hdr['Kmax'], hdr['Sz']),
    Height = {'data': np.arange(0, hdr['Kmax'] * hdr['Sz'], hdr['Sz']),
              'units': 'km',
              'long_name': 'Altitude'}
    if np.isscalar(Height['data']):
        Height['data'] = np.array(Height['data'])
        Height['data'].shape = (1, )
        
    return Height
    
def _get_x_y_from_header(hdr):
    """Calculated the X and Y arrays given windsyn binary header"""
    X = {'data': np.arange(0., hdr['Imax'] * hdr['Sx'], hdr['Sx']),
        'units': 'km',
        'long_name': 'East-West Distance'}
    Y = {'data': np.arange(0., hdr['Jmax'] * hdr['Sy'], hdr['Sy']),
        'units': 'km',
        'long_name': 'North-South Distance'}
    
    return X, Y
    
def _get_lon_lat_from_header(hdr):
    """Calculated the Lon and Lat arrays given windsyn binary header"""
    #Approximate radius of Earth
    Re=6371.
    
    X, Y = _get_x_y_from_header(hdr)
    
    # Close approximation of longitude and latitude
    Lon = {'data': (( X['data'] / Re ) * ( 180. / np.pi )) + float(hdr['Olon']),
          'units': 'degrees',
          'long_name': 'Longitude'}
    Lat = {'data': (( Y['data'] / Re ) * ( 180. / np.pi )) + float(hdr['Olat']),
              'units': 'degrees',
              'long_name': 'Latitude'}
    return Lon, Lat
    
def _get_datetime_start_from_header(hdr, radar_num=None):
    """Calculate the start time datetime instance from binary header"""
    if radar_num is None:
        radar_num = 1
    tstart = hdr['Itime_Limits'+str(radar_num)][0:5].copy()
    if tstart[0] < 50:
        tstart[0] = tstart[0] + 2000
    
    dt_start = datetime(tstart[0], tstart[1], tstart[2], tstart[3], tstart[4], 0)
    
    return dt_start
    
def _get_datetime_end_from_header(hdr, radar_num=None):
    if radar_num is None:
        radar_num = 1
    """Calculate the end time datetime instance from binary header"""
    tend = hdr['Itime_Limits'+str(radar_num)][6:11].copy()
    if tend[0] < 50:
        tend[0] = tend[0] + 2000
        
    dt_end = datetime(tend[0], tend[1], tend[2], tend[3], tend[4], 0)
     
    return dt_end
    
def _get_metadata_from_header(hdr):
    """Create a metadata dictionary with data details to be saved"""
    metadata = {'Flight_ID': hdr['Flid1'],
		        'Obrien_div_correction': hdr['Iobr'],
		        'number_smooth': hdr['Nsmth'],
		        'dual_Doppler_Analysis': hdr['Istyle'],
		        'Integration_direction': hdr['Idir_Int'],
		        'dBZ_calculation': hdr['IdBZ_parm'],
		        'Nmosm': hdr['Nmosm'],
		        'number_radars': hdr['Nrdrs'],
		        'creation_date': datetime.now(),
		        'source_file': hdr['filename'],
		        'title': "Windsyn computatons of P-3 tail radar data"
		        }
	 
    return metadata
    
###########################
# Create Variable methods #
###########################

def _ncvar_to_dict(ncvar):
    """ Convert a NetCDF Dataset variable to a dictionary. 
    Appropriated from PyArt package
    """
    d = dict((k, getattr(ncvar, k)) for k in ncvar.ncattrs())
    d['data'] = ncvar[:]
    if np.isscalar(d['data']):
        # netCDF4 1.1.0+ returns a scalar for 0-dim array, we always want
        # 1-dim+ arrays with a valid shape.
        d['data'] = np.array(d['data'])
        d['data'].shape = (1, )
    return d
    
def _get_tdr_name_map(mapping=None):
    '''Map out names used in NOAA P-3 TDR NetCDF files to AWOT'''
    name_map = {}
    
    # Radar reflectivity
    name_map['dBZ'] = 'dBZ'
    # U wind component (Along aircraft longitudinal axis)
    name_map['Uwind'] = 'U'
    # V wind component (Perpendicular aircraft longitudinal axis)
    name_map['Vwind'] = 'V'
    # Vertical wind component
    name_map['Wwind'] = 'W'
        
    if mapping == 'extended':
        name_map['divergence'] = 'Div'
        name_map['term_fall_speed'] = 'Vt'
        name_map['time_diff'] = 'Tdiff'
    return name_map
    
def _get_lf_name_map():
    '''Map out names used in NOAA P-3 LF NetCDF files to AWOT'''
    name_map = {}
    
    # Radar reflectivity
    name_map['dBZ'] = 'dBZ'

    
def _windsyn_var_to_dict(wsvar, hdr_dict, units=None, long_name=None):
    """ Convert a variable read in from the windsyn 
    binary files to a dictionary.
    """
    d = dict()
    
    # Apply a mask to data for "bad" values
    wsvar = np.ma.masked_equal(wsvar, hdr_dict['Baddata'])
    
    d['data'] = wsvar
    if np.isscalar(d['data']):
        # netCDF4 1.1.0+ returns a scalar for 0-dim array, we always want
        # 1-dim+ arrays with a valid shape.
        d['data'] = np.array(d['data'])
        d['data'].shape = (1, )
        
    d['units'] = units
    d['long_name'] = long_name
    d['_FillValue'] = hdr_dict['Baddata']
    
    return d

def _construct_windsyn_hdr_dict(filename):
    """Read through a windsyn generated header file and construct an
    output dictionary with data"""
    fi = open(filename+'.hdr', 'rb')
    
    # Skip through the initial padding
    fi.seek(4)
    
    # Now read individual values 
    # Adds up to 222 bytes = (4 byte pad front + 214 byte vars + 4 byte pad end)
    Real_Time = np.fromfile(fi, dtype='>c', count=24)
    Imax, Jmax, Kmax = np.fromfile(fi, dtype='>i4', count=3)
    Sx, Sy, Sz = np.fromfile(fi, dtype='>f', count=3)
    Olat, Olon, Z0 = np.fromfile(fi, dtype='>f', count=3)
    Nrdrs, Nmosm = np.fromfile(fi, dtype='>i4', count=2)
    Su, Sv = np.fromfile(fi, dtype='>f', count=2)
    IdBZ_Parm, Idir_Int, Istyle = np.fromfile(fi, dtype='>i4', count=3)
    Nsmth, Iobr = np.fromfile(fi, dtype='>i4', count=2)
    fname = np.fromfile(fi, dtype='>c', count=50)
    Baddata, Vt_Snow, Vt_Rain = np.fromfile(fi, dtype='>f', count=3)
    ThresH, ThresV = np.fromfile(fi, dtype='>f', count=2)
    IsmTyp, Ihole_Fill, Klim = np.fromfile(fi, dtype='>i4', count=3)
    Iw_at_top, I_w_0 = np.fromfile(fi, dtype='>i4', count=2)
    Top_Hgt = np.fromfile(fi, dtype='>f', count=1)
    junk = np.fromfile(fi, dtype='>i4', count=5)
    pad = np.fromfile(fi, dtype='>i4', count=1)
    
    # Establish dictionary
    names = ['Real_Time', 'Imax', 'Jmax', 'Kmax', 'Sx', 'Sy', 'Sz',
             'Olat', 'Olon', 'Z0', 'Nrdrs', 'Nmosm', 'Su', 'Sv',
             'IdBZ_parm', 'Idir_Int', 'Istyle', 'Nsmth', 'Iobr',
             'filename', 'Baddata', 'Vt_Snow', 'Vt_Rain', 'ThresH',
             'ThresV', 'IsmTyp', 'Ihole_Fill', 'Klim', 
             'Iw_at_top', 'I_w_0', 'Top_Hgt',
             ]
             
    values = [np.ndarray.tostring(Real_Time), Imax, Jmax, Kmax, Sx, Sy, Sz,
             Olat, Olon, Z0, Nrdrs, Nmosm, Su, Sv,
             IdBZ_Parm, Idir_Int, Istyle, Nsmth, Iobr,
             np.ndarray.tostring(fname), Baddata, Vt_Snow, Vt_Rain, ThresH,
             ThresV, IsmTyp, Ihole_Fill, Klim, 
             Iw_at_top, I_w_0, Top_Hgt,
             ]
             
    hdr_dict = dict(zip(names, values))
    
    # Loop through the multiple radars if there are more than one
    # Adds up to 319 bytes = ( 4 byte pad front + 311 byte vars + 4 byte pad end)
    Rnames = []
    Rvalues = []
    if Nrdrs > 1:
        for ii in range(Nrdrs):
            junk = np.fromfile(fi, '>i4', count=1)
            Real_Time = np.fromfile(fi, dtype='>c', count=24)
            Nameif = np.fromfile(fi, dtype='>c', count=50)
            Imax, Jmax, Kmax = np.fromfile(fi, dtype='>i4', count=3)
            Sx, Sy, Sz = np.fromfile(fi, dtype='>f', count=3)
            Flid = np.fromfile(fi, dtype='>c', count=8)
            Olat, Olon, Z0 = np.fromfile(fi, dtype='>f', count=3)
            Itime_Limits = np.fromfile(fi, dtype='>i4', count=12)
            Nmosm = np.fromfile(fi, dtype='>i4', count=1)
            Init_Time = np.fromfile(fi, dtype='>i4', count=6)
            Su, Sv = np.fromfile(fi, dtype='>f', count=2)
            Project = np.fromfile(fi, dtype='>c', count=16)
            Rotcr, Eloff, Cant = np.fromfile(fi, dtype='>f', count=3)
            Namdf = np.fromfile(fi, dtype='>c', count=50)
            Xniq, Rlat, Rlon = np.fromfile(fi, dtype='>f', count=3)
            Type = np.fromfile(fi, dtype='>c', count=3)
            RangeMax = np.fromfile(fi, dtype='>f', count=1)
            Alt_flag = np.fromfile(fi, dtype='>i4', count=1)
            junk = np.fromfile(fi, dtype='>i4', count=2)
            pad = np.fromfile(fi, dtype='>i4', count=1)

            Rnames = ['Real_Time'+str(ii+1), 'Nameif'+str(ii+1),
                    'Imax'+str(ii+1), 'Jmax'+str(ii+1), 'Kmax'+str(ii+1),
                    'Sx'+str(ii+1), 'Sy'+str(ii+1), 'Sz'+str(ii+1),
                    'Flid'+str(ii+1),
                    'Olat'+str(ii+1), 'Olon'+str(ii+1), 'Z0'+str(ii+1), 
                    'Itime_Limits'+str(ii+1), 'Nmosm'+str(ii+1), 'Init_Time'+str(ii+1),
                    'Su'+str(ii+1), 'Sv'+str(ii+1), 'Project'+str(ii+1),
                    'Rotcr'+str(ii+1), 'Eloff'+str(ii+1), 'Cant'+str(ii+1), 
                    'Namdf'+str(ii+1), 'Xniq'+str(ii+1),
                    'Rlat'+str(ii+1), 'Rlon'+str(ii+1), 'Type'+str(ii+1), 
                    'RangeMax'+str(ii+1), 'Alt_flag'+str(ii+1),
                    ]
              
            Rvalues = [np.ndarray.tostring(Real_Time), np.ndarray.tostring(Nameif),
                     Imax, Jmax, Kmax, Sx, Sy, Sz,
                     np.ndarray.tostring(Flid),
                     Olat, Olon, Z0, Itime_Limits, Nmosm,
                     Init_Time, Su, Sv,
                     np.ndarray.tostring(Project),
                     Rotcr, Eloff, Cant,
                     np.ndarray.tostring(Namdf),
                     Xniq, Rlat, Rlon,
                     np.ndarray.tostring(Type),
                     RangeMax, Alt_flag,
                     ]
            
            # Add the values for each radar to the header dictionary
            for var, val in zip(Rnames, Rvalues):
                hdr_dict[var] = val
    
    fi.close()
    
    return hdr_dict