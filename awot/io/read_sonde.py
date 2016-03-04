import numpy as np
import matplotlib.patches as patches
from .common import _build_dict


def read_sounding_data(filePath):

    '''
    method to retrieve standard sounding data from noaa radiosonde
    and U of Wyoming repository.

    Inputs:

    file path
    # ============

    output:

    # Dictionary with variable for use in plotting

    '''

    data = filePath
    fp = open(data, 'r')
    lines = fp.readlines()
    header = lines[1]
    fp.close()
    fp = open(data, 'r')
    sounding_data = fp
    p, h, T, Td, RH, MIXR, wd, ws = np.genfromtxt(
        sounding_data, skip_header=8, usecols=range(0, 8),
        dtype=float, delimiter=None, autostrip=True,
        missing_values='-9999.00', unpack=True, usemask=True)

    u = -ws*np.sin(np.radians(wd))
    v = -ws*np.cos(np.radians(wd))
    # RH = tC._dewpoint_to_RH(T+273.15, Td+273.15)

    mask = T.mask
    T = T[~mask]
    TD = Td[~mask]
    P = p[~mask]
    H = h[~mask]
    RH = RH[~mask]

    mask = u.mask
    Uwind = u[~mask]
    Vwind = v[~mask]

    data = dict()
    data['metadata'] = header
    data['temperature'] = _build_dict(T,
         'c', 'Temperature of ambient air', 'Temperature')
    data['dewpoint'] = _build_dict(TD,
        'c', 'Dewpoint temperature of ambient air', 'Dewpoint_Temperature')
    data['presssure'] = _build_dict(P,
        'hPa', 'Pressure of ambient air', 'Pressure')
    data['relative_humidity'] = _build_dict(RH,
        '%', 'Relative Humidity of ambient air’, ‘Relative_Humidity')
    data['u_component'] = _build_dict(Uwind,
        'm/s', 'u component of wind', 'U_component')
    data['v_component'] = _build_dict(Vwind,
        'm/s', 'v component of wind', 'V_component')
    data['height'] = _build_dict(H,
        'm', 'Geometric Height in meters', 'Height')
    data['data_format'] = 'radioSonde'

    fp.close()
    return data


def read_dropsonde_data(filePath, split_file=True):

    '''
    method to retrieve dropsonde data

    Inputs:

    file path
    #============

    output:

    #Dictionary with variable for use in plotting

    '''
    if split_file:
        data = split_cls_file(filePath)
        sounding_data = data
    else:
        data = open(filePath, 'r')
        lines = data.readlines()
        header = lines[1]
        data.close()
        data = open(data, 'r')
        sounding_data = data

    p, T, Td, RH, u, v, h = np.genfromtxt(
        sounding_data, skip_header=15, usecols=(1, 2, 3, 4, 5, 6, 14),
        dtype=float, missing_values='9999.0', unpack=True, usemask=True)

    # mask incoming T and dewpoint data

    T = np.ma.masked_greater_equal(T, 999.0)
    Td = np.ma.masked_greater_equal(Td, 999.0)
    RH = np.ma.masked_greater_equal(RH, 999.0)
    height = np.ma.masked_greater_equal(h, 99999.0)

    mask = T.mask
    T = T[~mask]
    TD = Td[~mask]
    P = p[~mask]
    H = h[~mask]
    RH = RH[~mask]

    mask = u.mask
    Uwind = u[~mask]
    Vwind = v[~mask]

    data = dict()
    data['metadata'] = header
    data['temperature'] = _build_dict(T,
        'c', 'Temperature of ambient air', 'Temperature')
    data['dewpoint'] = _build_dict(TD,
        'c', 'Dewpoint temperature of ambient air', 'Dewpoint Temperature')
    data['presssure'] = _build_dict(P, 'hPa',
        'Pressure of ambient air', 'Pressure')
    data['relative_humidity'] = _build_dict(RH,
        '%', 'Relative Humidity of ambient air’, ‘Relative Humidity')
    data['u_component'] = _build_dict(Uwind,
        'm/s', 'u component of wind', 'U component')
    data['v_component'] = _build_dict(Vwind,
        'm/s', 'v component of wind', 'V component')
    data['height'] = _build_dict(H,
        'm', 'Geometric Height in meters', 'Height')
    data['data_format'] = 'dropsonde'

    return data

def find_headers(filename):
    
    '''
    Finds the location of the headers of each data instance inside the CLS File
    
    Inputs:
    
    file path
    #============
    
    output:
    
    #list of locations
    
    
    '''
    headerID = 'Data Type:'
    
    f = open(filename, 'r')
    header_loc = []
    header_count = 0
    for linenum, line in enumerate(f):
        if  (headerID) in line:
            header_count +=1
            header_loc.append(linenum)
    
    return header_count, header_loc




def count_lines(filename):
    '''
        Returns the total number of lines in the file
        
        Inputs:
        
        file path
        #============
        
        output:
        
        #float
        '''
    f = open(filename, 'r')
    nlines = len(f.readlines())
    f.close()
    return nlines





def _get_header(f):
    '''
        data_type : Number of header lines
        project : NASA AMES FFI format number
        site : Originator/PI Name
        site_type : Name of organization
        site_id : Instrument/platform name
        lon0_dms : Project/mission name
        lat0_dms : Current volume number (almost always 1)
        lat0_dec : Number of volumes for data (almost always 1)
        alt : YYYY MM DD UTC begin date
        datetime0 : Reduction/revision UTC date
        launch_reference : Interval between successive values (data rate)
        sonde_id : Name/Description of DX variable above
        op_comments : Number of primary variables in file
        proc_comments : Scaling factor for each variable column
        varnames : Missing value for each variable column
        varunits : Name of first variable
        '''
    hdr = {}
    #List index out of range? See bototm of notebook
    hdr['data_type'] = f.readline().rstrip('\n').split(":")[1].strip()
    
    hdr['project'] = f.readline().rstrip('\n').split(":")[1].strip()
    hdr['site'], hdr['site_type'], hdr['site_id']  = f.readline().rstrip('\n').split(":")[1].split()
    lonlatalt= f.readline().rstrip('\n').split(":")[1].strip()
    hdr['lon0_dms'], hdr['lat0_dms'], hdr['lon0_dec'], hdr['lat0_dec'], hdr['alt'] = lonlatalt.split(',')
    # Turn some values to float
    hdr['lon0_dec'] = float(hdr['lon0_dec'])
    hdr['lat0_dec'] = float(hdr['lat0_dec'])
    hdr['alt'] = float(hdr['alt'])
    release_time = f.readline().rstrip('\n').split(":")
    yyyy, mm, dd, hh = release_time[1].strip().split(',')
    nn = release_time[2].strip()
    ss = release_time[2].strip()
    dstring = "%s-%s-%sT%s:%s:%s" % (int(yyyy), int(mm), int(dd), int(hh), int(nn), int(ss))
    hdr['datetime0'] = datetime(int(yyyy), int(mm), int(dd), int(hh), int(nn), int(ss))
    hdr['launch_reference'] = f.readline().rstrip('\n').split(":")[1].strip() # Not complete - Throwaway
    hdr['sonde_id'] = f.readline().rstrip('\n').split(":")[1].strip()
    hdr['op_comments'] = f.readline().rstrip('\n').split(":")[1].strip()
    hdr['proc_comments'] = f.readline().rstrip('\n').split(":")[1].strip()
    f.readline() # Empty line '/'
    f.readline() # Empty line '/'
    f.readline() # This is the Nominal Release Time - Repeate datetime0
    hdr['varnames'] = f.readline().rstrip('\n').strip().split()
    hdr['varunits'] = f.readline().rstrip('\n').strip().split()
    f.readline() # Divider line
    
    return hdr


def read_cls_file(filename):
    '''
    #Method designed to read and seperate CLS Files from Noaa p3 Dropsonde data.
    
    Input:
    #FilePath
    
    Output:
    Dictionary containing file instances as sorted by date of launch and variable recorded.
    
        
    '''
def read_cls_file(filename):
    # Find the total number of lines
    nfilelines = count_lines(filename)
    # Find the numbers of headers
    nhdr, hdrloc = find_headers(filename)
    # Calculate the line number for the start of data
    data_start = [x+15 for x in hdrloc]
    data_end = [x for x in hdrloc]
    data_end.pop(0)
    data_end.append(nfilelines)
    
    # Check to see if Header ID string if found
    if nhdr == 0:
        raise ValueError('No headers found in file')
        return

    print("%d sonde profiles in file, working to split..." % nhdr)

    # Create a dictionary to hold the cls data
    cls = {'time': [],
    }

    # Create a time variable
    hdrpos = 0
    for ii in range(nhdr):
        # Create profile dictionary
        prof = {}
        f = open(filename, 'r')
        f.seek(hdrpos)
        
        print(ii)
        #index oout of bounds???
        hdrdict = _get_header(f)
        
        
        # Reset the header position for next loop
        
        # Save some data to the profile dictionary
        prof['altitude'] = hdrdict['alt']
        prof['release_latitude'] = hdrdict['lat0_dec']
        prof['release_longitude'] = hdrdict['lon0_dec']
        prof['project'] = hdrdict['project']
        prof['platform'] = hdrdict['site_type']
        prof['sonde_id'] = hdrdict['sonde_id']
        prof['site_id'] = hdrdict['site_id']
        
        n_elements = (data_end[ii] - data_start[ii]) * len(hdrdict['varnames'])
        data = np.fromfile(f, count=n_elements, sep=' ')
        data = np.reshape(data, ((data_end[ii]-data_start[ii]), len(hdrdict['varnames'])))
        
        hdrpos = f.tell()
        # Create a fields dictionary
        fields = {}
        for nv, var in enumerate(hdrdict['varnames']):
            fields[var] = {'data': data[:, nv],
                'units': hdrdict['varunits'][nv]}
    
        # MODIFY THE TIME FIELD AS IN THE REAST OF AWOT
        #       fields['Time'] =
        # Save the fields to the profile
        prof['fields'] = fields
        
        # Save the profile into the larger CLS instance
        cls[hdrdict['datetime0']] = prof
        #        cls['time'].append(hdrdict['datetime0'].strftime("%Y%m%d%H%M"))
        #        print(hdrdict['datetime0'].strftime("%Y%m%d%H%M"))
        
        
        # Save the profile into the larger CLS instance
        cls['time'].append(hdrdict['datetime0'])
        cls[hdrdict['datetime0'].strftime("%Y%m%d%H%M")] = prof



    return cls




