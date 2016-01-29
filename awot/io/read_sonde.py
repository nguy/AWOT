import numpy as np
import matplotlib.patches as patches
from .common import _build_dict



# Need to change:
# ndarrays to dictionaries
# split file to work with current data gathering
#
#
#




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
    data['temperature'] = _build_dict(T, 'c', 'Temperature of ambient air', 'Temperature')
    data['dewpoint'] = _build_dict(TD, 'c', 'Dewpoint temperature of ambient air', 'Dewpoint_Temperature')
    data['presssure'] = _build_dict(P, 'hPa', 'Pressure of ambient air', 'Pressure')
    data['relative_humidity'] = _build_dict(RH, '%', 'Relative Humidity of ambient air’, ‘Relative_Humidity')
    data['u_component'] = _build_dict(Uwind, 'm/s', 'u component of wind', 'U_component')
    data['v_component'] = _build_dict(Vwind, 'm/s', 'v component of wind', 'V_component')
    data['height'] = _build_dict(H, 'm', 'Geometric Height in meters', 'Height')
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
    data['temperature'] = _build_dict(T, 'c', 'Temperature of ambient air', 'Temperature')
    data['dewpoint'] = _build_dict(TD, 'c', 'Dewpoint temperature of ambient air', 'Dewpoint Temperature')
    data['presssure'] = _build_dict(P, 'hPa', 'Pressure of ambient air', 'Pressure')
    data['relative_humidity'] = _build_dict(RH, '%', 'Relative Humidity of ambient air’, ‘Relative Humidity')
    data['u_component'] = _build_dict(Uwind, 'm/s', 'u component of wind', 'U component')
    data['v_component'] = _build_dict(Vwind, 'm/s', 'v component of wind', 'V component')
    data['height'] = _build_dict(H, 'm', 'Geometric Height in meters', 'Height')
    data['data_format'] = 'dropsonde'

    return data


#splits a cls file for use in the programm(needs fixing)



def split_cls_file(filename):
    
    splitString1 = 'Data Type:'
    splitString2 = 'Project ID:'
    splitString3 = '------'
    
    f = open(filename,'r')
    newString = ''
    splitFileList = []
    f.readline()
    s =0
    x=0
    for line in f:
        if  (splitString1 or splitString2) not in line:
            newString = newString + line
        
        else:
            splitFileList.append(newString)
            newString = ''

    splitFileList.append(newString)
    
    return splitFileList




