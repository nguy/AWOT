from matplotlib.ticker import ScalarFormatter, MultipleLocator
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
from .thermocalcs import ThermoCalcs
from .shearcalcs import ShearCalcs
from .skew import SkewXTick


# Need to change:
# ndarrays to dictionaries
# split file to work with current data gathering 
#
#
#




def get_sounding_data(filePath):

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
    data['Header'] = header
    data['temperature'] = T
    data['dewpoint'] = Td
    data['presssure'] = p
    data['relative_humidity'] = RH
    data['u_component'] = Uwind
    data['v_component'] = Vwind
    data['Height'] = h
    data['Type'] = 'radioSonde'

    fp.close()

    # type = 'radioSonde'

    return data

def get_dropsonde_data(filePath):

    '''
    method to retrieve dropsonde data

    Inputs:

    file path
    #============

    output:

    #Dictionary with variable for use in plotting

    '''

    data = filePath
    fp = open(data, 'r')
    lines = fp.readlines()
    header = lines[1]
    fp.close()
    fp = open(data, 'r')
    sounding_data = fp

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
    data['Header'] = header
    data['temperature'] = T
    data['dewpoint'] = TD
    data['presssure'] = P
    data['relative_humidity'] = RH
    data['u_component'] = Uwind
    data['v_component'] = Vwind
    data['Height'] = Height
    data['Type'] = 'dropsonde'

    return data

def plot_skewt_logp(data, **kwargs):

    '''
    method to plot sounding or dropsonde data

    Inputs:

    dictionary of sounding data

    keyword arguments min, max, titles, and label.

    #============

    output:

    #Image

    '''

    T = data['temperature']
    TD = data['dewpoint']
    P = data['presssure']

    fig = plt.figure(figsize=(10, 8))
    ax1 = fig.add_axes([0.05, 0.1, 0.6, 0.8],
                                 projection='skewx')
    plt.grid(True)

    plot_dryadiabats()

    ax1.semilogy(T, P, 'r-', linewidth=1.5)
    ax1.semilogy(TD, P, 'g-', linewidth=1.5)

    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dictated by the typical meteorological plot
    # Disables the log-formatting that comes with semilogy
    # set bounds for data on Skew T chart

    ax1.yaxis.set_major_formatter(ScalarFormatter())
    ax1.set_yticks(np.linspace(100, 1000, 10))
    ax1.set_ylim(1050, 100)
    ax1.set_ylabel('presssure mb')
    ax1.set_xlabel('temperature C')
    ax1.xaxis.set_major_locator(MultipleLocator(10))
    ax1.set_title('Skew T log P')

    y_min = kwargs.get('y_min')
    y_max = kwargs.get('y_max')

    x_min = kwargs.get('x_min')
    x_max = kwargs.get('x_max')

    # data label, title and font sizes value assignments
    # default values above.

    y_label = kwargs.get('y_label')
    x_label = kwargs.get('x_label')
    title = kwargs.get('title')
    font_size = kwargs.get('font_size')

    for item in ([ax1.xaxis.label, ax1.yaxis.label] +
                 ax1.get_xticklabels()+ax1.get_yticklabels()):
        item.set_fontsize(font_size)

    # set the labels and titles defaults to 'none'

    ax1.set_ylabel(y_label)
    ax1.set_xlabel(x_label)
    ax1.set_title(title)

    ax1.set_ylim(y_max, y_min)
    ax1.set_xlim(x_min, x_max)

def plot_hodograph(data):

    '''
    method to plot wind data on a hodograph

    Inputs:

    dictionary containing wind data.

    file path
    # ============

    output:

    # Image

    '''

    ax2 = fig.add_axes([.05, 0.6, 0.25, 0.3])

    # create axis and invert masks and assign values
    # for U, V, and h coordinates

    Uwind = data['u_component']
    Vwind = data['v_component']
    Height = data['Height']
    
    
    u_3km = []
    v_3km = []
    u_6km = []
    v_6km = []
    u_9km = []
    v_9km = []
    u_12km = []
    v_12km = []

    # =====================================================#
    # loop to place hodograph data into proper height array
    # 0-3km
    # 3-6km
    # 6-9km
    # 9-12km
    # ======#

    if data['Type'] == 'radioSonde':

        for unew, vnew, hnew in zip(U, V, H):
            if hnew <= 3001.0:
                u_3km.append(unew)
                v_3km.append(vnew)
            elif hnew <= 6001.0:
                u_6km.append(unew)
                v_6km.append(vnew)
            elif hnew <= 9000.0:
                u_9km.append(unew)
                v_9km.append(vnew)
            elif hnew <= 12000.0:
                u_12km.append(unew)
                v_12km.append(vnew)

    # print(len(u_3km))
    # print(len(v_3km))
    # =============================================================#
    # take the first point of the next height level and attach it
    # to the last point of the previous line segment to
    # connect hodograph segments
    # should only be done for soundings and non dense data
    # type determines how the data is placed in the array.
    # Soundings require endpoints be joined so that
    # the hodograph is continuous.
    # =============================================================#

    if data['Type'] == 'radioSonde':

        u_3km.append(u_6km[0])
        u_6km.append(u_9km[0])
        u_9km.append(u_12km[0])

        v_3km.append(v_6km[0])
        v_6km.append(v_9km[0])
        v_9km.append(v_12km[0])

    # ==============================#
    # Mask the arrays U and V coords
    # ==============================#

    v_3km = np.ma.asarray(v_3km)
    v_6km = np.ma.asarray(v_6km)
    v_9km = np.ma.asarray(v_9km)
    v_12km = np.ma.asarray(v_12km)

    u_3km = np.ma.asarray(u_3km)
    u_6km = np.ma.asarray(u_6km)
    u_9km = np.ma.asarray(u_9km)
    u_12km = np.ma.asarray(u_12km)

    # plotting of the hodograph information
    # Different colors are sued to represent different height levels
    # r = 0-3km
    # y = 3-6km
    # g = 6-9km
    # b = 9-12km

    ax2.plot(u_3km, v_3km, 'r-', linewidth=3)
    ax2.plot(u_6km, v_6km, 'y-', linewidth=3,)
    ax2.plot(u_9km, v_9km, 'g-', linewidth=3)
    ax2.plot(u_12km, v_12km, 'b-', linewidth=3)

    # set the default limits on the hodo axes
    # remove spines and set the axes to be invisible.

    ax2.spines['left'].set_position('zero')
    ax2.spines['right'].set_color('none')
    ax2.spines['bottom'].set_position('zero')
    ax2.spines['top'].set_color('none')
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    ax2.set_ylim(-50, 50)
    ax2.set_xlim(-50, 50)

    # Draw hodograph range rings
    # 10 20 30 40 m/s circles

    circ10 = patches.Circle((0, 0), 10, fc='white')
    circ20 = patches.Circle((0, 0), 20, fc='white')
    circ30 = patches.Circle((0, 0), 30, fc='white')
    circ40 = patches.Circle((0, 0), 40, fc='white')

    # Add circles to plot as patches#
    # descending order to make sure they layer ontop of each other.

    ax2.add_patch(circ40)
    ax2.add_patch(circ30)
    ax2.add_patch(circ20)
    ax2.add_patch(circ10)

def plot_aux_graph(x_value, y_value, **kwargs):

    '''
    method to plot an auxiliary graph of user defined data

    Inputs:

    # sounding data (User specified)
    # kwargs to adjust plot dimensions scales label, and title.
    # ============

    output:

    # Image

    '''

    # define axes fro figure.
    # set grid and plot user specified data

    ax3 = fig.add_axes([.7, 0.7, .29, .24])
    plt.grid(True)
    ax3.plot(x_value, y_value)

    # Max and min graph bonds value assignment. Defaults to autoscale.

    y_min = kwargs.get('y_min')
    y_max = kwargs.get('y_max')

    x_min = kwargs.get('x_min')
    x_max = kwargs.get('x_max')

    # data label, title and font size value assignments .

    y_lable = kwargs.get('y_lable')
    x_lable = kwargs.get('x_lable')
    title = kwargs.get('title')
    font_size = kwargs.get('font_size')

    # keyword arguments for minimum and maximum values

    ax3.set_ylim(y_min, y_max)
    ax3.set_xlim(x_min, x_max)

    # set label and tick fontsizes to user defined value. Defaults to 10pt.

    for item in ([ax3.xaxis.label, ax3.yaxis.label] +
                 ax3.get_xticklabels()+ax3.get_yticklabels()):
        item.set_fontsize(font_size)

    # set the labels and title from keyword arguments. Defaults to 'none'

    ax3.set_ylabel(y_lable)
    ax3.set_xlabel(x_lable)
    ax3.set_title(title)

def generate_parameter_list():

    '''
    method to generate a list of parameters calculated
    from the sounding data

    Inputs:

    # None
    # ============

    output:

    # Blank axis for plotting
    '''

    # Set axes position and set both axes invisible

    ax4 = fig.add_axes([.65, 0.1, .35, .54])
    ax4.xaxis.set_visible(False)
    ax4.yaxis.set_visible(False)

def run_thermo_calcs(data):

    '''
    method to calculate thermodynamic parameters from dropsonde data

    Inputs:

    # Dictionary of sounding data.
    # Uses calculations from thermocalcs.py
    # ============

    output:

    # multiple arrays containing thermodynamic information

    '''
    
    tC = ThermoCalcs()

    T = data['temperature']
    Td = data['dewpoint']
    p = data['presssure']
    RH = data['relative_humidity']
    u = data['u_component']
    v = data['v_component']
    h = data['Height']

    LCLT = round((
        tC._LCL_temperature(h, T+273.15, Td+273.15)-273.15), 2)
    LCLP = round((
        tC._LCL_presssure(h, p, T+273.15, Td+273.15)), 0)
    LCLZ = round((
        tC._LCL_Height(h, p, T+273.15, Td+273.15)), 0)
    THETA = tC._PTk2_Theta(p, T+273.15)
    MIXR = tC._RH_2_MixR(RH, p, T+273.15)
    THETAE = tC._Tk_RH_MixR_2_ThetaE(
        p, T+273.15, RH, MIXR/1000.)
    ESAT = tC._esat(T+273.15)
    
    thermoCalcData = dict()
    thermoCalcData['lclt'] = LCLT
    thermoCalcData['lclp'] = LCLP
    thermoCalcData['lclz'] = LCLZ
    thermoCalcData['theta'] = THETA
    thermoCalcData['mixr'] = MIXR
    thermoCalcData['thetae'] = THETAE
    thermoCalcData['esat'] = ESAT


    return thermoCalcData
    


def plot_thermo_calcs():

    '''
    method to plot the thermodynamic parameters on the parameter list.

    Inputs:

    # None
    # ============

    output:

    # Image
    '''

    # plot the parameters on the list generated.
    ax4.text(.01, .01, 'LCL presssure: ' + str(lclp) + (' hPa'))
    ax4.text(.01, .04, 'LCL Temp: ' + str(lclt) + ' c')
    ax4.text(.01, .07, 'LCL Height: ' + str(lclz) + ' m')

def plot_dryadiabats(**kwargs):

    '''
    method to plot the dry adibats. Used in the plotskewtlogp method.

    Inputs:

    # kwargs (does not function)
    # ============

    output:

    # Image

    '''
    # test = shear1km

    # temperature array and presssure array

    t0 = np.linspace(200, 430, 17)
    press = np.linspace(100, 1000.)

    # retrieve desired line style from user kwargs

    line_style = kwargs.get('line_style')
    if line_style is None:
        line_style = '-'

    # loop to calculate using Poisson's equation the
    # dry adiabats for the skew t.

    for temp in t0:
        theta = temp * (press/1000.)**(2./7.)

        # plot the dry adiabats given a specified color

        ax1.semilogy((theta-273.15), press, line_style,
                          color='#7F4B10', linewidth=0.5)

def plot_wind_barbs(data, **kwargs):

    P = data['presssure']
    Uwind = data['u_component']
    Vwind = data['v_component']
    Height = data['Height']

    # Copy y axis to plot wind barbs
    # Set x lim for windpbarbs
    # adjust location of barbs from x =0
    # plot every 30th windbarb
    # set axis to invisible
    ax1_copy = ax1.twiny()
    ax1_copy.set_xlim(x_min, x_max)
    ax1_copy.set_ylim(y_min, y_max)
    x_const = np.zeros(P.shape) + (x_max - 2)
    ax1_copy.xaxis.set_visible(False)
    
    
    #choses specific wind barb plotting pattern based on the density of wind observations.
    #Need to make this a function of the number of observations per second/vertical velocity 

    if data['Type'] == 'radioSonde':

        mask = U.mask
        Uwind = Uwind[~mask]
        Vwind = Vwind[~mask]
        P = P[~mask]

        ax1_copy.barbs(x_const[::3], P[::3], Uwind[::3], Vwind[::3])

    else:

        ax1_copy.barbs(x_const[::40], P[::40], Uwind[::40], Vwind[::40])

def run_shear_calcs(data):

    '''
    method to calculate thermodynamic parameters from dropsonde data

    Inputs:

    # Dictionary of sounding data.
    # Uses calculations from thermocalcs.py
    # ============

    output:

    # multiple arrays containing thermodynamic information

    '''

    T = data['temperature']
    Td = data['dewpoint']
    p = data['presssure']
    RH = data['relative_humidity']
    uwind = data['u_component']
    vwind = data['v_component']
    height = data['Height']

    mask = h.mask
    uwind = u[~mask]
    vwind = v[~mask]
    height = h[~mask]

    SHEAR1KM = sC._VertShear_Sfc_to_1km(h, u, v)
    SHEAR3KM = sC._VertShear_Sfc_to_3km(h, u, v)
    SHEAR6KM = sC._VertShear_Sfc_to_6km(h, u, v)
    BULKSHEAR1km = round(sC._bulkshear_sfc_1km(h, u, v), 2)
    BULKSHEAR3km = round(sC._bulkshear_sfc_3km(h, u, v), 2)
    BULKSHEAR6km = round(sC._bulkshear_sfc_6km(h, u, v), 2)


    shearCalcData = dict()
    shearCalcData['SHEAR1KM'] = SHEAR1KM
    shearCalcData['SHEAR3KM'] = SHEAR3KM
    shearCalcData['SHEAR6KM'] = SHEAR6KM
    shearCalcData['BULKSHEAR1km'] = BULKSHEAR1km
    shearCalcData['BULKSHEAR1km'] = BULKSHEAR3km
    shearCalcData['BULKSHEAR1km'] = BULKSHEAR6km

    return shearCalcData



def plot_shear_calcs():

    '''
    method to plot the thermodynamic parameters on the parameter list.

    Inputs:

    # None
    # ============

    output:

    # Image
    '''



    # plot the parameters on the list generated.
    ax4.text(.01, .1, '0-1 km shear: '+str(SHEAR1KM)+(' 1/s'))
    ax4.text(.01, .14, '0-3 km shear: '+str(SHEAR3KM)+' 1/s')
    ax4.text(.01, .18, '0-6 km shear: '+str(SHEAR6KM)+' 1/s')
    ax4.text(
        .01, .22, '0-1km Bulk Shear: '+str(BULKSHEAR1km)+' m/s')
    ax4.text(
        .01, .26, '0-3km Bulk Shear: '+str(BULKSHEAR3km)+' m/s')
    ax4.text(
        .01, .3, '0-6km Bulk Shear: '+str(BULKSHEAR6km)+' m/s')

#def dry_lift(data):
#
#    T = data['temperature']
#    Td = data['dewpoint']
#    p = data['presssure']
#    RH = data['relative_humidity']
#    u = data['u_component']
#    v = data['v_component']
#    h = data['Height']
#
#    t_parcel, p_parcel = tC.dry_lift(T, p, LCLT, LCLP)
#
#    ax1.semilogy(t_parcel, p_parcel, 'k--', ms=1)


def splitFile(filename):
    
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
