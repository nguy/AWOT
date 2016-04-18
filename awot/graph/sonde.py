from matplotlib.ticker import ScalarFormatter, MultipleLocator
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

from . import common
from .skew import SkewXTick

def plot_skewt_logp(data, instance,
                    temp_key=None, temp_color='r', temp_ls='-', temp_lw=1.5,
                    dewpt_key=None, dewpt_color='g', dewpt_ls='-', dewpt_lw=1.5,
                    press_key=None,
                    x_min=None, x_max=None, p_min=None, p_max=None,
                    title=None, titleFontSize=None,
                    xlab=None, xlabFontSize=None, xpad=None,
                    ylab=None, ylabFontSize=None, ypad=None,
                    ax=None, fig=None, **kwargs):

    '''
    Plot sounding or dropsonde data on Skew-T diagram.

    Parameters
    ----------
    data: dict
        Data dictionary of dictionaries of sounding data variables.
    instance: str
        A date time string used to pick the dropsonde to be plotted.
    temp_key: str
        Optional keyword for the temperature field.
        Default is 'Temp'.
    dewpt_key: str
        Optional keyword for the dewpoint temperature field.
        Default is 'Dewpt'.
    press_key: str
        Optional keyword for the pressure field.
        Default is 'Press'.
    x_min : float
        Minimum value for x-axis.
    x_max : float
        Maximum value for x-axis.
    p_min : float
        Minimum value for pressure-axis.
    p_max : float
        Maximum value for pressure-axis.
    title : str
        Plot title.
    titleFontSize : int
        Font size to use for Title label.
    xlab : str
        X-axis label.
    ylab : str
        Y-axis label.
    xpad : int
        Padding for X-axis label.
    ypad : int
        Padding for Y-axis label.
    xlabFontSize : int
        Font size to use for X-axis label.
    ylabFontSize : int
        Font size to use for Y-axis label.
    ax : Matplotlib axis instance
        Axis to plot. None will use the current axis.
    fig : Matplotlib figure instance
        Figure on which to add the plot. None will use the current figure.
    '''
    # pulling from new data structure
    if temp_key is None:
        try:
            T = data[instance]['fields']['Temp']['data'][:]
        except:
            T = data[instance]['fields'][temp_key]['data'][:]
            raise ValueError('Please check the temperature key!')
            return None
        try:
            TD = data[instance]['fields']['Dewpt']['data']
        except:
            TD = data[instance]['fields'][dewpt_key]['data']
            raise ValueError('Please check the dewpoint temperature key!')
            return None
        try:
            P = data[instance]['fields']['Press']['data'][:]
        except:
            P = data[instance]['fields'][press_key]['data'][:]
            raise ValueError('Please check the pressure key!')
            return None
    sub1 = (~T.mask & ~P.mask)
    sub2 = (~TD.mask & ~P.mask)
    T = T[sub1]
    P1 = P[sub1]
    TD = TD[sub2]
    P2 = P[sub2]
    if fig is None:
        fig = plt.figure(figsize=(10, 8))
    if ax is None:
        ax_skew = fig.add_axes([0.05, 0.1, 0.6, 0.8], projection='skewx')
    else:
        ax_skew = ax
    # Set up the axes
    if p_min is None:
        p_min = 100.
    if p_max is None:
        p_max = 1050.
    if xlab is None:
        xlab = 'Temperature (C)'
    if ylab is None:
        ylab = 'Pressure (hPa)'
    if title is None:
        title = 'Skew T - log P'
    common._set_axes(ax_skew, x_min=x_min, x_max=x_max,
                     y_min=p_max, y_max=p_min,
                     title=title, titleFontSize=titleFontSize,
                     xlab=xlab, ylab=ylab, xpad=xpad, ypad=ypad,
                     xlabFontSize=xlabFontSize, ylabFontSize=ylabFontSize)
    plt.grid(True)

    plot_dryadiabats(ax_skew)
    l = ax_skew.axvline(0, color='k', ls='-', lw=2)

    Tp = ax_skew.semilogy(T, P1, color=temp_color, ls=temp_ls, lw=1.5)
    TDp = ax_skew.semilogy(TD, P2, color=dewpt_color, ls=dewpt_ls, lw=1.5)

    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dictated by the typical meteorological plot
    # Disables the log-formatting that comes with semilogy
    # set bounds for data on Skew T chart

    ax_skew.yaxis.set_major_formatter(ScalarFormatter())
#    ax_skew.set_yticks(np.linspace(100, 1000, 10))
    ax_skew.set_yticks(np.linspace(p_min, 1000, 10))
    ax_skew.xaxis.set_major_locator(MultipleLocator(10))

    for item in ([ax_skew.xaxis.label, ax_skew.yaxis.label] +
                 ax_skew.get_xticklabels() + ax_skew.get_yticklabels()):
        item.set_fontsize(xlabFontSize)
    return


def plot_hodograph(data, instance):

    '''
    method to plot wind data on a hodograph

    Parameters
    ----------
    data: dict
        Data dictionary of dictionaries of sounding data variables.
    instance: string
        A date time string used to pick the dropsonde to be plotted.

    '''
    ax2 = fig.add_axes([.05, 0.6, 0.25, 0.3])
    # create axis and invert masks and assign values
    # for U, V, and h coordinates

    Uwind = data[instance]['fields']['Ucmp']['data']
    Vwind = data[instance]['fields']['Vcmp']['data']
    Height = data[instance]['fields']['Alt']['data']

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

    Parameters
    ----------
    xvalue : arraylike
        x value of auxillary sounding data to be plotted.
    yvalue : arraylike
        Y value of auxillary sounding data to be plotted.

    **kwargs
        Adjust plot dimensions scales label, and title.
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


def plot_parameter_list():

    '''
    Method to generate a list of parameters calculated
    from the sounding data.
    Blank axis for plotting.

    Parameters
    ---------

    '''
    # Set axes position and set both axes invisible

    ax4 = fig.add_axes([.65, 0.1, .35, .54])
    ax4.xaxis.set_visible(False)
    ax4.yaxis.set_visible(False)


def plot_thermo_calcs():

    return
#
#    '''
#    method to plot the thermodynamic parameters on the parameter list.
#
#    Inputs
#    '''
#
#    # plot the parameters on the list generated.
#    ax4.text(.01, .01, 'LCL presssure: ' + str(lclp) + (' hPa'))
#    ax4.text(.01, .04, 'LCL Temp: ' + str(lclt) + ' c')
#    ax4.text(.01, .07, 'LCL Height: ' + str(lclz) + ' m')


def plot_dryadiabats(ax, **kwargs):

    '''
    method to plot the dry adibats. Used in the plotskewtlogp method.

    Parameters
    ----------
    kwargs (does not function).

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

        ax.semilogy((theta - 273.15), press, line_style,
                    color='#7F4B10', linewidth=0.5)


def plot_wind_barbs(data, instance, **kwargs):
    '''
    Method to plot wind barbs on the skewT

    Parameters
    ----------
    data: dict
        data dictionary of dictionaries of sounding data variables.
    instance: string
        a date time string used to pick the dropsonde to be plotted.
    kwargs: string max,min.
    '''
    P = data[instance]['fields']['Press']['data']
    Uwind = data[instance]['fields']['Ucmp']['data']
    Vwind = data[instance]['fields']['Vcmp']['data']
    Height = data[instance]['fields']['Alt']['data']

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

    if data['Type'] == 'radioSonde':

        mask = U.mask
        Uwind = Uwind[~mask]
        Vwind = Vwind[~mask]
        P = P[~mask]

        ax1_copy.barbs(x_const[::3], P[::3], Uwind[::3], Vwind[::3])

    else:

        ax1_copy.barbs(x_const[::40], P[::40], Uwind[::40], Vwind[::40])


def plot_shear_calcs():
    #
    #    '''
    #    method to plot the thermodynamic parameters on the parameter list.
    #
    #    Input
    #
    #
    #
    #    # plot the parameters on the list generated.
    #    ax4.text(.01, .1, '0-1 km shear: '+str(SHEAR1KM)+(' 1/s'))
    #    ax4.text(.01, .14, '0-3 km shear: '+str(SHEAR3KM)+' 1/s')
    #    ax4.text(.01, .18, '0-6 km shear: '+str(SHEAR6KM)+' 1/s')
    #    ax4.text(
    #        .01, .22, '0-1km Bulk Shear: '+str(BULKSHEAR1km)+' m/s')
    #    ax4.text(
    #        .01, .26, '0-3km Bulk Shear: '+str(BULKSHEAR3km)+' m/s')
    #    ax4.text(
    #        .01, .3, '0-6km Bulk Shear: '+str(BULKSHEAR6km)+' m/s')

    # def dry_lift(data):
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
        return
