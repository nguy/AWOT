"""
awot.util.write_kmz
========================

Functions to save AWOT data into KML/KMZ files. These files may be displayed
for example with Google Earth.

Some code was directely adapted from the NASA PyAMPR package by Timothy Lang.
https://github.com/nasa/PyAMPR/blob/master/pyampr/pyampr.py

This present method is proof of concept and is expected to expand over time.
"""
from __future__ import absolute_import
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import os

from ..graph import common as gcommon
from . import helper
from .google_earth_tools import gearth_fig, make_kml
try:
    import simplekml
except:
    raise ValueError("This module requires installation of simplekml...")
#    return


def write_track_kmz(awot, field, lat_name=None, lon_name=None,
                    time_name=None, start_time=None, end_time=None,
                    latrange=None, lonrange=None,
                    cmap=None, track_color='k', track_lw=2.5,
                    file_path=None, file_name=None,
                    show_legend=True, legend_label=None):
    """
    This method plots geolocated AWOT track data as a filled color Google Earth
    kmz.
    Will produce overlay.png and, if a legend is created,
    legend.png as temporary image files in the current working
    directory.

    Parameters
    ----------
    awot : dict
        AWOT flight data instance.
    field : str
        Name of variable to use for track data.
    lat_name : str
        Key in radar instance for latitude variable.
        None uses AWOT default.
    lon_name : str
        Key in radar instance for longitude variable.
        None uses AWOT default.
    time_name : str
        Key in radar instance for time variable.
        None uses AWOT default.
    start_time : str
        UTC time to use as start time for subsetting in datetime format.
        (e.g. 2014-08-20 12:30:00)
    end_time : str
        UTC time to use as an end time for subsetting in datetime format.
        (e.g. 2014-08-20 16:30:00)
    latrange : 2-tuple
        List with lat range defined.
    lonrange : 2-tuple
        List with lon range defined.
    cmap : Matplotlib colormap instance
        Colormap desired.
        See http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
        and dir(cm) for more.
    track_color : str
        If color_by_altitude False, this color (see matplotlib)
        is used for track.
    track_lw : float or int
        Line width to use for track.
    file_path : str
        Path to kmz file. Defaults to current working directory.
    file_name : str
        Desired file name. If None specified, AWOT attemps
        to build using dictionary information.
    show_legend : bool
        False to suppress the color bar.
    legend_label : str
        Label to display in legend. If None, AWOT attempts to build
        using field dictionary information.
    """
    plt.close()  # mpl seems buggy if multiple windows are left open
    _method_printout()
    print('Writing AWOT track KMZ file:')

    # Check to see if field exists
    gcommon._check_field(awot, field)
    # Get the field dictionary
    var = gcommon._get_variable_dict(awot, field)
    # Get the lat/lon/time dictionaries
    if lat_name is None:
        latdict = gcommon._get_variable_dict(awot, 'latitude')
    else:
        latdict = gcommon._get_variable_dict(awot, lat_name)
    if lon_name is None:
        londict = gcommon._get_variable_dict(awot, 'longitude')
    else:
        londict = gcommon._get_variable_dict(awot, lon_name)
    if time_name is None:
        timedict = gcommon._get_variable_dict(awot, 'time')
    else:
        timedict = gcommon._get_variable_dict(awot, time_name)

    dt_start = gcommon._get_start_datetime(timedict, start_time)
    dt_end = gcommon._get_start_datetime(timedict, end_time)
    datasub = helper.time_subset_awot_dict(timedict, var,
                                           start_time, end_time)
    lonsub = helper.time_subset_awot_dict(timedict, londict,
                                          start_time, end_time)
    latsub = helper.time_subset_awot_dict(timedict, latdict,
                                          start_time, end_time)
    timesub = helper.time_subset_awot_dict(timedict, timedict,
                                           start_time, end_time)

    # Filter any bad geolocation data
    latd, lond, data, time = _filter_bad_geolocations(
        latsub['data'], lonsub['data'], timesub['data'], datasub['data'])

    # Get the lat/lon range
    latrange, lonrange = _get_latrange_lonrange(
        latd, lond, latrange, lonrange)

    # Set the beginning and ending times
    times = [dt_start, dt_end]

    # Set file info
    if file_path is None:
        file_path = os.getcwd()
    if file_name is None:
        file_name = ('awot_' + awot['platform'] + '_' + awot['flight_number'] +
                     '_' + field + '.kmz')
    longname = os.path.join(file_path, file_name)

    # Google Earth image production
    # If no cmap is specified, grab current
    if cmap is None:
        cmap = plt.get_cmap()

    print(lonrange)
    print(latrange)
    fig, ax = gearth_fig(np.min(lonrange), np.min(latrange),
                         np.max(lonrange), np.max(latrange))

    cs = ax.plot(lond, latd, color=track_color, lw=track_lw)
    ax.set_axis_off()
    fig.savefig('overlay.png', transparent=True, format='png')

    # Now we convert to KMZ
    if show_legend is True:
        fig = plt.figure(figsize=(1.0, 4.0), facecolor=None, frameon=False)
        ax = fig.add_axes([0.0, 0.05, 0.2, 0.9])
        cb = fig.colorbar(cs, cax=ax)
        cbytick_obj = plt.getp(cb.ax.axes, 'yticklabels')
        plt.setp(cbytick_obj, color='w', weight='bold')
        if legend_label is None:
            ptitle = var['standard_name'] + ' (' + var['units'] + ')'
        else:
            ptitle = legend_label
        cb.set_label(ptitle, rotation=-90, color='w', labelpad=20,
                     weight='bold')
        fig.savefig('legend.png', transparent=True, format='png')
        make_kml(np.min(lonrange), np.min(latrange), np.max(lonrange),
                 np.max(latrange), figs=['overlay.png'],
                 kmzfile=longname, colorbar='legend.png',
                 times=times)
        os.remove('overlay.png')
        os.remove('legend.png')
    else:
        make_kml(np.min(lonrange), np.min(latrange), np.max(lonrange),
                 np.max(latrange), figs=['overlay.png'],
                 kmzfile=longname, times=times)
        os.remove('overlay.png')

    print('Google Earth image saved to: %s' % longname)
    _method_printout()
    return


def write_line_kml(awot, field, lat_name=None, lon_name=None,
                    time_name=None, start_time=None, end_time=None,
                    latrange=None, lonrange=None,
                    cmap=None, color=None, lw=None,
                    file_path=None, file_name=None,
                    line_name=None, join_to_ground=False):
    """
    Write a field of AWOT track data as a KML file for Google Earth.

    Parameters
    ----------
    awot : dict
        AWOT flight data instance.
    field : str
        Name of variable to use for track data.
    lat_name : str
        Key in radar instance for latitude variable.
        None uses AWOT default.
    lon_name : str
        Key in radar instance for longitude variable.
        None uses AWOT default.
    time_name : str
        Key in radar instance for time variable.
        None uses AWOT default.
    start_time : str
        UTC time to use as start time for subsetting in datetime format.
        (e.g. 2014-08-20 12:30:00)
    end_time : str
        UTC time to use as an end time for subsetting in datetime format.
        (e.g. 2014-08-20 16:30:00)
    latrange : 2-tuple
        List with lat range defined.
    lonrange : 2-tuple
        List with lon range defined.
    cmap : Matplotlib colormap instance
        Colormap desired.
        See http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
        and dir(cm) for more.
    color : str
        Hex value color string to use for track.
    lw : float
        Line width to use for track.
    file_path : str
        Path to kmz file. Defaults to current working directory.
    file_name : str
        Desired file name. If None specified, AWOT attemps
        to build using dictionary information.
    line_name : str
        Name of line string to create.
    join_to_ground : bool
        True to create a line that is joined to the ground. False does not.
    """
    plt.close()  # mpl seems buggy if multiple windows are left open
    _method_printout()
    print('Writing AWOT track KML file:')

    # Check to see if field exists
    gcommon._check_field(awot, field)
    # Get the field dictionary
    var = gcommon._get_variable_dict(awot, field)
    # Get the lat/lon/time dictionaries
    if lat_name is None:
        latdict = gcommon._get_variable_dict(awot, 'latitude')
    else:
        latdict = gcommon._get_variable_dict(awot, lat_name)
    if lon_name is None:
        londict = gcommon._get_variable_dict(awot, 'longitude')
    else:
        londict = gcommon._get_variable_dict(awot, lon_name)
    if time_name is None:
        timedict = gcommon._get_variable_dict(awot, 'time')
    else:
        timedict = gcommon._get_variable_dict(awot, time_name)

    dt_start = gcommon._get_start_datetime(timedict, start_time)
    dt_end = gcommon._get_start_datetime(timedict, end_time)
    datasub = helper.time_subset_awot_dict(timedict, var,
                                           start_time, end_time)
    lonsub = helper.time_subset_awot_dict(timedict, londict,
                                          start_time, end_time)
    latsub = helper.time_subset_awot_dict(timedict, latdict,
                                          start_time, end_time)
    timesub = helper.time_subset_awot_dict(timedict, timedict,
                                           start_time, end_time)

    # Filter any bad geolocation data
    latd, lond, data, time = _filter_bad_geolocations(
        latsub['data'], lonsub['data'], timesub['data'], datasub['data'])

    # Get the lat/lon range
    latrange, lonrange = _get_latrange_lonrange(
        latd, lond, latrange, lonrange)

    # Set the beginning and ending times
    times = [dt_start, dt_end]

    # Set file info
    if file_path is None:
        file_path = os.getcwd()
    if file_name is None:
        file_name = ('awot_' + awot['platform'] + '_' + awot['flight_number'] +
                     '_' + field + '.kml')
    longname = os.path.join(file_path, file_name)

    if line_name is None:
        line_name = field

    # Now we convert to KMZ
    kml = simplekml.Kml()
    linestr = kml.newlinestring(name=line_name)
    linestr.coords = zip(lond, latd)

    # Set properties according to keywords
    if lw is not None:
        linestr.style.linestyle.width = lw
    if color is not None:
        linestr.style.linestyle.color = color
    if join_to_ground:
        linestr.altitudemode = simplekml.AltitudeMode.relativetoground
        linestr.extrude = 1

    # Save the file
    kml.save(longname)
    print('KML file saved to: %s' % longname)
    _method_printout()
    return


def write_poly_kml(name="Polygon", innerboundary=None, outerboundary=None,
                   lw=4, color=None, fill=False, extrude=False,
                   join_to_ground=False,
                   timestamp=False, outputdir=None, filename=None):
    '''
    Write a polygon KML file for Google Earth.

    Parameters
    ----------
    name : str
        Name to assign KML polygon instance.
    innerboundary : tuple or array
        Longitude/latitude coordinate pairs of inner boundary.
    outerboundary : tuple or array
        Longitude/latitude coordinate pairs of outer boundary.
    lw : int
        Width of line connecting polygon points.
    color : str
        Hex string color of line connection polygon points.
    fill : bool
        True to fill in polygon, False (default) for no fill.
    extrude : bool
        True to connect line to ground, Fales (default) does not.
    join_to_ground : bool
        True to create a line that is joined to the ground. False does not.
    timestamp : bool
        True adds UTC time at processing.
    outputdir : str
        Directory to save KML file. Defaults to current working directory.
    filename : str
        KML filename to save. Defualts to polygon.kml
    '''
    # Create a kml instance
    kml = simplekml.Kml()

    # Instantiate a polygon instance
    pol = kml.newpolygon(name=name)

    # Set the boundaries
    if outerboundary is not None:
        pol.outerboundaryis = outerboundary
    if innerboundary is not None:
        pol.innerboundaryis = innerboundary

    if fill:
        pol.polystyle.fill = 1
    else:
        pol.polystyle.fill = 0
    if extrude:
        pol.extrude = 1
    else:
        pol.extrude = 0
    if join_to_ground:
        pol.altitudemode = simplekml.AltitudeMode.relativetoground
        pol.extrude = 1


    # Set styling properties according to keywords
    if lw is not None:
        pol.linestyle.width = lw
    if color is not None:
        pol.linestyle.color = simplekml.Color.hex(color)

    if timestamp:
        pol.timestamp.when = datetime.strftime(datetime.utcnow(),
                                               '%Y-%m-%d %H:%M:%SZ')
    # Save the file
    if outputdir is None:
        outputdir = os.getcwd()
    if filename is None:
        filename = "polygon"
    outfname = os.path.join(outputdir, filename + ".kml")
    kml.save(outfname)
    return


def write_kmz(fig, ax, plot, lonrange, latrange, times,
              file_path=None, file_name=None,
              show_legend=True, legend_label=None):
    """
    This method plots geolocated AWOT track data as a filled color Google Earth
    kmz.
    Will produce overlay.png and, if a legend is created,
    legend.png as temporary image files in the current working
    directory.

    Parameters
    ----------
    latrange : 2-tuple
        List with lat range defined.
    lonrange : 2-tuple
        List with lon range defined.
    file_path : str
        Path to kmz file. Defaults to current working directory.
    file_name : str
        Desired file name. If None specified, AWOT attemps
        to build using dictionary information.
    show_legend : bool
        False to suppress the color bar.
    legend_label : str
        Label to display in legend. If None, AWOT attempts to build
        using field dictionary information.
    """
    plt.close()  # mpl seems buggy if multiple windows are left open
    _method_printout()
    print('Writing AWOT track KMZ file:')

    longname = os.path.join(file_path, file_name)

    print(lonrange)
    print(latrange)
    ax.set_axis_off()
    fig.savefig('overlay.png', transparent=True, format='png')

    # Now we convert to KMZ
    if show_legend is True:
        fig = plt.figure(figsize=(1.0, 4.0), facecolor=None, frameon=False)
        ax = fig.add_axes([0.0, 0.05, 0.2, 0.9])
        cb = fig.colorbar(plot, cax=ax)
        cbytick_obj = plt.getp(cb.ax.axes, 'yticklabels')
        plt.setp(cbytick_obj, color='w', weight='bold')
        cb.set_label(legend_label, rotation=-90, color='w', labelpad=20,
                     weight='bold')
        fig.savefig('legend.png', transparent=True, format='png')
        make_kml(np.min(lonrange), np.min(latrange), np.max(lonrange),
                 np.max(latrange), figs=['overlay.png'],
                 kmzfile=longname, colorbar='legend.png',
                 times=times)
        os.remove('overlay.png')
        os.remove('legend.png')
    else:
        make_kml(np.min(lonrange), np.min(latrange), np.max(lonrange),
                 np.max(latrange), figs=['overlay.png'],
                 kmzfile=longname, times=times)
        os.remove('overlay.png')

    print('Google Earth image saved to: %s' % longname)
    _method_printout()
    return

###################
#   Get methods   #
###################


def _get_latrange_lonrange(lats=None, lons=None, latrange=None, lonrange=None):
    """ Determine domain of plot based on what user provided. """
    if latrange is None:
        latrange = [np.min(lats), np.max(lats)]
    else:
        latrange = np.sort(latrange)
    if lonrange is None:
        lonrange = [np.min(lons), np.max(lons)]
    else:
        lonrange = np.sort(lonrange)
    return latrange, lonrange

####################
#   Data methods   #
####################


def _filter_bad_geolocations(lats, lons, data, time):
    """ Internal method to filter bad geolocation data. """
    # Attempt to deal with bad geolocation data
    # (e.g., Latitude/Longitude=bad_data)
    cond1 = np.logical_or(lats < -90, lats > 90)
    cond2 = np.logical_or(lons < -180, lons > 180)
    condition = np.logical_or(cond1, cond2)
    indices = np.where(condition)
    if np.shape(indices)[1] > 0:
        print("Removing bad data")
        data = np.delete(data, indices[0], axis=0)
        lons = np.delete(lons, indices[0], axis=0)
        lats = np.delete(lats, indices[0], axis=0)
        time = np.delete(time, indices[0], axis=0)
    return lats, lons, data, time

######################
#   Helper methods   #
######################


def _method_printout():
    """ Helps clarify text output. """
    print('\n********************\n')
    print('')
