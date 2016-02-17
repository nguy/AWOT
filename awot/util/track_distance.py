"""
awot.util.track_distance
========================

Convert existing data to track distances.


Copyright (c) 2006-2015 geopy authors (see AUTHORS)

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
"""
from __future__ import division

from netCDF4 import date2num
import numpy as np

from ..graph.common import _get_start_datetime, _get_end_datetime
from ..io.common import _build_dict

# From http://www.movable-type.co.uk/scripts/LatLongVincenty.html:
#   The most accurate and widely used globally-applicable model for the earth
#   ellipsoid is WGS-84, used in this script. Other ellipsoids offering a
#   better fit to the local geoid include Airy (1830) in the UK, International
#   1924 in much of Europe, Clarke (1880) in Africa, and GRS-67 in South
#   America. America (NAD83) and Australia (GDA) use GRS-80, functionally
#   equivalent to the WGS-84 ellipsoid.
ELLIPSOIDS = {
    # model           major (km)   minor (km)     flattening
    'WGS-84':        (6378.137, 6356.7523142, 1 / 298.257223563),
    'GRS-80':        (6378.137, 6356.7523141, 1 / 298.257222101),
    'Airy (1830)':   (6377.563396, 6356.256909, 1 / 299.3249646),
    'Intl 1924':     (6378.388, 6356.911946, 1 / 297.0),
    'Clarke (1880)': (6378.249145, 6356.51486955, 1 / 293.465),
    'GRS-67':        (6378.1600, 6356.774719, 1 / 298.25)
}


def calc_ground_distance(awot, method='great circle',
                         radius=None, ellipsoid_key=None, iterations=None,
                         lon_key=None, lat_key=None,
                         add_to_dict=False, keyname=None, units=None,
                         longname=None, stdname=None):
    '''
    Calculate the distance travelled in reference to the ground.

    The calculation assumes decimal latitude and longitude coordinates.
    The calculation is adapted from the geopy package, see copyright at the
    top of this module.
    https://github.com/geopy/geopy

    Parameters
    ----------
    awot : dict
        AWOT data dictionary instance.
    method : str
        Either 'great circle' or 'vincentry' are accepted.
    radius : float
        Radius of earth. Default is to use the average great-circle radius
        of a sphere in meters. According to geopy documentation, using a
        sphere with this radius results in an error of up to about 0.5%.
        This keyword is passed to great_circle calculation.
    ellipsoid_key : str
        Global ellipsoid model to use in calculations. 'WGS-84' is default.
        This keyword is passed to vincentry calculation.
    iterations : int
        Number ov iterations to try for convergent solution.
        This keyword is passed to vincentry calculation.
    lon_key : str
        Key to use for longitude variable in AWOT dictionary.
    lat_key : str
        Key to use for latitude variable in AWOT dictionary.
    add_to_dict : bool
        True to add the calculated distance to the AWOT dictionary.
    keyname : str
        Key to be used when adding to the AWOT dictionary.
    units : str
        Units associated with keyname dictionary. Optional.
    longname : str
        Long name associated with keyname dictionary. Optional.
    stdname : str
        Standard name associated with keyname dictionary. Optional.
    '''
    gcnames = ['gc', 'great_circle', 'great circle', 'g circle', 'great circ',
               'greatcircle', 'greatcirc']

    if lon_key is None:
        try:
            lon = awot['longitude']['data'][:]
        except:
            print("WARNING: Cannot find suitable longitude variable in file")
            return
    else:
        lon = awot[lon_key]['data'][:]

    if lat_key is None:
        try:
            lat = awot['latitude']['data'][:]
        except:
            print("WARNING: Cannot find suitable latitude variable in file")
            return
    else:
        lat = awot[lat_key]['data'][:]

    if method.lower().replace(" ", "") in gcnames:
        d = great_circle(lat, lon, radius=radius)
    else:
        d = vincentry(lat, lon, ellipsoid_key=ellipsoid_key,
                             iterations=iterations)
    # Mask any invalid entries
##    distance = np.ma.masked_invalid(distance)

    # Sum the distance over track
    dcum = np.cumsum(d)
    dcum = np.insert(dcum, 0, 0.)

    # Build a dictionary for output
    if units is None:
        units = 'meters'
    if longname is None:
        longname = "Track distance traveled over ground"
    if stdname is None:
        stdname = "Track ground distance"
    newdict = _build_dict(dcum, units, longname, stdname)

    # Add dictionary to AWOT dictionary if indicated
    if keyname is None:
        keyname = 'track_distance_ground'
    if add_to_dict:
        awot[keyname] = newdict
    return newdict


def calc_air_distance(awot, airspeed_key=None, time_key=None, add_to_dict=False,
                      keyname=None, units=None, longname=None,
                      stdname=None):
    '''
    Calculate the distance travelled in the air by multiplying air speed
    by time elapsed.

    The calculation assumes the SI standard units of meters per second for
    air speed and AWOT standard Epoch time.

    Parameters
    ----------
    awot : dict
        AWOT data dictionary instance.
    airspeed_key : str
        The key to use for air speed variable in AWOT dictionary.
    time_key : str
        The key to use for time variable in AWOT dictionary.
    add_to_dict : bool
        True to add the calculated distance to the AWOT dictionary.
    keyname : str
        The key to be used when adding to the AWOT dictionary.
    units : str
        The units associated with keyname dictionary. Optional.
    longname : str
        The long name associated with keyname dictionary. Optional.
    stdname : str
        The standard name associated with keyname dictionary. Optional.
    '''
    if airspeed_key is None:
        try:
            airspeed = awot['tas']['data'][:]
        except:
            print("WARNING: Cannot find suitable air speed variable in file")
            return
    else:
        airspeed = awot[airspeed_key]['data'][:]

    if time_key is None:
        try:
            timesec = date2num(awot['time']['data'][:], awot['time']['units'])
        except:
            print("WARNING: Cannot find suitable time variable in file")
            return
    else:
        timesec = date2num(awot[time_key]['data'][:], awot[time_key]['units'])
    time_elapsed = timesec[1:] - timesec[:-1]
    time_elapsed = np.insert(time_elapsed, 0, 0.)

    d = airspeed * time_elapsed
    # Mask any invalid entries
    d = np.ma.masked_invalid(d)

    # Sum the distance over track
    dcum = np.cumsum(d)

    # Build a dictionary for output
    if units is None:
        units = 'meters'
    if longname is None:
        longname = "Track distance traveled through air"
    if stdname is None:
        stdname = "Track air distance"
    newdict = _build_dict(dcum, units, longname, stdname)

    # Add dictionary to AWOT dictionary if indicated
    if keyname is None:
        keyname = 'track_distance_air'
    if add_to_dict:
        awot[keyname] = newdict
    return newdict


def great_circle(latitude, longitude, radius=None):
    """
    Use spherical geometry to calculate the surface distance between two
    geodesic points. This formula can be written many different ways,
    including just the use of the spherical law of cosines or the haversine
    formula.

    Parameters
    ----------
    latitude : array
        Array of floating point decimal latitude values.
        If want the distance between two points use (lat_start, lat_end).
    longitude : array
        Array of floating point decimal longitude values.
        If want the distance between two points use (lon_start, lon_end).
    radius : float
        Radius of earth. Default is to use the average great-circle radius
        of a sphere in meters. According to geopy documentation, using a
        sphere with this radius results in an error of up to about 0.5%.
    """
    if radius is None:
        radius = 6372795.
    if len(latitude) == 2 & len(longitude) == 2:
        latr1, lonr1 = np.radians(latitude[0]), np.radians(longitude[0])
        latr2, lonr2 = np.radians(latitude[1]), np.radians(longitude[1])
    else:
        latr1, lonr1 = np.radians(latitude[1:]), np.radians(longitude[1:])
        latr2, lonr2 = np.radians(latitude[:-1]), np.radians(longitude[:-1])

    sin_lat1, cos_lat1 = np.sin(latr1), np.cos(latr1)
    sin_lat2, cos_lat2 = np.sin(latr2), np.cos(latr2)

    delta_lon = lonr2 - lonr1
    cos_delta_lon, sin_delta_lon = np.cos(delta_lon), np.sin(delta_lon)

    d = np.arctan2(np.sqrt((cos_lat2 * sin_delta_lon) ** 2 +
                           (cos_lat1 * sin_lat2 -
                            sin_lat1 * cos_lat2 * cos_delta_lon) ** 2),
                   sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_delta_lon)

    return radius * d


def vincentry(latitude, longitude, ellipsoid_key=None, iterations=None):
    """
    Calculate the geodesic distance between two points using the formula
    devised by Thaddeus Vincenty, with an accurate ellipsoidal model of the
    earth.

    Set which ellipsoidal model of the earth to use by specifying an
    ``ellipsoid`` keyword argument. The default is 'WGS-84', which is the
    most globally accurate model.  If ``ellipsoid`` is a string, it is
    looked up in the `ELLIPSOIDS` dictionary to obtain the major and minor
    semiaxes and the flattening. Otherwise, it should be a tuple with those
    values.  See the comments above the `ELLIPSOIDS` dictionary for
    more information.

    Note: This implementation of Vincenty distance fails to converge for
    some valid points. In some cases, a result can be obtained by increasing
    the number of iterations (`iterations` keyword argument, given in the
    class `__init__`, with a default of 20). It may be preferable to use
    :class:`.great_circle`, which is marginally less accurate, but always
    produces a result.

    Parameters
    ----------
    latitude : array
        Array of floating point decimal latitude values.
        If want the distance between two points use (lat_start, lat_end).
    longitude : array
        Array of floating point decimal longitude values.
        If want the distance between two points use (lon_start, lon_end).
    ellipsoid_key : str
        Global ellipsoid model to use in calculations. 'WGS-84' is default.
    iterations : int
        Number ov iterations to try for convergent solution.
    """
    if ellipsoid_key is None:
        ellipsoid_key = 'WGS-84'
    if iterations is None:
        iterations = 20
    if len(latitude) == 2 & len(longitude) == 2:
        latr1, lonr1 = np.radians(latitude[0]), np.radians(longitude[0])
        latr2, lonr2 = np.radians(latitude[1]), np.radians(longitude[1])
    else:
        latr1, lonr1 = np.radians(latitude[0:-2]), np.radians(longitude[0:-2])
        latr2, lonr2 = np.radians(latitude[1:-1]), np.radians(longitude[1:-1])

    major, minor, f = ELLIPSOIDS[ellipsoid_key]

    dist = np.empty(len(latr1))

    delta_lon = lonr2 - lonr1

    # Compute reduced variables
    red_lat1 = np.arctan((1 - f) * np.tan(latr1))
    red_lat2 = np.arctan((1 - f) * np.tan(latr2))

    sin_red1, cos_red1 = np.sin(red_lat1), np.cos(red_lat1)
    sin_red2, cos_red2 = np.sin(red_lat2), np.cos(red_lat2)

    lam_lon = delta_lon
    lam_prime = 2 * np.pi

    for nn in range(len(latr1)):
        if np.isfinite(lam_lon[nn]):
            # Begin iterative calculation
            i = 0
            while abs(lam_lon[nn] - lam_prime) > 10e-12 and i <= iterations:
                i += 1

                sin_lam_lon = np.sin(lam_lon[nn])
                cos_lam_lon = np.cos(lam_lon[nn])

                sin_sig = np.sqrt((cos_red2[nn] * sin_lam_lon) ** 2 +
                                  (cos_red1[nn] * sin_red2[nn] -
                                   sin_red1[nn] * cos_red2[nn] *
                                   cos_lam_lon) ** 2)

                if sin_sig == 0:
                    print("WARNING: Coincident points found")
                    return 0

                cos_sig = (sin_red1[nn] * sin_red2[nn] + cos_red1[nn] *
                           cos_red2[nn] * cos_lam_lon)

                sig = np.arctan2(sin_sig, cos_sig)

                sin_alpha = (cos_red1[nn] * cos_red2[nn] *
                             sin_lam_lon / sin_sig)
                cos_sq_alpha = 1 - sin_alpha ** 2

                if cos_sq_alpha != 0:
                    cos2_sig_m = cos_sig - 2 * (sin_red1[nn] *
                                                sin_red2[nn] / cos_sq_alpha)
                else:
                    cos2_sig_m = 0.0 # Equatorial line

                C = f / 16. * cos_sq_alpha * (4 + f * (4 - 3 * cos_sq_alpha))

                lam_prime = lam_lon[nn]
                lam_lon[nn] = (delta_lon[nn] + (1 - C) * f * sin_alpha *
                           (sig + C * sin_sig * (cos2_sig_m + C * cos_sig *
                            (-1 + 2 * cos2_sig_m ** 2))))

            if i > iterations:
                raise ValueError("Vincenty formula failed to converge!")

            u_sq = cos_sq_alpha * (major ** 2 - minor ** 2) / minor ** 2

            A = 1 + u_sq / 16384. * (4096 + u_sq *
                                     (-768 + u_sq * (320 - 175 * u_sq)))

            B = u_sq / 1024. * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))

            delta_sig = (B * sin_sig * (cos2_sig_m + B / 4. *
                                        (cos_sig * (-1 + 2 * cos2_sig_m ** 2) -
                                         B / 6. * cos2_sig_m *
                                         (-3 + 4 * sin_sig ** 2) *
                                         (-3 + 4 * cos2_sig_m ** 2))))

            dist[nn] = minor * A * (sig - delta_sig)
        else:
            dist[nn] = np.nan
    return dist
