"""
awot.util.matcher
=================

Matching function for airborne and ground-based radar data.

Original version provided by George Duffy 'APRmacher2.py'.
This code looks to generalize it to a structure that works with
AWOT objects.

This code should be extendable to any input data set.
"""
import numpy as np
import scipy
import datetime
from netCDF4 import date2num, num2date
import time as timer

from ..io import common
from ..graph import common as gcommon


class TrackMatch(object):
    """Class for matching flight level data to other observations."""

    def __init__(self, flight, datavolume,
                 start_time=None, end_time=None,
                 data_lon=None, data_lat=None,
                 data_height=None, data_time=None,
                 flight_lon_name=None, flight_lat_name=None,
                 flight_alt_name=None, flight_time_name=None,
                 flight_uwind_name=None, flight_vwind_name=None,
                 field_match_dict=None,
                 ):
        '''
        Intitialize the class to create plots.

        Parameters
        ----------
        flight : dict
            AWOT flight data dictionary object.
        datavolume : dict
            AWOT dictionary object of some data field. The nearest points
            in this field are matched to flight track points.

        Optional
        --------
        start_time : str
            UTC time to use as start time for subsetting in datetime format.
            (e.g. 2014-08-20 12:30:00)
        end_time : str
            UTC time to use as an end time for subsetting in datetime format.
            (e.g. 2014-08-20 16:30:00)
        data_lon : flt array
            Ndarray longitude [deg] array same size as data.
        data_lat : flt array
            Ndarray latitude [deg] array same size as data.
        data_height : flt array
            Ndarray altitude [meters] array same size as data.
        data_time : dict
            Time array in datetime format.
        flight_lon_name : str
            Variable name to use for longitude array.
            None uses AWOT default mapping.
        flight_lat_name : str
            Variable name to use for latitude array.
            None uses AWOT default mapping.
        flight_alt_name : str
            Variable name to use for altitude array.
            None uses AWOT default mapping.
        flight_time_name : str
            Variable name to use for time array.
            None uses AWOT default mapping.
        flight_uwind_name : str
            Variable name to use for zonal wind array.
            None uses AWOT default mapping.
        flight_vwind_name : str
            Variable name to use for meridional wind array.
            None uses AWOT default mapping.
        field_match_dict : list
            List of field names to match for data object.
            None loops through existing fields.
        '''
        appbegin = timer.time() ##NG

        self.flight_data = flight
        self.datavol = datavolume
        self.time = flight['time']

        # Check if keyword names were given
        if flight_lon_name is None:
            lonkey = 'longitude'
        else:
            lonkey = lon_name
        if flight_lat_name is None:
            latkey = 'latitude'
        else:
            latkey = lat_name
        if flight_alt_name is None:
            altkey = 'altitude'
        else:
            altkey = alt_name
        if flight_time_name is None:
            timekey = 'time'
        else:
            timekey = time_name
        if flight_uwind_name is None:
            ukey = 'Uwind'
        else:
            ukey = uwind_name
        if flight_vwind_name is None:
            vkey = 'Vwind'
        else:
            vkey = vwind_name

        # Set the required flight fields
        self.flight_lon = self._get_var_time_subset(
                             flight[lonkey].copy(),
                             start_time, end_time)
        self.flight_lat = self._get_var_time_subset(
                             flight[latkey].copy(),
                             start_time, end_time)
        self.flight_alt = self._get_var_time_subset(
                             flight[altkey].copy(),
                             start_time, end_time)
        if flight[ukey] is not None:
            self.flight_Uwind = self._get_var_time_subset(
                                 flight[ukey].copy(),
                                 start_time, end_time)
        if flight[vkey] is not None:
            self.flight_Vwind = self._get_var_time_subset(
                                 flight[vkey].copy(),
                                 start_time, end_time)

        # Since time could come in in different datetime units,
        # we need to convert to common epoch times
        flight_epoch = common.convert_to_epoch_dict(flight[timekey])
        self.flight_time = self._get_var_time_subset(
                             flight_epoch,
                             start_time, end_time)
        self.flight_numtime = date2num(
            self.flight_time['data'], self.flight_time['units'])

        # Create a field to store flight-matched data
        self.flight_matchdata = {}
        self.flight_matchdata['longitude'] = self.flight_lon
        self.flight_matchdata['latitude'] = self.flight_lat
        self.flight_matchdata['altitude'] = self.flight_alt
        self.flight_matchdata['time'] = self.flight_time

        for field in flight.keys():
            try:
                self.flight_matchdata[field] = self._get_var_time_subset(
                                          flight[field].copy(),
                                          start_time, end_time)
            except:
                self.flight_matchdata[field] = None

        # Create a field to store data field to matching data
        self.data_fields = {}
        if field_match_dict is None:
            for field in self.datavol['fields'].keys():
                self.data_fields[field] = self.datavol['fields'][field].copy()
        else:
            for field in field_match_dict:
                self.data_fields[field] = self.datavol['fields'][field].copy()
        dshape = self.data_fields[self.data_fields.keys()[0]]['data'].shape

        # Run checks to make sure array lengths are the same
        warntxt = "Check that data is an AWOT object!"
        if data_lon is None:
            try:
                self.data_lon = self.datavol['longitude']['data']
            except:
                print("Cannot find longitude, %s" % warntxt)
        else:
            self.data_lon = np.array(data_lon)

        if data_lat is None:
            try:
                self.data_lat = self.datavol['latitude']['data']
            except:
                print("Cannot find latitude, %s" % warntxt)
        else:
            self.data_lat = np.array(data_lat)

        if data_height is None:
            try:
                self.data_height = self.datavol['height']['data']
            except:
                print("Cannot find height, %s" % warntxt)
        else:
            self.data_height = np.array(data_height)

        if data_time is None:
            try:
                self.data_time = common.convert_to_epoch_dict(self.datavol['time'])
            except:
                self.data_time = common._get_epoch_dict(
                    self.datavol['time']['data'],
                    self.datavol['time']['units'])
            else:
                print("Cannot find time, %s" % warntxt)
        else:
            try:
                self.data_time = common.convert_to_epoch_dict(data_time)
            except:
                self.data_time = common._get_epoch_dict(
                                    data_time['data'],
                                    data_time['units'])
        self.data_numtime = date2num(self.data_time['data'][:], self.time['units'])
        print("A --- %s seconds ---" % (timer.time() - appbegin)) ##NG

        # Convert time to same shape as fields if not already
        if len(self.data_numtime.shape) == 1:
            self.data_numtime = np.resize(self.data_numtime, dshape)
        print("B --- %s seconds ---" % (timer.time() - appbegin)) #NG

        self.start_time = start_time
        self.end_time = end_time

        # Convert latitude and longitude to radians
        self.latvalsr = np.radians(self.data_lat)
        self.lonvalsr = np.radians(self.data_lon)
        self.lat0valsr = np.radians(self.flight_lat['data'][:])
        self.lon0valsr = np.radians(self.flight_lon['data'][:])

    def kdtree(self, leafsize=16,
               use_time=False, print_match_pairs=False,
               query_k=1, query_eps=0, query_p=2,
               query_distance_upper_bound=np.inf,
               query_n_jobs=1):
        '''
        Find the closest point using a K-Dimensional Tree method.
        The tree building can take a long time depending on how
        large the input data array.

        Parameters
        ----------
        leafsize : int
            Positive number of points at which time
            scipy.spatial.cKDTree switches to brute force.
            Default same as scipy default.
        use_time: bool
            True includes the time variables in KD-Tree.
            Default is False.
        print_match_pairs: bool
            True returns screen print out of pair results.
            Default is False.
        query_?? : See scipy.spatial.cKDTree.query
            Default keyword values used in the cKDTree.query.
        '''
        kdbegin = timer.time() ##NG
        # Create lists to contain the indices for each flight point
##        ind1d = []
##        distance = []

        # Calculate sines and cosines of lats/lon
        clat, clon = np.cos(self.latvalsr), np.cos(self.lonvalsr)
        slon, slat = np.sin(self.lonvalsr), np.sin(self.latvalsr)
        clats0, clons0 = np.cos(self.lat0valsr), np.cos(self.lon0valsr)
        slons0, slats0 = np.sin(self.lon0valsr), np.sin(self.lat0valsr)

        dlon, dlat = np.ravel(self.data_lon), np.ravel(self.data_lat)
        dht, dti = np.ravel(self.data_height), np.ravel(self.data_numtime)
        print("C --- %s seconds ---" % (timer.time() - kdbegin)) ##NG

        # Build kd-tree from big arrays of 3D coordinates
#         if use_time:
#             triples = list(zip(np.ravel(self.data_lon),
#                                np.ravel(self.data_lat),
#                                np.ravel(self.data_height),
#                                np.ravel(self.data_numtime)))
# ##                               date2num(np.ravel(self.data_time['data']),
# ##                                        self.data_time['units'])))
#         else:
#             triples = list(zip(np.ravel(self.data_lon),
#                                np.ravel(self.data_lat),
#                                np.ravel(self.data_height)))
#         kdt = scipy.spatial.cKDTree(triples, leafsize=leafsize)
        if use_time:
            kdt = scipy.spatial.cKDTree(zip(dlon, dlat, dht, dti), leafsize=leafsize)
        else:
            kdt = scipy.spatial.cKDTree(zip(dlon, dlat, dht), leafsize=leafsize)
        print("D --- %s seconds ---" % (timer.time() - kdbegin)) ##NG
        if use_time:
            distance, ind1d = kdt.query(
                zip(self.flight_lon['data'], self.flight_lat['data'],
                self.flight_alt['data'], self.flight_numtime),
                k=query_k, eps=query_eps, p=query_p,
                distance_upper_bound=query_distance_upper_bound,
                n_jobs=query_n_jobs)
        else:
            distance, ind1d = kdt.query(
                zip(self.flight_lon['data'], self.flight_lat['data'],
                self.flight_alt['data']),
                k=query_k, eps=query_eps, p=query_p,
                distance_upper_bound=query_distance_upper_bound,
                n_jobs=query_n_jobs)
        print("E --- %s seconds ---" % (timer.time() - kdbegin)) ##NG

#         for ii in range(len(self.lat0valsr)):
#             if use_time:
#                 dist_sq_min, minindex_1d = kdt.query(
#                     [self.flight_lon['data'][ii], self.flight_lat['data'][ii],
#                      self.flight_alt['data'][ii],
#                      self.flight_numtime[ii]],
# ##                     date2num(self.flight_time['data'][ii],
# ##                              self.flight_time['units'])],
#                     k=query_k, eps=query_eps, p=query_p,
#                     distance_upper_bound=query_distance_upper_bound,
#                     n_jobs=query_n_jobs)
#             else:
#                 dist_sq_min, minindex_1d = kdt.query(
#                     [self.flight_lon['data'][ii], self.flight_lat['data'][ii],
#                      self.flight_alt['data'][ii]],
#                     k=query_k, eps=query_eps, p=query_p,
#                     distance_upper_bound=query_distance_upper_bound,
#                     n_jobs=query_n_jobs)
#             distance.append(dist_sq_min)
#             ind1d.append(minindex_1d)
        indnd = np.unravel_index(ind1d, self.latvalsr.shape)

        if print_match_pairs:
            self._print_pair_by_index(ind1d)

##        matchdata = self._get_data_by_index(indy, indx)
        matchdata = self._get_data_by_index(ind1d)
        return MatchData(self.flight_matchdata, matchdata, distance,
                         ind1d, indnd,
                         start_time=self.start_time, end_time=self.end_time)

    def near_neighbor_tunnel(self, use_time=False, print_match_pairs=False):
        '''
        Find closest point to the set of (lat,lon) points
        provided by the flight data object to data object points.
        This routine was adapted from:
        http://nbviewer.ipython.org/github/Unidata/tds-python-workshop/blob/master/netcdf-by-coordinates.ipynb

        Returns iy,ix such that the square of the tunnel distance
        between (latval[it,ix],lonval[iy,ix]) and (lat0,lon0)
        is minimum.

        Parameters
        ----------
        use_time: bool
            True to limit results to closest time value.
            Default is False.
        '''
        nntbegin = timer.time() ##NG
        # Create lists to contain the indices for each flight point
        ind1d = []
##        indx = []
##        indy = []
        distance = []

        # Calculate sines and cosines of lats/lon
        clat, clon = np.cos(self.latvalsr), np.cos(self.lonvalsr)
        slon, slat = np.sin(self.lonvalsr), np.sin(self.latvalsr)
        clats0, clons0 = np.cos(self.lat0valsr), np.cos(self.lon0valsr)
        slons0, slats0 = np.sin(self.lon0valsr), np.sin(self.lat0valsr)

        for ii in range(len(self.lat0valsr)):
            delX = clats0[ii] * clons0[ii] - clat * clon
            delY = clats0[ii] * slons0[ii] - clat * slon
            delZ = slats0[ii] - slat
            # Calculate the distance squared
            dist_sq = delX**2 + delY**2 + delZ**2
            # Find the 1D index of the minimum element
            minindex_1d = dist_sq.argmin()
            distance.append(dist_sq.min())
##            iy_min, ix_min = np.unravel_index(minindex_1d, self.latvalsr.shape)
            ind1d.append(minindex_1d)
##            indy.append(iy_min)
##            indx.append(ix_min)
        indnd = np.unravel_index(ind1d, self.latvalsr.shape)

        if print_match_pairs:
            self._print_pair_by_index(ind1d)

        matchdata = self._get_data_by_index(ind1d)
        return MatchData(self.flight_matchdata, matchdata, distance,
                         ind1d, indnd,
                         start_time=self.start_time, end_time=self.end_time)
##        matchdata = self._get_data_by_index(indy, indx)
##        return MatchData(self.flight_matchdata, matchdata, distance,
##                         ind1d, indx, indy,
##                         start_time=self.start_time, end_time=self.end_time)

    def near_neighbor_pyart(self, radar, basemap=None):
        """
        Calculate variables using a nearest neighbor approach when the input
        data is a Py-ART radar object.

        Parameters
        ----------
        radar : object or list of objects
            A Py-ART radar object dervied using Py-ART read statement.
            If a list is of Py-ART radar objects is supplied the function
            loops through each instance looking for nearest neighbor points.
        basemap : object

        """
        nnpbegin = timer.time() ##NG
        # Create lists to contain the indices for each flight point
        indaz = []
        indrng = []
        indel = []
        self.matchdata = {}

        # Check if list is passed, if single build a list of 1
        if type(radar).__name__ == 'list':
            rlist = radar
        else:
            rlist = [radar]

        # First let's save the fields data into this class structure
        # We assume that all files have the same setup
        for field in rlist[0].fields.keys():
                self.matchdata[field] = rlist[0].fields[field].copy()
                self.matchdata[field]['data'] = np.ma.empty(len(self.lat0valsr))

        for num, pr in enumerate(rlist):
            rlat, rlon = pr.latitude['data'][0], pr.longitude['data'][0]

            # Calculate the distance to the aircraft from radar
#            ac_dist = self.distance_to_point(
#                rlat, rlon, self.flight_lat['data'], self.flight_lon['data'])
            ac_dist = self.distance_to_point2(
                rlat, rlon, self.flight_lat['data'], self.flight_lon['data'])

            # Calculate the elevation angle corresponding to aircraft position
##            ac_elev = np.degrees(np.arctan2(self.flight_alt['data'][:],
##                                            ac_dist))
##            ac_az = self.azimuth_to_point(
##                rlat, rlon, self.flight_lat['data'],
##                self.flight_lon['data'], ac_dist))
            ac_az = self.bearing_to_point(
                rlat, rlon, self.flight_lat['data'],
                self.flight_lon['data'])
            ac_rng, ac_elev = self.slant_range_and_elev(
                ac_dist, self.flight_alt['data'][:])

            # Convert the radar time to AWOT epoch
            rtime = common.convert_to_epoch_dict(pr.time.copy())
            rnumtime = date2num(rtime['data'], rtime['units'])
##            pr_start = datetime.datetime.strptime(
##                pr.time['units'][14:], '%Y-%m-%dT%H:%M:%SZ')
##        # Subtracting 35 seconds (half the volume update time)
##        # Treating volume start as midpoint for aircraft comparison
##    std = pr_start + dt.timedelta(seconds=pr.time['data'][0]-35)
##    edd = pr_start + dt.timedelta(seconds=pr.time['data'][-1]-35)
##    cond = np.logical_and(fl['time']['data'] >= std, fl['time']['data'] < edd)

            cond = np.logical_and(self.flight_numtime >= rnumtime[0],
                                  self.flight_numtime <= rnumtime[-1])
            indices = np.where(cond)

            if np.size(indices[0]) > 0:
##                for i, index in enumerate(indices[0]):
##                    elin = np.argmin(np.abs(kel[index]-radar2.fixed_angle['data']))
##                    radar_sweep = radar2.extract_sweeps([elin])
##                    azin = np.argmin(np.abs(kbear[index]-radar_sweep.azimuth['data']))
##                    rgin = np.argmin(np.abs(ksr[index]-radar_sweep.range['data']/1000.0))

                for index in indices[0]:
                    elin = np.argmin(np.abs(ac_elev[index] - pr.fixed_angle['data']))
                    pr_sweep = pr.extract_sweeps([elin])
                    azin = np.argmin(np.abs(ac_az[index]-pr_sweep.azimuth['data']))
#                    rgin = np.argmin(np.abs(ac_dist[index] -
#                                    pr_sweep.range['data']))
                    rgin = np.argmin(np.abs(ac_rng[index] -
                                    pr_sweep.range['data']))
#                    print(azin, rgin)
#                    print("AC/Radar lat/lon/alt: %g/%g, %g/%g, %g/%g" % (
#                        self.flight_lat['data'][index],
#                        pr.gate_latitude['data'][azin, rgin],
#                        self.flight_lon['data'][index],
#                        pr.gate_longitude['data'][azin, rgin],
#                        self.flight_alt['data'][index],
#                        pr.gate_altitude['data'][azin, rgin]))
                    for field in pr_sweep.fields.keys():
                        self.matchdata[field]['data'][index] = pr_sweep.fields[field]['data'][azin, rgin]
            print("B --- %s seconds for loop %d---" % (timer.time() - nnpbegin, num)) ##NG

        return self.matchdata

###################
#  Calculations   #
###################

    def distance_to_point(self, lat0, lon0, aclats, aclons, R=6371100.):
        '''
        Calculate the distance to the aircraft from some point using
        the Haversine formula.

        Parameters
        ----------
        lat0 : float
            Latitude [deg] value of origin point from which to calculate.
        lon0 : float
            Longitude [deg] value of origin point from which to calculate.
        aclats : float array
            Array of aircraft latitude [deg] values.
        aclons : float array
            Array of aircraft longitude [deg] values.
        R : float, optional
            Earth radius in meters.
        '''
        lat0r = np.radians(lat0)
        lon0r = np.radians(lon0)
        aclatsr = np.radians(aclats)
        aclonsr = np.radians(aclons)
        dlat = lat0r - aclatsr
        dlon = lon0r - aclonsr

        a = ((np.sin(dlat/2.))**2 + np.cos(lat0r) *
             np.cos(aclatsr) * (np.sin(dlon/2.))**2)
        c = 2. * np.arctan2(np.sqrt(a), np.sqrt(1. - a))
        return R * c

    def distance_to_point2(self, lat0, lon0, aclats, aclons, R=6371100.):
        '''
        Calculate the distance to the aircraft.

        Parameters
        ----------
        lat0 : float
            Latitude [deg] value of origin point from which to calculate.
        lon0 : float
            Longitude [deg] value of origin point from which to calculate.
        aclats : float array
            Array of aircraft latitude [deg] values.
        aclons : float array
            Array of aircraft longitude [deg] values.
        R : float, optional
            Earth radius in meters.
        '''
        lat0r = np.radians(lat0)
        lon0r = np.radians(lon0)
        aclatsr = np.radians(aclats)
        aclonsr = np.radians(aclons)
        dlat = lat0r - aclatsr
        dlon = lon0r - aclonsr
        return (np.arccos(np.sin(aclatsr) * np.sin(lat0r) +
            np.cos(aclatsr) * np.cos(lat0r) * np.cos(lon0r - aclonsr)) * R)

    def azimuth_to_point(self, lat0, lon0, aclats, aclons, rng, R=6371100.):
        '''
        Calculate the azimuth to the aircraft.

        Parameters
        ----------
        lat0 : float
            Latitude value of origin point from which to calculate.
        lon0 : float
            Longitude value of origin point from which to calculate.
        aclats : float array
            Array of aircraft latitude values.
        aclons : float array
            Array of aircraft longitude values.
        R : float, optional
            Earth radius in meters.
        '''
        factor = 2. * np.pi * R / 360.
        top = np.sin(np.radians(aclats)) - np.cos(rng/factor) * np.sin(np.radians(lat0))
        bottom = np.sin(rng/factor) * np.cos(np.radians(lat0))
        azimuth = np.degrees(np.arccos(top / bottom))
        try:
            azimuth[azimuth < 0.] = azimuth[azimuth < 0.] + 360.
        except:
            pass
        return azimuth

    def bearing_to_point(self, lat0, lon0, aclats, aclons):
        '''
        Calculate the bearing to the aircraft.
        Method was adopted from PyTDA by Timothy Lang.

        Parameters
        ----------
        lat0 : float
            Latitude value of origin point from which to calculate.
        lon0 : float
            Longitude value of origin point from which to calculate.
        aclats : float array
            Array of aircraft latitude values.
        aclons : float array
            Array of aircraft longitude values.
        '''
        aclonsr = np.radians(aclons)
        aclatsr = np.radians(aclats)
        lon0r = np.radians(lon0)
        lat0r = np.radians(lat0)

        bear = np.arctan2(
            (np.sin(aclonsr - lon0r) * np.cos(aclatsr)),
            (np.cos(lat0r) * np.sin(aclatsr) -
             np.sin(lat0r) * np.cos(aclatsr) * np.cos(aclonsr - lon0r)))
        return np.degrees(bear)

    def slant_range_and_elev(self, ground_range, height, R=6371100.):
        '''
        Calculate slant range and elevation.
        Method was adopted from PyTDA by Timothy Lang.

        Parameters
        ----------
        ground_range: float or array
            Ground range [m]
        height: float or array
            Height [m]
        R : float, optional
            Earth radius in meters.
        '''
        Re = 4.0/3.0 * R  # Effective earth radius.
        rh = height + Re
        slantrsq = Re**2 + rh**2 - (2 * Re * rh * np.cos(ground_range/Re))
        slantr = np.sqrt(slantrsq)
        elev = np.arccos((Re**2 + slantrsq - rh**2)/(2 * Re * slantr))
        elev = elev * 180.0/np.pi
        elev = elev - 90.0
        return slantr, elev

####################
#  Print methods   #
####################

    def _print_pair_by_index(self, ind1d):
        for ii, ind in enumerate(ind1d):
            print(("AC Lat: %g, Lon: %g, Alt: %g | "
                   "Rad Lat: %g, Lon: %g, Alt: %g") % (
                  self.flight_lon['data'][ii],
                  self.flight_lat['data'][ii],
                  self.flight_alt['data'][ii],
                  np.ravel(self.data_lon)[ind],
                  np.ravel(self.data_lat)[ind],
                  np.ravel(self.data_height)[ind]))
            print(
            "AC/Rad Turb: %g / %g" % (self.flight_matchdata['turb']['data'][ii],
                                      np.ravel(self.data_fields['turbulence']['data'])[ind]))

##################
#  Get methods   #
##################

    def _get_datetime(self, time_string, get_start=False, get_end=False):
        '''Get a start time as datetime instance for subsetting.'''
        # Check to see if time is subsetted
        if time_string is None:
                if get_start is True:
                    dt = self.time['data'][:].min()
                if get_end is True:
                    dt = self.time['data'][:].max()
        else:
            tStr = [time_string[0:4], time_string[5:7], time_string[8:10],
                    time_string[11:13], time_string[14:16],
                    time_string[17:19], '0']
            tInt = [int(s) for s in tStr]
            try:
                dt = datetime.datetime(tInt[0], tInt[1], tInt[2], tInt[3],
                                       tInt[4], tInt[5], tInt[6])
            except:
                print("Check the format of date string "
                      "(e.g. '2014-08-20 12:30:00')")
                return

        return dt

    def _get_var_time_subset(self, var, start_time, end_time):
        '''Variable subsetted by time.'''
        # Check to see if time is subsetted
        dt_start = self._get_datetime(start_time, get_start=True)
        dt_end = self._get_datetime(end_time, get_end=True)
        vsub = var.copy()
        vsub['data'] = vsub['data'][(self.time['data'][:] >= dt_start) &
                                    (self.time['data'][:] <= dt_end)]
        return vsub

#     def _get_data_by_index(self, indlon, indlat):
#         '''Retrieve data fields using calculated indices.'''
#         data_by_ind = {}
#         for field in self.data_fields.keys():
#             dfield = self.data_fields[field]
#             data_by_ind[field] = dfield.copy()
#             data_by_ind[field]['data'] = dfield['data'][indlon, indlat]
#         return data_by_ind

    def _get_data_by_index(self, indices_1d):
        '''Retrieve data fields using calculated indices.'''
        data_by_ind = {}
        for field in self.data_fields.keys():
            dfield = self.data_fields[field]
            data_by_ind[field] = dfield.copy()
            data_by_ind[field]['data'] = np.ravel(dfield['data'])[indices_1d]
        return data_by_ind

    def _get_data_by_pyart_index(self, radar, indel, indaz, indrng):
        '''Retrieve data fields using calculated Py-ART indices.'''
        data_by_ind = {}
        for field in radar.fields.keys():
            rfield = radar.fields[field]
            data_by_ind[field] = rfield.copy()
            data_by_ind[field]['data'] = rfield['data'][azin, rgin]
        return data_by_ind


class MatchData(object):
    """Class for storing matched flight level data to other observations."""

    def __init__(self, flight, data, distance_to_point, indices_1d, indices_nd,
                 start_time=None, end_time=None):
##                 indices_x, indices_y, start_time=None, end_time=None):
        '''

        Parameters
        ----------
        flight: dict
            Dictionary of flight data variables.
        data: dict
            Dictionary of matched data points.
        distance_to_point: array
            Distance of nearest neighbor to flight point.
        indices_1d: list
            List of 1D indices of nearest neighbor.
        indices_x: array
            Longitudinal indices of nearest neighbor.
        indices_y: array
            Latitudinal indices of nearest neighbor.
        start_time: datetime
            Start of nearest neighbor times.
        end_time: datetime
            End of nearest neighbor times.
        '''
        self.flight = flight
        self.data = data
        self.distance_to_point = distance_to_point
        self.indices_1d = indices_1d
        self.indices_nd = indices_nd
##        self.lon_indices = lon_indices
##        self.lat_indices = lat_indices
        self.start_time = start_time
        self.end_time = end_time
        return
