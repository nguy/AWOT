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
from netCDF4 import date2num
from ..io import common


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
        import time as timer
        appbegin = timer.time()

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

        # Create a field to store flight-matched data
        self.flight_matchdata = {}
        self.flight_matchdata['longitude'] = self.flight_lon
        self.flight_matchdata['latitude'] = self.flight_lat
        self.flight_matchdata['altitude'] = self.flight_alt
        self.flight_matchdata['time'] = self.flight_time
##        print("A --- %s seconds ---" % (timer.time() - appbegin))

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
##        print("B --- %s seconds ---" % (timer.time() - appbegin))

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
        print("C --- %s seconds ---" % (timer.time() - appbegin))

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
        # Create lists to contain the indices for each flight point
        ind1d = []
##        indy = []
##        indx = []
        distance = []

        # Calculate sines and cosines of lats/lon
        clat, clon = np.cos(self.latvalsr), np.cos(self.lonvalsr)
        slon, slat = np.sin(self.lonvalsr), np.sin(self.latvalsr)
        clats0, clons0 = np.cos(self.lat0valsr), np.cos(self.lon0valsr)
        slons0, slats0 = np.sin(self.lon0valsr), np.sin(self.lat0valsr)

        # Build kd-tree from big arrays of 3D coordinates
        if use_time:
            triples = list(zip(np.ravel(self.data_lon),
                               np.ravel(self.data_lat),
                               np.ravel(self.data_height),
                               date2num(np.ravel(self.data_time['data']),
                                        self.data_time['units'])))
        else:
            triples = list(zip(np.ravel(self.data_lon),
                               np.ravel(self.data_lat),
                               np.ravel(self.data_height)))
        kdt = scipy.spatial.cKDTree(triples, leafsize=leafsize)

        for ii in range(len(self.lat0valsr)):
            if use_time:
                dist_sq_min, minindex_1d = kdt.query(
                    [self.flight_lon['data'][ii], self.flight_lat['data'][ii],
                     self.flight_alt['data'][ii],
                     date2num(self.flight_time['data'][ii],
                              self.flight_time['units'])],
                    k=query_k, eps=query_eps, p=query_p,
                    distance_upper_bound=query_distance_upper_bound,
                    n_jobs=query_n_jobs)
            else:
                dist_sq_min, minindex_1d = kdt.query(
                    [self.flight_lon['data'][ii], self.flight_lat['data'][ii],
                     self.flight_alt['data'][ii]],
                    k=query_k, eps=query_eps, p=query_p,
                    distance_upper_bound=query_distance_upper_bound,
                    n_jobs=query_n_jobs)
            distance.append(dist_sq_min)
##            iy_min, ix_min = np.unravel_index(minindex_1d, self.data_lat.shape)
            ind1d.append(minindex_1d)
##            indy.append(iy_min)
##            indx.append(ix_min)
        indnd = np.unravel_index(ind1d, self.latvalsr.shape)

        if print_match_pairs:
            self._print_pair_by_index(indy, indx)

##        matchdata = self._get_data_by_index(indy, indx)
        matchdata = self._get_data_by_index(ind1d)
        return MatchData(self.flight_matchdata, matchdata, distance,
                         ind1d, indnd,
                         start_time=self.start_time, end_time=self.end_time)

    def near_neighbor_tunnel(self, use_time=False):
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
        radar : object
            A Py-ART radar object dervied using Py-ART read statement.
        basemap : object

        """
        # Create lists to contain the indices for each flight point
        indaz = []
        indrng = []
        indel = []
        self.matchdata = {}

        # First let's save the fields data into this class structure
        for field in radar.fields.keys():
                self.matchdata[field] = radar.fields[field].copy()
                self.matchdata[field]['data'] = np.ma.empty(len(self.lat0valsr))

        rlat, rlon = radar.latitude['data'][0], radar.longitude['data'][0]

        # Calculate the distance to the aircraft from radar
        dist_to_ac = self.calc_dist_to_aircraft(rlat, rlon,
                                                self.flight_lat['data'],
                                                self.flight_lon['data'])

        # Calculate the elevation angle corresponding to aircraft position
        ac_elev = np.degrees(np.arctan2(self.flight_alt['data'][:],
                                        dist_to_ac))

        ac_az = self.calc_az_to_aircraft(rlat, rlon,
                                         self.flight_lat['data'],
                                         self.flight_lon['data'], dist_to_ac)
        for ii in range(len(self.lat0valsr)):
            elin = np.argmin(np.abs(ac_elev[ii] - radar.fixed_angle['data']))
            radar_sweep = radar.extract_sweeps([elin])
            azin = np.argmin(np.abs(ac_az[ii]-radar_sweep.azimuth['data']))
            rgin = np.argmin(np.abs(dist_to_ac[ii] -
                                    radar_sweep.range['data']))#/1000.0))
            for field in self.data_fields.keys():
                self.matchdata[field]['data'][ii] = radar_sweep.fields[field]['data'][azin, rgin]

        return self.matchdata
#         tpos1 = sweep_latlon_to_flat_xy(radar, elin)
#         swpfields = radar_sweep.fields
#         ts_dz[index] = swpfields['reflectivity']['data'][azin, rgin]
#         ts_edr[index] = swpfields['turbulence']['data'][azin, rgin]
#         ts_vel[index] = swpfields['velocity']['data'][azin, rgin]
#         ts_sw[index] = swpfields['spectrum_width']['data'][azin, rgin]

###################
#  Calculations   #
###################

    def calc_dist_to_aircraft(self, lat0, lon0, aclats, aclons, R=6371000.):
        '''
        Calculate the distance to the aircraft from some point using
        the Haversine formula.

        Parameters
        ----------
        lat0 : float
            Latitude value of origin point to calculate from.
        lon0 : float
            Longitude value of origin point to calculate from.
        aclats : float array
            Array of aircraft latitude values.
        aclons : float array
            Array of aircraft longitude values.
        R : float, optional
            Earth radius in meters.
        '''
        delLat = lat0 - aclats
        delLon = lon0 - aclons

        a = ((np.sin(np.radians(delLat)/2))**2 + np.cos(np.radians(lat0)) *
             np.cos(np.radians(aclats)) * (np.sin(np.radians(delLon)/2))**2)
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
        distance = R * c
        return distance

    def calc_az_to_aircraft(self, lat0, lon0, aclats, aclons, rng, R=6371000.):
        '''
        Calculate the azimuth to the aircraft.

        Parameters
        ----------
        lat0 : float
            Latitude value of origin point to calculate from.
        lon0 : float
            Longitude value of origin point to calculate from.
        aclats : float array
            Array of aircraft latitude values.
        aclons : float array
            Array of aircraft longitude values.
        R : float, optional
            Earth radius in meters.
        '''
        factor = 2 * np.pi * R / 360.

        top = np.sin(aclats) - np.cos(rng/factor) * np.sin(lat0)
        bottom = np.sin(rng/factor) * np.cos(lat0)

        azimuth = np.arccos(top / bottom)
        return azimuth

####################
#  Print methods   #
####################

    def _print_pair_by_index(self, indlon, indlat):
        for ii in range(len(indlon)):
            print(("AC Lat: %g, Lon: %g, Alt: %g | "
                   "Rad Lat: %g, Lon: %g, Alt: %g") % (
                  self.flight_lon['data'][ii], self.flight_lat['data'][ii],
                  self.flight_alt['data'][ii], self.data_lon[iy_min, ix_min],
                  self.data_lat[iy_min, ix_min],
                  self.data_height[iy_min, ix_min]))

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
