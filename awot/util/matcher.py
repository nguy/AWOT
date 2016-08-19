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
from scipy.spatial import cKDTree
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
                 data_time=None,
                 flight_lon_name=None, flight_lat_name=None,
                 flight_alt_name=None, flight_time_name=None,
                 flight_uwind_name=None, flight_vwind_name=None,
                 field_match_dict=None,
                 ):
        '''
        Intitialize the class to create plots.
        It is assumed that data fields have the same shape.

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
        # Set variables to be used later by any function
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
        try:
            if flight[ukey] is not None:
                self.flight_Uwind = self._get_var_time_subset(
                    flight[ukey].copy(), start_time, end_time)
        except:
            pass
        try:
            if flight[vkey] is not None:
                self.flight_Vwind = self._get_var_time_subset(
                flight[vkey].copy(), start_time, end_time)
        except:
            pass

        # Since time could come in in different datetime units,
        # we need to convert to AWOT common epoch times
        flight_epoch = common.convert_to_epoch_dict(flight[timekey])
        self.flight_time = self._get_var_time_subset(
                             flight_epoch,
                             start_time, end_time)
        self.flight_numtime = date2num(
            self.flight_time['data'], self.flight_time['units'])

        # Create a field to store flight-matched data
        self.fldata = {}
        self.fldata['longitude'] = self.flight_lon
        self.fldata['latitude'] = self.flight_lat
        self.fldata['altitude'] = self.flight_alt
        self.fldata['time'] = self.flight_time

        for field in flight.keys():
            try:
                self.fldata[field] = self._get_var_time_subset(
                    flight[field].copy(), start_time, end_time)
            except:
                self.fldata[field] = None

        # Create a field to store data field to matching data
        self.data_fields = {}
        if field_match_dict is None:
            for field in self.datavol['fields'].keys():
                self.data_fields[field] = self.datavol['fields'][field].copy()
        else:
            for field in field_match_dict:
                self.data_fields[field] = self.datavol['fields'][field].copy()
        self.data_shape = self.data_fields[
            self.data_fields.keys()[0]]['data'].shape

        # Run checks to make sure array lengths are the same
        warntxt = "Check that data is an AWOT style object or specify keyword!"
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

        if data_time is None:
            if self.datavol['time'] is None:
                self.data_time = None
                print("Cannot find time, %s" % warntxt)
            else:
                try:
                    self.data_time = common.convert_to_epoch_dict(
                        self.datavol['time'])
                except:
                    self.data_time = common._get_epoch_dict(
                        self.datavol['time']['data'],
                        self.datavol['time']['units'])
                else:
                    self.data_time = None
                    print("Cannot find time, %s" % warntxt)
        else:
            try:
                self.data_time = common.convert_to_epoch_dict(data_time)
            except:
                self.data_time = common._get_epoch_dict(
                                    data_time['data'],
                                    data_time['units'])
        if self.data_time is not None:
            self.data_numtime = date2num(
                self.data_time['data'][:], self.time['units'])

            # Convert time to same shape as fields if not already
            if len(self.data_numtime.shape) == 1:
                self.data_numtime = np.resize(
                    self.data_numtime, self.data_shape)

        self.start_time = start_time
        self.end_time = end_time

        # Convert latitude and longitude to radians
        self.latvalsr = np.radians(self.data_lat)
        self.lonvalsr = np.radians(self.data_lon)
        self.lat0valsr = np.radians(self.flight_lat['data'][:])
        self.lon0valsr = np.radians(self.flight_lon['data'][:])

    def kdtree(self, leafsize=16,
               use_time=False, verbose=False, timeit=True,
               query_k=1, query_eps=0, query_p=2,
               query_distance_upper_bound=np.inf,
               query_n_jobs=1):
        '''
        Find the closest point using a K-Dimensional Tree method.
        The tree building can take a long time depending on how
        large the input data array. Subsetting is attemped to increase
        execution speed.

        Parameters
        ----------
        leafsize : int
            Positive number of points at which time
            scipy.spatial.cKDTree switches to brute force.
            Default same as scipy default.
        use_time: bool
            True includes the time variables in KD-Tree. Defaults False.
        verbose: bool
            True returns screen print out of pair results. Defaults False.
        timeit: bool
            True prints elapsed time of operation.
        query_?? : See scipy.spatial.cKDTree.query
            Default keyword values used in the cKDTree.query.
        '''
        opbegin = timer.time()
        # Create 1-D variables to use
        dlon, dlat = np.ravel(self.data_lon), np.ravel(self.data_lat)

        # Subset input data to only include relevant points within
        condlon = np.logical_and(dlon >= np.min(self.flight_lon['data']),
                                 dlon <= np.max(self.flight_lon['data']))
        condlat = np.logical_and(dlat >= np.min(self.flight_lat['data']),
                                 dlat <= np.max(self.flight_lat['data']))

        if use_time:
            dti = np.ravel(self.data_numtime)
            condti = np.logical_and(dti >= np.min(self.flight_numtime),
                                    dti <= np.max(self.flight_numtime))

            cond = ((dlon >= np.min(self.flight_lon['data'])) &
                    (dlon <= np.max(self.flight_lon['data'])) &
                    (dlat >= np.min(self.flight_lat['data'])) &
                    (dlat <= np.max(self.flight_lat['data'])) &
                    (dti >= np.min(self.flight_numtime)) &
                    (dti <= np.max(self.flight_numtime)))
        else:
            cond = ((dlon >= np.min(self.flight_lon['data'])) &
                    (dlon <= np.max(self.flight_lon['data'])) &
                    (dlat >= np.min(self.flight_lat['data'])) &
                    (dlat <= np.max(self.flight_lat['data'])))

        indices = np.where(cond)
        print(len(indices[0]))

        if len(indices[0]) > 1:
            dlon = dlon[indices[0]]
            dlat = dlat[indices[0]]
            if use_time:
                dti = dti[indices[0]]
        else:
            print("No points found for condition")

        # Build and query kd-tree
        if use_time:
            kdt = cKDTree(zip(dlon, dlat, dti), leafsize=leafsize)
        else:
            kdt = cKDTree(zip(dlon, dlat), leafsize=leafsize)
        if use_time:
            distance, ind1d = kdt.query(
                zip(self.flight_lon['data'], self.flight_lat['data'],
                    self.flight_numtime),
                k=query_k, eps=query_eps, p=query_p,
                distance_upper_bound=query_distance_upper_bound,
                n_jobs=query_n_jobs)
        else:
            distance, ind1d = kdt.query(
                zip(self.flight_lon['data'], self.flight_lat['data']),
                k=query_k, eps=query_eps, p=query_p,
                distance_upper_bound=query_distance_upper_bound,
                n_jobs=query_n_jobs)
#        indnd = np.unravel_index(ind1d, self.latvalsr.shape)

        if verbose:
            self._print_pair_by_index(ind1d)

        matchdata = self._get_data_by_index(ind1d)
        if timeit:
            print("--- Elapsed time: %s sec ---" % (timer.time() - opbegin))
        return MatchData(self.fldata, matchdata, distance_to_point=distance,
                         indices_1d=ind1d[0], #indices_nd=indnd,
                         start_time=self.start_time, end_time=self.end_time)

    def near_neighbor_tunnel(self, use_time=False, verbose=False, timeit=True):
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
            True to limit results to closest time value. Defaults False.
        verbose: bool
            True returns screen print out of pair results. Defaults False.
        timeit: bool
            True prints elapsed time of operation.
        '''
        opbegin = timer.time()
        # Create lists to contain the indices for each flight point
        ind1d = []
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
            ind1d.append(minindex_1d)
        indnd = np.unravel_index(ind1d, self.latvalsr.shape)

        if verbose:
            self._print_pair_by_index(ind1d)

        matchdata = self._get_data_by_index(ind1d, )
        if timeit:
            print("--- Elapsed time: %s sec ---" % (timer.time() - opbegin))
        return MatchData(self.fldata, matchdata, distance_to_point=distance,
                         indices_1d=ind1d[0], indices_nd=indnd,
                         start_time=self.start_time, end_time=self.end_time)

########################
#  Get/print methods   #
########################

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

    def _get_data_by_index(self, indices_1d, orig_indices):
        '''Retrieve data fields using calculated indices.'''
        data_by_ind = {}
        for field in self.data_fields.keys():
            dfield = self.data_fields[field]
            indices_nd = np.unravel_index(orig_indices[indices_1d], dfield['data'].shape)
            data_by_ind[field] = dfield.copy()
            data_by_ind[field]['data'] = np.ravel(dfield['data'])[indices_1d]
            data_by_ind[field]['original_indices'] = indices_nd
        return data_by_ind

    def _print_pair_by_index(self, ind1d):
        for ii, ind in enumerate(ind1d):
            print(("AC Lat: %g, Lon: %g, Alt: %g | "
                   "Rad Lat: %g, Lon: %g, Alt: %g") % (
                  self.flight_lon['data'][ii],
                  self.flight_lat['data'][ii],
                  np.ravel(self.data_lon)[ind],
                  np.ravel(self.data_lat)[ind],
                  ))
            print("AC/Rad Turb: %g / %g" % (
                self.fldata['turb']['data'][ii],
                np.ravel(self.data_fields['turbulence']['data'])[ind]))


class FlightLevelMatch(object):
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
        # Set variables to be used later by any function
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
        # we need to convert to AWOT common epoch times
        flight_epoch = common.convert_to_epoch_dict(flight[timekey])
        self.flight_time = self._get_var_time_subset(
                             flight_epoch,
                             start_time, end_time)
        self.flight_numtime = date2num(
            self.flight_time['data'], self.flight_time['units'])

        # Create a field to store flight-matched data
        self.fldata = {}
        self.fldata['longitude'] = self.flight_lon
        self.fldata['latitude'] = self.flight_lat
        self.fldata['altitude'] = self.flight_alt
        self.fldata['time'] = self.flight_time

        for field in flight.keys():
            try:
                self.fldata[field] = self._get_var_time_subset(
                    flight[field].copy(), start_time, end_time)
            except:
                self.fldata[field] = None

        # Create a field to store data field to matching data
        self.data_fields = {}
        if field_match_dict is None:
            for field in self.datavol['fields'].keys():
                self.data_fields[field] = self.datavol['fields'][field].copy()
        else:
            for field in field_match_dict:
                self.data_fields[field] = self.datavol['fields'][field].copy()
        self.data_shape = self.data_fields[
            self.data_fields.keys()[0]]['data'].shape

        # Run checks to make sure array lengths are the same
        warntxt = "Check that data is an AWOT style object or specify keyword!"
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
                self.data_time = common.convert_to_epoch_dict(
                    self.datavol['time'])
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
        self.data_numtime = date2num(
            self.data_time['data'][:], self.time['units'])

        # Convert time to same shape as fields if not already
        if len(self.data_numtime.shape) == 1:
            self.data_numtime = np.resize(self.data_numtime, self.data_shape)

        self.start_time = start_time
        self.end_time = end_time

        # Convert latitude and longitude to radians
        self.latvalsr = np.radians(self.data_lat)
        self.lonvalsr = np.radians(self.data_lon)
        self.lat0valsr = np.radians(self.flight_lat['data'][:])
        self.lon0valsr = np.radians(self.flight_lon['data'][:])

    def kdtree(self, leafsize=16,
               use_time=False, verbose=False, timeit=True,
               query_k=1, query_eps=0, query_p=2,
               query_distance_upper_bound=np.inf,
               query_n_jobs=1):
        '''
        Find the closest point using a K-Dimensional Tree method.
        The tree building can take a long time depending on how
        large the input data array. Subsetting is attemped to increase
        execution speed.

        Parameters
        ----------
        leafsize : int
            Positive number of points at which time
            scipy.spatial.cKDTree switches to brute force.
            Default same as scipy default.
        use_time: bool
            True includes the time variables in KD-Tree. Defaults False.
        verbose: bool
            True returns screen print out of pair results. Defaults False.
        timeit: bool
            True prints elapsed time of operation.
        query_?? : See scipy.spatial.cKDTree.query
            Default keyword values used in the cKDTree.query.
        '''
        opbegin = timer.time()
        # Create 1-D variables to use
        dlon, dlat = np.ravel(self.data_lon), np.ravel(self.data_lat)
        dht, dti = np.ravel(self.data_height), np.ravel(self.data_numtime)

        # Subset input data to only include relevant points within
        condlon = np.logical_and(dlon >= np.min(self.flight_lon['data']),
                                 dlon <= np.max(self.flight_lon['data']))
        condlat = np.logical_and(dlat >= np.min(self.flight_lat['data']),
                                 dlat <= np.max(self.flight_lat['data']))
        condht = np.logical_and(dht >= np.min(self.flight_alt['data']),
                                dht <= np.max(self.flight_alt['data']))
        condti = np.logical_and(dti >= np.min(self.flight_numtime),
                                dti <= np.max(self.flight_numtime))

        cond = ((dlon >= np.min(self.flight_lon['data'])) &
                (dlon <= np.max(self.flight_lon['data'])) &
                (dlat >= np.min(self.flight_lat['data'])) &
                (dlat <= np.max(self.flight_lat['data'])) &
                (dht >= np.min(self.flight_alt['data'])) &
                (dht <= np.max(self.flight_alt['data'])) &
                (dti >= np.min(self.flight_numtime)) &
                (dti <= np.max(self.flight_numtime)))
        indices = np.where(cond)
        print(len(indices[0]))

        if len(indices[0]) > 1:
            dlon = dlon[indices[0]]
            dlat = dlat[indices[0]]
            dht = dht[indices[0]]
            dti = dti[indices[0]]
        else:
            print("No points found for condition")

        # Build and query kd-tree
        if use_time:
            kdt = cKDTree(zip(dlon, dlat, dht, dti), leafsize=leafsize)
        else:
            kdt = cKDTree(zip(dlon, dlat, dht), leafsize=leafsize)
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
        indnd = np.unravel_index(ind1d, self.latvalsr.shape)

        if verbose:
            self._print_pair_by_index(ind1d)

        matchdata = self._get_data_by_index(ind1d)
        if timeit:
            print("--- Elapsed time: %s sec ---" % (timer.time() - opbegin))
        return MatchData(self.fldata, matchdata, distance_to_point=distance,
                         indices_1d=ind1d[0], indices_nd=indnd,
                         start_time=self.start_time, end_time=self.end_time)

    def near_neighbor_tunnel(self, use_time=False, verbose=False, timeit=True):
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
            True to limit results to closest time value. Defaults False.
        verbose: bool
            True returns screen print out of pair results. Defaults False.
        timeit: bool
            True prints elapsed time of operation.
        '''
        opbegin = timer.time()
        # Create lists to contain the indices for each flight point
        ind1d = []
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
            ind1d.append(minindex_1d)
        indnd = np.unravel_index(ind1d, self.latvalsr.shape)

        if verbose:
            self._print_pair_by_index(ind1d)

        matchdata = self._get_data_by_index(ind1d)
        if timeit:
            print("--- Elapsed time: %s sec ---" % (timer.time() - opbegin))
        return MatchData(self.fldata, matchdata, distance_to_point=distance,
                         indices_1d=ind1d[0], indices_nd=indnd,
                         start_time=self.start_time, end_time=self.end_time)

########################
#  Get/print methods   #
########################

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

    def _get_data_by_index(self, indices_1d):
        '''Retrieve data fields using calculated indices.'''
        data_by_ind = {}
        for field in self.data_fields.keys():
            dfield = self.data_fields[field]
            data_by_ind[field] = dfield.copy()
            data_by_ind[field]['data'] = np.ravel(dfield['data'])[indices_1d]
        return data_by_ind

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
            print("AC/Rad Turb: %g / %g" % (
                self.fldata['turb']['data'][ii],
                np.ravel(self.data_fields['turbulence']['data'])[ind]))


class RadarMatch(object):
    """Class for matching flight level data to Py-ART radar objects."""

    def __init__(self, flight, radar, volume_update_time=0,
                 start_time=None, end_time=None,
                 mask_below=None, mask_above=None,
                 flight_lon_name=None, flight_lat_name=None,
                 flight_alt_name=None, flight_time_name=None,
                 flight_uwind_name=None, flight_vwind_name=None,
                 field_match_dict=None,
                 ):
        '''
        Intitialize the class with provided data.

        Parameters
        ----------
        flight : dict
            AWOT flight data dictionary object.
        radar : object or list of objects
            A Py-ART radar object dervied using Py-ART read statement.
            If a list is of Py-ART radar objects is supplied the function
            loops through each instance looking for nearest neighbor points.

        Optional
        --------
        volume_update_time : int
            Update time [in seconds] for volume scan. If set half of the
            time will be subtracted from from the beginning and end, thus
            treating the volume start time as the mid-point for comparison
            to aircraft track.
        start_time : str
            UTC time to use as start time for subsetting in datetime format.
            (e.g. 2014-08-20 12:30:00)
        end_time : str
            UTC time to use as an end time for subsetting in datetime format.
            (e.g. 2014-08-20 16:30:00)
        mask_below : float
            Values of matched data less than mask_below will be masked.
        mask_above : float
            Values of matched data greater than mask_above will be masked.
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
        # Set data to be used later by any function
        self.time = flight['time']
        self.volup = volume_update_time
        self.masklo = mask_below
        self.maskhi = mask_above

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
        # we need to convert to AWOT common epoch times
        flight_epoch = common.convert_to_epoch_dict(flight[timekey])
        self.flight_time = self._get_var_time_subset(
                             flight_epoch,
                             start_time, end_time)
        self.flight_numtime = date2num(
            self.flight_time['data'], self.flight_time['units'])

        # Create a field to store flight-matched data
        self.fldata = {}
        self.fldata['longitude'] = self.flight_lon
        self.fldata['latitude'] = self.flight_lat
        self.fldata['altitude'] = self.flight_alt
        self.fldata['time'] = self.flight_time

        for field in flight.keys():
            try:
                self.fldata[field] = self._get_var_time_subset(
                    flight[field].copy(), start_time, end_time)
            except:
                self.fldata[field] = None

        self.start_time = start_time
        self.end_time = end_time

        # Convert latitude and longitude to radians
        self.lat0valsr = np.radians(self.flight_lat['data'][:])
        self.lon0valsr = np.radians(self.flight_lon['data'][:])

        # Check if list is passed, if single build a list of 1
        if type(radar).__name__ == 'list':
            self.rlist = radar
        else:
            self.rlist = [radar]

    def kdtree(self, leafsize=16,
               verbose=False, timeit=True,
               query_k=1, query_eps=0, query_p=2,
               query_distance_upper_bound=np.inf,
               query_n_jobs=1):
        '''
        Find the closest point using a K-Dimensional Tree method.
        The tree building can take a long time depending on how
        large the input data array. Subsetting is attemped to increase
        execution speed.

        Parameters
        ----------
        leafsize : int
            Positive number of points at which time
            scipy.spatial.cKDTree switches to brute force.
            Default same as scipy default.
        verbose: bool
            True returns screen print out of pair results. Defaults False.
        timeit: bool
            True prints elapsed time of operation.
        query_?? : See scipy.spatial.cKDTree.query
            Default keyword values used in the cKDTree.query.
        '''
        opbegin = timer.time()
        # Create processing variables
        ind1d = []
        distance = []
        prdata = {}
        # We assume that all files have the same data field dimensions
        for field in self.rlist[0].fields.keys():
            prdata[field] = self.rlist[0].fields[field].copy()
            prdata[field]['data'] = np.full_like(self.flight_numtime, np.nan)

        # Loop through the list or radar instances
        for num, pr in enumerate(self.rlist):
            dfield = {}
            # Create 1-D variables to use
            dlon = np.ravel(pr.gate_longitude['data'].copy())
            dlat = np.ravel(pr.gate_latitude['data'].copy())
            dht = np.ravel(pr.gate_altitude['data'].copy())
            for field in pr.fields.keys():
                dfield[field] = np.ravel(pr.fields[field]['data'])
            dshape = pr.fields[field]['data'].shape

            # Convert the radar time to AWOT epoch
            rtime = common.convert_to_epoch_dict(pr.time.copy())
            rtime['data'] = (
                rtime['data'] -
                datetime.timedelta(seconds=self.volup/2.))
            rnumtime = date2num(rtime['data'], rtime['units'])
            if len(rnumtime.shape) == 1:
                rnumtime = np.resize(rnumtime, dshape)
            dti = np.ravel(rnumtime)

            # Subset input data to relevant points
            cond = ((dlon >= np.min(self.flight_lon['data'])) &
                    (dlon <= np.max(self.flight_lon['data'])) &
                    (dlat >= np.min(self.flight_lat['data'])) &
                    (dlat <= np.max(self.flight_lat['data'])) &
                    (dht >= np.min(self.flight_alt['data'])) &
                    (dht <= np.max(self.flight_alt['data'])) &
                    (dti >= np.min(self.flight_numtime)) &
                    (dti <= np.max(self.flight_numtime)))
            indices = np.where(cond)

            if len(indices[0]) > 1:
                dlon = dlon[indices[0]]
                dlat = dlat[indices[0]]
                dht = dht[indices[0]]
                dti = dti[indices[0]]
                for field in pr.fields.keys():
                    dfield[field] = dfield[field][indices[0]]
            else:
                print("No points found for condition on file %d" % num)

            # Subset for matching times
            tcond = np.logical_and(self.flight_numtime >= np.min(dti),
                                   self.flight_numtime <= np.max(dti))
            indt = np.where(tcond)

            # Build and query kd-tree
            kdt = cKDTree(zip(dlon, dlat, dht, dti), leafsize=leafsize)
            if np.size(indt[0]) > 0:
                prdistance, prind1d = kdt.query(
                    zip(self.flight_lon['data'][indt[0]],
                        self.flight_lat['data'][indt[0]],
                        self.flight_alt['data'][indt[0]],
                        self.flight_numtime[indt[0]]),
                    k=query_k, eps=query_eps, p=query_p,
                    distance_upper_bound=query_distance_upper_bound,
                    n_jobs=query_n_jobs)
                distance.append(prdistance)
                ind1d.append(prind1d)

                for field in pr.fields.keys():
                    prdata[field]['data'][indt[0]] = dfield[field][prind1d]

                if verbose:
                    for ii, ind in enumerate(ind1d[0]):
                        print(("AC Lat: %g, Lon: %g, Alt: %g | "
                               "Rad Lat: %g, Lon: %g, Alt: %g") % (
                              self.flight_lon['data'][ii],
                              self.flight_lat['data'][ii],
                              self.flight_alt['data'][ii],
                              dlon[ind],
                              dlat[ind],
                              dht[ind]))

        # Apply masks to data
        for field in prdata.keys():
            prdata[field]['data'] = np.ma.masked_invalid(prdata[field]['data'])
            if self.masklo is not None:
                prdata[field]['data'] = np.ma.masked_less(
                    prdata[field]['data'], self.masklo)
            if self.maskhi is not None:
                prdata[field]['data'] = np.ma.masked_greater(
                    prdata[field]['data'], self.maskhi)
        if timeit:
            print("--- Elapsed time: %s sec ---" % (timer.time() - opbegin))
        return MatchData(self.fldata, prdata,
                         distance_to_point=distance, indices_1d=ind1d[0],
                         start_time=self.start_time, end_time=self.end_time)

    def near_neighbor_tunnel(self, verbose=False, timeit=True):
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
        verbose: bool
            True returns screen print out of pair results. Defaults False.
        timeit: bool
            True prints elapsed time of operation.
        '''
        opbegin = timer.time()
        # Create processing variables
        ind1d = []
        distance = []
        prdata = {}
        # We assume that all files have the same data field dimensions
        for field in self.rlist[0].fields.keys():
            prdata[field] = self.rlist[0].fields[field].copy()
            prdata[field]['data'] = np.full_like(self.flight_numtime, np.nan)

        # Calculate sines and cosines of flight lats/lon
        clats0, clons0 = np.cos(self.lat0valsr), np.cos(self.lon0valsr)
        slons0, slats0 = np.sin(self.lon0valsr), np.sin(self.lat0valsr)

        # Loop through radar instances
        for num, pr in enumerate(self.rlist):
            dfield = {}
            prind1d = []
            dlon = np.ravel(pr.gate_longitude['data'].copy())
            dlat = np.ravel(pr.gate_latitude['data'].copy())
            dht = np.ravel(pr.gate_altitude['data'].copy())
            for field in pr.fields.keys():
                dfield[field] = np.ravel(pr.fields[field]['data'])
            dshape = pr.fields[field]['data'].shape

            # Convert the radar time to AWOT epoch
            rtime = common.convert_to_epoch_dict(pr.time.copy())
            rtime['data'] = (
                rtime['data'] -
                datetime.timedelta(seconds=self.volup/2.))
            rnumtime = date2num(rtime['data'], rtime['units'])
            if len(rnumtime.shape) == 1:
                rnumtime = np.resize(rnumtime, dshape)
            dti = np.ravel(rnumtime)

            # Subset input data to relevant points
            cond = ((dlon >= np.min(self.flight_lon['data'])) &
                    (dlon <= np.max(self.flight_lon['data'])) &
                    (dlat >= np.min(self.flight_lat['data'])) &
                    (dlat <= np.max(self.flight_lat['data'])) &
                    (dht >= np.min(self.flight_alt['data'])) &
                    (dht <= np.max(self.flight_alt['data'])) &
                    (dti >= np.min(self.flight_numtime)) &
                    (dti <= np.max(self.flight_numtime)))
            indices = np.where(cond)

            if len(indices[0]) > 1:
                dlon = dlon[indices[0]]
                dlat = dlat[indices[0]]
                dht = dht[indices[0]]
                dti = dti[indices[0]]
                for field in pr.fields.keys():
                    dfield[field] = dfield[field][indices[0]]
            else:
                print("No points found for condition on file %d" % num)

            clat, clon = np.cos(np.radians(dlat)), np.cos(np.radians(dlon))
            slon, slat = np.sin(np.radians(dlat)), np.sin(np.radians(dlon))

            # Subset for matching times
            tcond = np.logical_and(self.flight_numtime >= np.min(dti),
                                   self.flight_numtime <= np.max(dti))
            indt = np.where(tcond)

            if np.size(indt[0]) > 0:
                for ii in range(len(self.flight_numtime)):
                    delX = clats0[ii] * clons0[ii] - clat * clon
                    delY = clats0[ii] * slons0[ii] - clat * slon
                    delZ = slats0[ii] - slat
                    # Calculate the distance squared
                    dist_sq = delX**2 + delY**2 + delZ**2
                    # Find the 1D index of the minimum element
#                    prind1d.append(dist_sq.argmin())
#                    distance.append(dist_sq.min())
                    prind1d.append(np.ma.argmin(dist_sq))
                    distance.append(np.ma.min(dist_sq))
                ind1d.append(prind1d)

                for field in pr.fields.keys():
                    prdata[field]['data'][indt[0]] = dfield[field][prind1d]

                if verbose:
                    for ii, ind in enumerate(ind1d[0]):
                        print(("AC Lat: %g, Lon: %g, Alt: %g | "
                               "Rad Lat: %g, Lon: %g, Alt: %g") % (
                              self.flight_lon['data'][ii],
                              self.flight_lat['data'][ii],
                              self.flight_alt['data'][ii],
                              dlon[ind],
                              dlat[ind],
                              dht[ind]))

        # Apply masks to data
        for field in prdata.keys():
            prdata[field]['data'] = np.ma.masked_invalid(prdata[field]['data'])
            if self.masklo is not None:
                prdata[field]['data'] = np.ma.masked_less(
                    prdata[field]['data'], self.masklo)
            if self.maskhi is not None:
                prdata[field]['data'] = np.ma.masked_greater(
                    prdata[field]['data'], self.maskhi)
        if timeit:
            print("--- Elapsed time: %s sec ---" % (timer.time() - opbegin))
        return MatchData(self.fldata, prdata,
                         distance_to_point=distance, indices_1d=ind1d[0],
                         start_time=self.start_time, end_time=self.end_time)

    def near_neighbor(self, verbose=False, timeit=True):
        """
        Calculate variables using a nearest neighbor approach when the input
        data is a Py-ART radar object.

        Parameters
        ----------
        verbose: bool
            True returns screen print out of pair results. Defaults False.
        timeit: bool
            True prints elapsed time of operation.
        """
        opbegin = timer.time()
        # Create processing variables
        indaz = []
        indrng = []
        indel = []
        prdata = {}
        # We assume that all files have the same data field dimensions
        for field in self.rlist[0].fields.keys():
            prdata[field] = self.rlist[0].fields[field].copy()
            prdata[field]['data'] = np.full_like(self.flight_numtime, np.nan)

        # Loop through radar instances
        for num, pr in enumerate(self.rlist):
            rlat, rlon = pr.latitude['data'][0], pr.longitude['data'][0]

            # Calculate the distance to the aircraft from radar
            ac_dist_hav = self.ground_distance_haversine(
                rlat, rlon, self.flight_lat['data'], self.flight_lon['data'])
            ac_elev2 = np.degrees(np.arctan2(self.flight_alt['data'][:],
                                            ac_dist_hav))
            # Calculate the elevation angle corresponding to aircraft position
#            ac_rng2, ac_elev2 = self.slant_range_and_elev(
#                ac_dist_hav, self.flight_alt['data'][:])

            # Calculate the azimuth angle corresponding to aircraft position
            ac_az2 = self.azimuth_to_point(
                rlat, rlon, self.flight_lat['data'],
                self.flight_lon['data'], ac_dist_hav)

            ac_dist_cos = self.ground_distance_cosines(
                rlat, rlon, self.flight_lat['data'], self.flight_lon['data'])
            # Calculate the azimuth angle corresponding to aircraft position
            ac_az = self.bearing_to_point(
                rlat, rlon, self.flight_lat['data'],
                self.flight_lon['data'])
            # Calculate the elevation angle corresponding to aircraft position
            ac_rng, ac_elev = self.slant_range_and_elev(
                ac_dist_cos, self.flight_alt['data'][:])

            # Convert the radar time to AWOT epoch
            rtime = common.convert_to_epoch_dict(pr.time.copy())
            rtime['data'] = (
                rtime['data'] -
                datetime.timedelta(seconds=self.volup/2.))
            rnumtime = date2num(rtime['data'], rtime['units'])

            # Subset for matching times
            tcond = np.logical_and(self.flight_numtime >= rnumtime[0],
                                   self.flight_numtime <= rnumtime[-1])
            indt = np.where(tcond)

            if np.size(indt[0]) > 0:
                for index in indt[0]:
                    elin = np.argmin(
                        np.abs(ac_elev[index] - pr.fixed_angle['data']))
                    pr_sweep = pr.extract_sweeps([elin])
                    azin = np.argmin(
                        np.abs(ac_az[index] - pr_sweep.azimuth['data']))
                    rgin = np.argmin(
                        np.abs(ac_rng[index] - pr_sweep.range['data']))

                    elin2 = np.argmin(
                        np.abs(ac_elev2[index] - pr.fixed_angle['data']))
                    pr_sweep2 = pr.extract_sweeps([elin2])
                    azin2 = np.argmin(
                        np.abs(ac_az2[index] - pr_sweep2.azimuth['data']))
                    rgin2 = np.argmin(
                        np.abs(ac_dist_hav[index] - pr_sweep2.range['data']))

                    if verbose:
                        print("AC/Radar lat/lon/alt: %g/%g, %g/%g, %g/%g" % (
                            self.flight_lat['data'][index],
                            pr.gate_latitude['data'][azin, rgin],
                            self.flight_lon['data'][index],
                            pr.gate_longitude['data'][azin, rgin],
                            self.flight_alt['data'][index],
                            pr.gate_altitude['data'][azin, rgin]))

                    for field in pr_sweep.fields.keys():
                        prdata[field]['data'][index] = pr_sweep.fields[field]['data'][azin, rgin]

        # Apply masks to data
        for field in prdata.keys():
            prdata[field]['data'] = np.ma.masked_invalid(prdata[field]['data'])
            if self.masklo is not None:
                prdata[field]['data'] = np.ma.masked_less(
                    prdata[field]['data'], self.masklo)
            if self.maskhi is not None:
                prdata[field]['data'] = np.ma.masked_greater(
                    prdata[field]['data'], self.maskhi)
        if timeit:
            print("--- Elapsed time: %s sec ---" % (timer.time() - opbegin))
        return MatchData(self.fldata, prdata,
                         start_time=self.start_time, end_time=self.end_time)

#################################
#  RadarMatch Helper functions  #
#################################

    def ground_distance_haversine(self, lat0, lon0, aclats, aclons, R=6371100.):
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

    def ground_distance_cosines(self, lat0, lon0, aclats, aclons, R=6371100.):
        '''
        Calculate the distance to the aircraft using the spherical law
        of cosines.

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
        top = (np.sin(np.radians(aclats)) -
               np.cos(rng/factor) * np.sin(np.radians(lat0)))
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
        Method was adopted from PyTDA by Timothy Lang from RSL library.
        Great circle formulation.

        Parameters
        ----------
        lat0 : float
            Latitude [degrees] value of origin point from which to calculate.
        lon0 : float
            Longitude [degrees] value of origin point from which to calculate.
        aclats : float array
            Array of aircraft latitude [degrees] values.
        aclons : float array
            Array of aircraft longitude [degrees] values.
        '''
        aclonsr = np.radians(aclons)
        aclatsr = np.radians(aclats)
        lon0r = np.radians(lon0)
        lat0r = np.radians(lat0)

        bear = np.arctan2(
            (np.sin(aclonsr - lon0r) * np.cos(aclatsr)),
            (np.cos(lat0r) * np.sin(aclatsr) -
             np.sin(lat0r) * np.cos(aclatsr) * np.cos(aclonsr - lon0r)))
        return np.degrees(bear) % 360.

    def slant_range_and_elev(self, ground_range, height, R=6371100.):
        '''
        Calculate slant range and elevation.
        Method was adopted from PyTDA by Timothy Lang from RSL library.

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
                      "(e.g. '2014-08-20T12:30:00')")
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
    def __init__(self, flight, data, distance_to_point=None,
                 indices_1d=None, indices_nd=None,
                 start_time=None, end_time=None):
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
        self.start_time = start_time
        self.end_time = end_time
        return
