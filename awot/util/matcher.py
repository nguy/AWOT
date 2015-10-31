"""
awot.util.matcher
================

Matching function for airborne and ground-based radar data.

Original version provided by George Duffy 'APRmacher2.py'.
This code looks to generalize it to a formation that works with
AWOT objects.

"""
import numpy as np
import datetime
import calendar
import scipy.io
import scipy.spatial
import pdb
from operator import itemgetter
from mpl_toolkits.basemap import Basemap

class TrackMatch(object):
    """Class for matching flight level data to other observations."""

    def __init__(self, flight, data, basemap=None,
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
        data : dict
            AWOT dictionary object.

        Optional
        --------
        start_time : str
            UTC time to use as start time for subsetting in datetime format.
            (e.g. 2014-08-20 12:30:00)
        end_time : str
            UTC time to use as an end time for subsetting in datetime format.
            (e.g. 2014-08-20 16:30:00)
        basemap : basemap instance
            Basemap instance to use for plotting.
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
        self.flight_data = flight
        self.data = data
        self.time = flight['time']

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

        self.flight_lon = self._get_var_time_subset(
                             flight[lonkey].copy(),
                             start_time, end_time)
        self.flight_lat = self._get_var_time_subset(
                             flight[latkey].copy(),
                             start_time, end_time)
        self.flight_alt = self._get_var_time_subset(
                             flight[altkey].copy(),
                             start_time, end_time)
        self.flight_Uwind = self._get_var_time_subset(
                             flight[ukey].copy(),
                             start_time, end_time)
        self.flight_Vwind = self._get_var_time_subset(
                             flight[vkey].copy(),
                             start_time, end_time)
        self.flight_time = self._get_var_time_subset(
                             flight[timekey].copy(),
                             start_time, end_time)

        # Create a field to store matched data
        self.matchdata = {}
        self.matchdata['flight_longitude'] = self.flight_lon
        self.matchdata['flight_latitude'] = self.flight_lat
        self.matchdata['flight_altitude'] = self.flight_alt
        self.matchdata['flight_time'] = self.flight_time

        for field in flight.keys():
            try:
                self.matchdata[field] = self._get_var_time_subset(
                                          flight[field].copy(),
                                          start_time, end_time)
            except:
                self.matchdata[field] = None

        self.data_fields = {}
        if field_match_dict is None:
            for field in self.data['fields'].keys():
                self.data_fields[field] = self.data['fields'][field]
        else:
            for field in field_match_dict:
                self.data_fields[field] = self.data['fields'][field]

        # Run checks to make sure array lengths are the same
        if data_lon is None:
            try:
                self.data_lon = self.data['longitude']['data']
            except:
                print("Check that data is an AWOT object!")
        else:
            self.data_lon = np.array(data_lon)

        if data_lat is None:
            try:
                self.data_lat = self.data['latitude']['data']
            except:
                print("Check that data is an AWOT object!")
        else:
            self.data_lat = np.array(data_lat)

        if data_height is None:
            try:
                self.data_height = self.data['height']['data']
            except:
                print("Check that data is an AWOT object!")
        else:
            self.data_height = np.array(data_height)

        if data_time is None:
            try:
                self.data_time = self.data['time']
            except:
                print("Check that data is an AWOT object!")
        else:
            self.data_time = data_time

        self.basemap = basemap

        # Calculate x,y map position coordinates
        if self.basemap is not None:
            self.flight_x, self.flight_y = self.basemap(
                 self.flight_lon['data'][:].ravel(),
                 self.flight_lat['data'][:].ravel())

            self.data_x, self.data_y = self.basemap(self.data_lon, self.data_lat)

        # Convert latitude and longitude to radians
        self.latvals = np.radians(self.data_lat)
        self.lonvals = np.radians(self.data_lon)
        self.lat0vals = np.radians(self.flight_lat['data'][:])
        self.lon0vals = np.radians(self.flight_lon['data'][:])

    def kdtree(self, leafsize=16,
               print_match_pairs=False):
        '''
        Find the closest point using a KD Tree method.

        Parameters
        ----------
        leafsize : int
            Positive number of points at which time
            scipy.spatial.cKDTree switches to brute force.
            Default same as scipy default.
        print_match_pairs: bool
            True returns screen print out of pair results.
            Default is False.
        '''
        # Create lists to contain the indices for each flight point
        indlat = []
        indlon = []

        # Calculate sines and cosines of lats/lon
        clat, clon = np.cos(self.latvals), np.cos(self.lonvals)
        slon, slat = np.sin(self.lonvals), np.sin(self.latvals)
        clats0, clons0 = np.cos(self.lat0vals), np.cos(self.lon0vals)
        slons0, slats0 = np.sin(self.lon0vals), np.sin(self.lat0vals)

        # Build kd-tree from big arrays of 3D coordinates
        triples = list(zip(np.ravel(self.data_lon), np.ravel(self.data_lat), np.ravel(self.data_height)))
        kdt = scipy.spatial.cKDTree(triples)

        for ii in range(len(self.lat0vals)):
            dist_sq_min, minindex_1d = kdt.query([self.flight_lon['data'][ii],
                                                  self.flight_lat['data'][ii],
                                                  self.flight_alt['data'][ii]])
            iy_min, ix_min = np.unravel_index(minindex_1d, self.data_lat.shape)
            indlon.append(iy_min)
            indlat.append(ix_min)

        if print_match_pairs:
            self._print_pair_by_index(indlon, indlat)

        self._get_matchdata_by_index(indlon, indlat)
        return self.matchdata

    def near_neighbor_tunnel(self, start_time=None, end_time=None):
        '''
        Find closest point to the set of (lat,lon) points
        provided by the flight data object to data object points.
        This routine was adapted from:
        http://nbviewer.ipython.org/github/Unidata/tds-python-workshop/blob/master/netcdf-by-coordinates.ipynb

        Returns iy,ix such that the square of the tunnel distance
        between (latval[it,ix],lonval[iy,ix]) and (lat0,lon0)
        is minimum.
        '''
        # Create lists to contain the indices for each flight point
        indlat = []
        indlon = []

        # Calculate sines and cosines of lats/lon
        clat, clon = np.cos(self.latvals), np.cos(self.lonvals)
        slon, slat = np.sin(self.lonvals), np.sin(self.latvals)
        clats0, clons0 = np.cos(self.lat0vals), np.cos(self.lon0vals)
        slons0, slats0 = np.sin(self.lon0vals), np.sin(self.lat0vals)

        for ii in range(len(self.lat0vals)):
            delX = clats0[ii] * clons0[ii] - clat * clon
            delY = clats0[ii] * slons0[ii] - clat * slon
            delZ = slats0[ii] - slat
            # Calculate the distance squared
            dist_sq = delX**2 + delY**2 + delZ**2
            # Find the 1D index of the minimum element
            minindex_1d = dist_sq.argmin()
            iy_min, ix_min = np.unravel_index(minindex_1d, self.latvals.shape)
            indlon.append(iy_min)
            indlat.append(ix_min)

        self._get_matchdata_by_index(indlon, indlat)
        return self.matchdata

#### CODE SNIPPET FROM GEORGE DUFFY - EQUIVALENT TO KDTREE ABOVE? ###
#### DIFFERENCE IS THAT TIME IS ALSO SEARCHED IN THESE ROUTINES ###
#     def weighted_stats(threshold_time=180, threshold_hor=250,
#                      threshold_vert=125, threshold_ac=6):
#         """
#         Match the radar field values to each aicraft time
#         """
#         dt, dr, timestamp, actime, height, longheight, matchtimes = [], [], [], [], [], [], []
#         for i in range(len(self.flight_time['data'][:])):
#             timespan = range(int(self.flight_time['data'][i] - threshold_time),
#                                  int(self.flight_time['data'][i] + threshold_time))
#             rdex = np.in1d(self.data_time['data'][:], timespan)
#             if np.mod(i, 1000) == 0:
#                 print(i)
#             if np.any(self.data_time['data'][:][rdex]):
#                 field = np.dstack([self.data_x[rdex, :, :].ravel(),
#                                    self.data_y[rdex, :, :].ravel(),
#                                    self.data_height[rdex, :, :].ravel()])[0]
#                 tree = scipy.spatial.cKDTree(field)
#                 goahead = True
#                 predist, predex = tree.query([self.flight_x[i],
#                                               self.flight_y[i],
#                                               self.flight_alt['data'][i]], 1)
#                 if np.any(predex):
# #***                    loctime = itemgetter(predex)(radtime_3D[rdex,:,:].ravel())
#                     dtime = self.flight_time['data'][i] - loctime
#                     trux = self.flight_x[i] - self.flight_Uwind['data'][i] * dtime
#                     truy = self.flight_y[i] - self.flight_Vwind['data'][i] * dtime
#                     truz = self.flight_alt['data'][i]
#                     closest, matchdexa = tree.query([trux[0], truy[0], truz], 10)
#                     if np.max(closest) > 300:
#                         goahead = False
#                     else:
#                         matchtimes = np.hstack((matchtimes, self.flight_time['data'][:][i]*np.ones(np.shape(matchdexa))))
#                     ###Retrieve reflectivities with weighted statistics
#                     if goahead:
#                         matchdex = matchdexa.astype(int)
#                         branch14 = itemgetter(matchdex)(flec14[rdex, :, :].ravel())
#                         branch35 = itemgetter(matchdex)(flec35[rdex, :, :].ravel())
#                         timebranch = itemgetter(matchdex)(radtime_3D[rdex, :, :].ravel())
#                         if np.all((branch14 > -100) & (branch35 > -100)):
#                             if np.any(branch14 > 30):
#                                 branch14 = branch14[branch14 < 30]
#                                 branch35 = branch14[branch14 < 30]
#                             print('match')
#                             timestamp.append(np.mean(timebranch))
#                             actime.append(self.flight_time['data'][i])
#                             dt.append(np.mean(np.abs(matchtimes - timestamp[-1])))
#                             dt_std.append(np.std(np.abs(matchtimes - timestamp[-1])))
#         return match_data

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

        # First let's save the fields data into this class structure
        for field in radar.fields.keys():
                self.matchdata[field] = radar.fields[field].copy()
                self.matchdata[field]['data'] = np.ma.empty(len(self.lat0vals))

        rlat, rlon = radar.latitude['data'][0], radar.longitude['data'][0]

        # Calculate the distance to the aircraft from radar
        dist_to_ac = self.calc_dist_to_aircraft(rlat, rlon,
                                           self.flight_lat['data'],
                                           self.flight_lon['data'])

        # Calculate the elevation angle corresponding to aircraft position
        ac_elev = np.degrees(np.arctan2(self.flight_alt['data'][:], dist_to_ac))

        ac_az = self.calc_az_to_aircraft(rlat, rlon,
                                           self.flight_lat['data'],
                                           self.flight_lon['data'], dist_to_ac)
        for ii in range(len(self.lat0vals)):
            elin = np.argmin(np.abs(ac_elev[ii] - radar.fixed_angle['data']))
            radar_sweep = radar.extract_sweeps([elin])
            azin = np.argmin(np.abs(ac_az[ii]-radar_sweep.azimuth['data']))
            rgin = np.argmin(np.abs(dist_to_ac[ii] - radar_sweep.range['data']/1000.0))
            for field in self.data_fields.keys():
                self.data_fields[field]['data'][ii] = radar_sweep.fields[field]['data'][azin, rgin]

        return self.matchdata
#            tpos1 = sweep_latlon_to_flat_xy(radar, elin)
#            ts_dz[index] = radar_sweep.fields['reflectivity']['data'][azin, rgin]
#            ts_edr[index] = radar_sweep.fields['turbulence']['data'][azin, rgin]
#            ts_vel[index] = radar_sweep.fields['velocity']['data'][azin, rgin]
#            ts_sw[index] = radar_sweep.fields['spectrum_width']['data'][azin, rgin]

###################
#  Calculations   #
###################

    def calc_dist_to_aircraft(self, lat0, lon0, aclats, aclons):
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
        '''
        R = 6371000. # Earth radius
        delLat = lat0 - aclats
        delLon = lon0 - aclons

        a = (np.sin(np.radians(delLat)/2))**2 + np.cos(np.radians(lat0)) * np.cos(np.radians(aclats)) * (np.sin(np.radians(delLon)/2))**2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
        distance = R * c
        return distance

    def calc_az_to_aircraft(self, lat0, lon0, aclats, aclons, rng):
        R = 6371000. # Earth radius
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
            print("AC Lat: %g, Lon: %g, Alt: %g | Rad Lat: %g, Lon: %g, Alt: %g"%(
                  self.flight_lon['data'][ii], self.flight_lat['data'][ii],
                  self.flight_alt['data'][ii], self.data_lon[iy_min, ix_min],
                  self.data_lat[iy_min, ix_min], self.data_height[iy_min, ix_min]))

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
        '''Get a subsetted time and Variable.'''
        # Check to see if time is subsetted
        dt_start = self._get_datetime(start_time, get_start=True)
        dt_end = self._get_datetime(end_time, get_end=True)
        vsub = var.copy()
        vsub['data'] = vsub['data'][(self.time['data'][:] >= dt_start) & (self.time['data'][:] <= dt_end)]
        return vsub

    def _get_matchdata_by_index(self, indlon, indlat):
        '''Retrieve data fields using calculated indices.'''
        for field in self.data_fields.keys():
            self.matchdata[field] = self.data_fields[field].copy()
            self.matchdata[field]['data'] = self.data_fields[field]['data'][indlon, indlat]
        return

    def _get_matchdata_by_pyart_index(self, radar, indel, indaz, indrng):
        '''Retrieve data fields using calculated Py-ART indices.'''
        for field in radar.fields.keys():
            self.matchdata[field] = radar.fields[field].copy()
            self.matchdata[field]['data'] = radar.fields[field]['data'][azin, rgin]