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
    """Class for flight level plots."""

    def __init__(self, flightdata, radar, basemap=None,
                 flight_lon_name=None, flight_lat_name=None,
                 flight_alt_name=None, flight_time_name=None,
                 flight_uwind_name=None, flight_vwind_name=None,
                 radar_lon_name=None, radar_lat_name=None,
                 radar_height_name=None, radar_time_name=None,
                 radar_match_dict=None,
                 ):
        '''
        Intitialize the class to create plots.

        Parameters
        ----------
        flightdata : dict
            AWOT flight data dictionary instance.
        radar : dict
            AWOT radar instance.
        basemap : basemap instance
            Basemap instance to use for plotting.
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
        radar_lon_name : str
            Key in radar instance for longitude variable.
            None uses AWOT default.
        radar_lat_name : str
            Key in radar instance for latitude variable.
            None uses AWOT default.
        radar_height_name : str
            Key in radar instance for height variable.
            None uses AWOT default.
        radar_time_name : str
            Key in radar instance for time variable.
            None uses AWOT default.
        radar_match_dict : list
            List of field names to match for radar object.
            None loops through existing fields.
        '''
        if flight_lon_name is None:
            self.flight_longitude = flightdata['longitude']
        else:
            self.flight_longitude = flightdata[flight_lon_name]
        if flight_lat_name is None:
            self.flight_latitude = flightdata['latitude']
        else:
            self.flight_latitude = flightdata[flight_lat_name]
        if flight_alt_name is None:
            self.flight_altitude = flightdata['altitude']
        else:
            self.flight_atlitude = flightdata[flight_alt_name]
        if flight_time_name is None:
            self.flight_time = flightdata['time']
        else:
            self.flight_time = flightdata[flight_time_name]
        if flight_uwind_name is None:
            self.flight_Uwind = flightdata['Uwind']
        else:
            self.flight_Uwind = flightdata[flight_uwind_name]
        if flight_vwind_name is None:
            self.flight_Vwind = flightdata['Vwind']
        else:
            self.flight_Vwind = flightdata[flight_vwind_name]
        self.flight_number = flightdata['flight_number']
        self.project = flightdata['project']
        self.platform = flightdata['platform']
        self.flight_data = flightdata

        self.radar = radar
        self.radar_fields = self.radar['fields']

        if radar_lon_name is None:
            self.radar_longitude = self.radar['longitude']
        else:
            self.radar_longitude = self.radar[radar_lon_name]
        if radar_lat_name is None:
            self.radar_latitude = self.radar['latitude']
        else:
            self.radar_latitude = self.radar[radar_lat_name]
        if radar_height_name is None:
            self.radar_height = self.radar['height']
        else:
            self.radar_height = self.radar[radar_height_name]
        if radar_time_name is None:
            self.radar_time = self.radar['time']
        else:
            self.radar_time = self.radar[radar_time_name]
        self.basemap = basemap
        _check_basemap(self, strong=True)

        # Calculate x,y map position coordinates
        self.flight_x, self.flight_y = self.basemap(
             self.flight_longitude['data'][:], self.flight_latitude['data'][:])

        # Convert lats/lons to 2D grid
        Lon2D, Lat2D = np.meshgrid(self.radar_longitude['data'][
                                   :], self.radar_latitude['data'][:])
        # Convert lats/lons to map projection coordinates
        self.radar_x, self.radar_y = self.basemap(Lon2D, Lat2D)


    def weighted_stats(threshold_time=180, threshold_hor=250,
                     threshold_vert=125, threshold_ac=6):
        """
        Match the radar field values to each aicraft time
        """
        Z14_wmean, Z35_wmean, DFR_wmean = [], [], []
        Z14_wstd, Z35_wstd, DFR_wstd, dr_std, dt_std = [], [], [], [], []
        dt, dr, timestamp, actime, height, longheight, matchtimes = [], [], [], [], [], [], []
        twctrail, lwctrail = [],[]
        for i in range(len(self.flight_time['data'][:])):
            timespan = range(int(self.flight_time['data'][i] - threshold_time),
                                 int(self.flight_time['data'][i] + threshold_time))
            rdex = np.in1d(self.radar_time['data'][:], timespan)
            if np.mod(i, 1000) == 0:
                print(i)
            if np.any(self.radar_time['data'][:][rdex]):
                field = np.dstack([self.radar_x[rdex, :, :].ravel(),
                                   self.radar_y[rdex, :, :].ravel(),
                                   self.radar_height['data'][rdex, :, :].ravel()])[0]
                tree = scipy.spatial.cKDTree(field)
                goahead = True
                predist, predex = tree.query([self.flight_x[i],
                                              self.flight_y[i],
                                              self.flight_altitude['data'][i]], 1)
                if np.any(predex):
#***                    loctime = itemgetter(predex)(radtime_3D[rdex,:,:].ravel())
                    dtime = self.flight_time['data'][i] - loctime
                    trux = self.flight_x[i] - self.flight_Uwind['data'][i] * dtime
                    truy = self.flight_y[i] - self.flight_Vwind['data'][i] * dtime
                    truz = self.flight_altitude['data'][i]
                    closest, matchdexa = tree.query([trux[0], truy[0], truz], 10)
                    if np.max(closest) > 300:
                        goahead = False
                    else:
                        matchtimes = np.hstack((matchtimes, self.flight_time['data'][:][i]*np.ones(np.shape(matchdexa))))
#***BELOW HERE WOULD BE WHERE LOOP THROUGH LIST OF FIELDS WOULD OCCUR??
                    ###Retrieve reflectivities with weighted statistics
                    if goahead:
                        matchdex = matchdexa.astype(int)
                        branch14 = itemgetter(matchdex)(flec14[rdex, :, :].ravel())
                        branch35 = itemgetter(matchdex)(flec35[rdex, :, :].ravel())
                        timebranch = itemgetter(matchdex)(radtime_3D[rdex, :, :].ravel())
                        if np.all((branch14 > -100) & (branch35 > -100)):
                            if np.any(branch14 > 30):
                                branch14 = branch14[branch14 < 30]
                                branch35 = branch14[branch14 < 30]
                            DFRbranch = branch14-branch35
                            print('match')
                            Z14_wmean.append(np.mean(branch14))
                            Z35_wmean.append(np.mean(branch35))
                            DFR_wmean.append(np.mean(DFRbranch))
                            Z14_wstd.append(np.sqrt(np.mean((branch14-Z14_wmean[-1])**2)))
                            Z35_wstd.append(np.sqrt(np.mean((branch35-Z35_wmean[-1])**2)))
                            DFR_wstd.append(np.sqrt(np.mean((DFRbranch-DFR_wmean[-1])**2)))
                            timestamp.append(np.mean(timebranch))
                            actime.append(self.flight_time['data'][i])
                            dt.append(np.mean(np.abs(matchtimes - timestamp[-1])))
                            dt_std.append(np.std(np.abs(matchtimes - timestamp[-1])))
                            twctrail.append(TWC[i])
                            lwctrail.append(LWC[i])

        try:
            finaldex = np.in1d(self.flight_time['data'][:], actime)
            outSD = SD[finaldex,:]
        except:
            pdb.set_trace()
        return match_data

    def nearest_neighbor_pyart():
        """
        """
#        condition = np.logical_and(asod >= start_time, asod < end_time)
#    indices = np.where(condition)
#    print start_time, end_time, np.size(indices[0])
#    if np.size(indices[0]) > 0:
#        fname = os.path.basename(files[iii])
#        print fname
#        radar = pyart.io.read(files[iii])
            for index in indices[0]:
                elin = np.argmin(np.abs(kel_kgwx[index]-radar.fixed_angle['data']))
                radar_sweep = radar.extract_sweeps([elin])
                azin = np.argmin(np.abs(kbear_kgwx[index]-radar_sweep.azimuth['data']))
                rgin = np.argmin(np.abs(ksr_kgwx[index]-radar_sweep.range['data']/1000.0))
                tpos1 = sweep_latlon_to_flat_xy(radar, elin)
                cond = get_circle_condition(radar, tpos1, atlat[index], atlon[index], 8.0)
                ts_dz_max[index], ts_dz_min[index] = max_min_given_condition(radar_sweep, 'reflectivity', cond)
                ts_edr_max[index], ts_edr_min[index] = max_min_given_condition(radar_sweep, 'turbulence', cond)
                print asod[index]/3600.0, kel_kgwx[index], ksr_kgwx[index], kbear_kgwx[index],\
                      radar_sweep.fields['reflectivity']['data'][azin, rgin],\
                      radar_sweep.fields['turbulence']['data'][azin, rgin],\
                      radar_sweep.fields['velocity']['data'][azin, rgin],\
                      radar_sweep.fields['spectrum_width']['data'][azin, rgin],\
                      ts_dz_max[index], ts_dz_min[index], ts_edr_max[index], ts_edr_min[index]
                ts_dz[index] = radar_sweep.fields['reflectivity']['data'][azin, rgin]
                ts_edr[index] = radar_sweep.fields['turbulence']['data'][azin, rgin]
                ts_vel[index] = radar_sweep.fields['velocity']['data'][azin, rgin]
                ts_sw[index] = radar_sweep.fields['spectrum_width']['data'][azin, rgin]

    def spline_interp():
        """
        """

        return match_data