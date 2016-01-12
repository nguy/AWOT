"""
awot.io.read_t28
================

Function for reading T-28 data.
See awot.io.flight.read_netcdf for a generalized version.

"""

from __future__ import absolute_import, print_function
import numpy as np
from netCDF4 import num2date
from .flight import _winduv
from . import common

def read_t28_netcdf(ncFile):
    """
    T-28 has enough differences from other aircraft netCDF files to
    warrant its own read script.

    Key differences are the large number of non-standard variables
    (e.g., from the electric field meters) and how time is handled by
    in the netCDF file.
    """
    name_map = _t28_flight_namemap()
    data = _make_data_dictionary(ncFile, name_map)
    data['time'] = _get_time(ncFile)
    if 'turb' in data.keys():
        CM2M_TWOTHIRDS = 100**0.6666667
        data['turb'] /= CM2M_TWOTHIRDS

    if 'Uwind' not in name_map:
        Uwind, Vwind = _winduv(data)
        # Add to the dictionary
        data['Uwind'] = Uwind
        data['Vwind'] = Vwind

    data['project'] = ''
    data['platform'] = 'SDSM&T T-28'
    if hasattr(ncFile, 'FlightNumber'):
        data['flight_number'] = ncFile.FlightNumber
    else:
        data['flight_number'] = ''

    return data


def _make_data_dictionary(ncFile, name_map):
    """
    T-28 provides 2D arrays (20-Hz sampling reporting every second),
    and these need to be flattened.
    """
    data = {}
    for var in name_map:
        tmp = np.array(ncFile.variables[name_map[var]]).ravel()
        data[var] = tmp.ravel()
        try:
            np.ma.masked_invalid(data[var])
        except:
            data[var] = None
        try:
            mask = data[var].mask
        except:
            data[var] = np.ma.masked_array(data[var], mask=False)
    return data


def _get_time(ncFile):
    # Now convert the time array into a datetime instance
    dtHrs = num2date(np.array(ncFile.variables['Time'][:]).ravel(),
                     ncFile.variables['Time'].units)

    # Now convert this datetime instance into a number of seconds since Epoch
    TimeSec = date2num(dtHrs, common._get_epoch_units())

    Time_unaware = num2date(TimeSec, common._get_epoch_units())
    Time = {'data': Time_unaware, 'units': common._get_epoch_units(),
            'title': 'Time', 'full_name': 'Time (UTC)'}
    return Time_unaware


def _t28_flight_namemap():
    '''
    Map SDSM&T T-28 variables to AWOT structure

    Variables that will be mapped to new names (not full list):
     LONGITUDE_DECIMAL_DEG_20Hz = GPS Longitude
     LATITUDE_DECIMAL_DEG_20Hz = GPS Latitude
     GPS_ALTITUDE = GPS Altitude [m]
     GPS_GROUND_TRACK_ANGLE = Track [deg]
     PRESSURE_ALTITUDE = Pressure altitude [m]
     TEMPERATURE_EQUIVALENT_POTENTIAL = Equivalent Pot. Temperature [K]
     INDICATED_AIRSPEED = Indicated Airspeed [m/s]
     TURBULENCE = Cubic Root of Eddy Dissipation Rate [cm^(2/3) s^-1]
     ROLL = Roll angle [deg]
     PITCH = Pitch angle [deg]
     TRUE_AIRSPEED_CALCULATED = True Airspeed [m/s]
    '''
    name_map = {
        'latitude': 'LATITUDE_DECIMAL_DEG_20Hz',
        'longitude': 'LONGITUDE_DECIMAL_DEG_20Hz',
        'altitude': 'GPS_ALTITUDE',
        'pressure_altitude': 'PRESSURE_ALTITUDE',
        'track': 'GPS_GROUND_TRACK_ANGLE',
        'equiv_potential_temp': 'TEMPERATURE_EQUIVALENT_POTENTIAL',
        'roll_angle': 'ROLL',
        'pitch': 'PITCH',
        'ias': 'INDICATED_AIRSPEED',
        'turb': 'TURBULENCE',
        'tas': 'TRUE_AIRSPEED_CALCULATED',
        'UPDRAFT': 'UPDRAFT',
        'PRESSURE_STATIC_1': 'PRESSURE_STATIC_1',
        'PRESSURE_STATIC_2': 'PRESSURE_STATIC_2',
        'TEMPERATURE_ROSEMOUNT_SENSOR': 'TEMPERATURE_ROSEMOUNT_SENSOR',
        'TEMPERATURE_REVERSE_FLOW_SENSOR': 'TEMPERATURE_REVERSE_FLOW_SENSOR',
        'DENSITY_AIR': 'DENSITY_AIR',
        'NO_CONCENTRATION': 'NO_CONCENTRATION',
        'FSSP_LIQUID_WATER': 'FSSP_LIQUID_WATER',
        'FSSP_TOTAL_COUNTS': 'FSSP_TOTAL_COUNTS',
        'FSSP_AVERAGE_DIAMETER': 'FSSP_AVERAGE_DIAMETER',
        'FSSP_TOTAL_PARTICLE_CONCENTRATION':
            'FSSP_TOTAL_PARTICLE_CONCENTRATION',
        'FSSP_EQUIVALENT_DIAMETER': 'FSSP_EQUIVALENT_DIAMETER',
        'FSSP_EQUIVALENT_DIAMETER_VARIANCE':
            'FSSP_EQUIVALENT_DIAMETER_VARIANCE',
        'FSSP_LIQUID_WATER_MIXING_RATIO': 'FSSP_LIQUID_WATER_MIXING_RATIO',
        'FSSP_GATED_STROBES': 'FSSP_LIQUID_WATER_MIXING_RATIO',
        'FSSP_TOTAL_STROBES': 'FSSP_TOTAL_STROBES',
        'FSSP_ACTIVITY': 'FSSP_ACTIVITY',
        'HAIL_WATER': 'HAIL_WATER',
        'HAIL_TOTAL_COUNTS': 'HAIL_TOTAL_COUNTS',
        'HAIL_AVERAGE_DIAMETER': 'HAIL_AVERAGE_DIAMETER',
        'HAIL_CONCENTRATION': 'HAIL_CONCENTRATION',
        'HAIL_LIQUID_WATER_MIXING_RATIO': 'HAIL_LIQUID_WATER_MIXING_RATIO',
        'SHADOW_OR_PMS': 'SHADOW_OR_PMS',
        'PMS_END_ELEMENT_1': 'PMS_END_ELEMENT_1',
        'PMS_END_ELEMENT_2': 'PMS_END_ELEMENT_2',
        'LWC_DMT': 'LWC_DMT',
        'PRESSURE_DYNAMIC_1': 'PRESSURE_DYNAMIC_1',
        'PRESSURE_DYNAMIC_2': 'PRESSURE_DYNAMIC_2',
        'ACCELERATION_X': 'ACCELERATION_X',
        'ACCELERATION_Y': 'ACCELERATION_Y',
        'ACCELERATION_Z': 'ACCELERATION_Z',
        'ACCELERATION_VERTICAL': 'ACCELERATION_VERTICAL',
        'DZDT_POINT': 'DZDT_POINT',
        'TIME_GPS_DECIMAL': 'TIME_GPS_DECIMAL',
        'GPS_RATE_OF_CLIMB': 'GPS_RATE_OF_CLIMB',
        'FIELD_METER_TOP': 'FIELD_METER_TOP',
        'FIELD_METER_BOTTOM': 'FIELD_METER_BOTTOM',
        'FIELD_METER_LEFT': 'FIELD_METER_LEFT',
        'FIELD_METER_RIGHT': 'FIELD_METER_RIGHT',
        'FIELD_METER_5': 'FIELD_METER_5',
        'FIELD_METER_6': 'FIELD_METER_6',
        'FIELD_METER_TOP': 'FIELD_METER_TOP',
        'EZ_AIRCRAFT': 'EZ_AIRCRAFT',
        'EY_AIRCRAFT': 'EY_AIRCRAFT',
        'EX_AIRCRAFT': 'EX_AIRCRAFT',
        'EQ_AIRCRAFT_CHARGE': 'EQ_AIRCRAFT_CHARGE',
        'GPS_GROUNDSPEED': 'GPS_GROUNDSPEED',
        'TRUE_AIRSPEED_NCAR': 'TRUE_AIRSPEED_NCAR',
        'EZ_PATH': 'EZ_PATH',
        'EY_PATH': 'EY_PATH',
        'EX_PATH': 'EX_PATH',
        'EZ_EARTH': 'EZ_EARTH',
        'EY_EARTH': 'EY_EARTH',
        'EX_EARTH': 'EX_EARTH',
        'GPS_GROUND_TRACK_ANGLE': 'GPS_GROUND_TRACK_ANGLE',
        'HEADING_MAGNETIC': 'HEADING_MAGNETIC',
        'GPS_MAGNETIC_DEVIATION': 'GPS_MAGNETIC_DEVIATION',
        'MANIFOLD_PRESSURE': 'MANIFOLD_PRESSURE',
        'VOLTAGE_REGULATOR': 'VOLTAGE_REGULATOR',
        'INTERIOR_TEMPERATURE': 'INTERIOR_TEMPERATURE',
        'HEATER_CURRENT': 'HEATER_CURRENT',
        'EVENT_MARKERS': 'EVENT_MARKERS'}
    return name_map
