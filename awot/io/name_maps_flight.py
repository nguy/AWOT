"""
awot.io.name_maps_flight
========================

Routines for mapping flight level data to
AWOT objects. These are default name maps
and are not guaranteed to work as they are
based on specific project data generally.
No real standardization currently exists in
the airborne research community.

Thse maps may be a good starting point for creating a
custom name_map.
"""

from __future__ import print_function


def _get_name_map(platform):
    """
    Retrieve a name_map used to map file variables to
    AWOT object.
    """
    # Set up lists of potential names for various aircraft
    p3_names = ['p3', 'p-3', 'noaap3', 'aoc',
                'noaa42', 'noaa43', 'gv']

    uwka_names = ['uwka', 'uwkingair', 'kingair',
                  'n2uw', 'raf']

    latmos_names = ['safire', 'latmos', 'french_falcon']

    t28_names = ['t28', 't-28', 'sdsmt', 'sdsmtt28', 'sdsm&t',
                 'sdsm&tt-28', 'sdsmtt-28', 'sdsm&tt28']

    und_citation_names = ['citation', 'und_citation',
                          'und', 'undcitation']

    if platform is not None:
        if platform.lower().replace(" ", "") in p3_names:
            name_map = _p3_flight_namemap()

        elif platform.lower().replace(" ", "") in uwka_names:
            name_map = _uwka_name_map()

        elif platform.lower().replace(" ", "") in latmos_names:
            name_map = _latmos_name_map()

        elif platform.lower().replace(" ", "") in t28_names:
            name_map = _t28_flight_namemap()

        elif platform.lower().replace(" ", "") in und_citation_names:
            name_map = _und_citation_name_map()
    return name_map


###########################
#  NAME MAP DICTIONARIES  #
###########################


def _p3_flight_namemap():
    '''
    Map NOAA P3 variables to AWOT structure.
    This map was developed based upon data from the
    DYNAMO 2011 project.

    Available variables (not full list) :
     LonGPS.3      = Novatel GPS Longitude
     LatGPS.3      = Novatel GPS Latitude
     AltGPS.3      = Novatel GPS Altitude [m]
     THdgI-GPS.1   = True heading [deg]
     TRK.1         = Track [deg]
     AltPaADDU.1   = Pressure altitude [m]
     WSZ_DPJ       = Vertical wind via D Jorgensen calculation [m/s]
     TA.1          = Ambient Temperature [C]
     TD.1          = Dewpoint Temperature [C]
     TVIRT.1       = Virtual Temperature [K]
     THETA.1       = Potential Temperature [K]
     THETAE.1      = Equivalent Potential Temperature [K]
     THETAV.1      = Virtual Potential Temperature [K]
     WS.1          = Wind Speed [m/s]
     WD.1          = Wind Direction [deg]
     HUM_REL.1     = Relative Humidity [%]
     HUM_SPEC.1    = Specific Humidity [g/kg]
     MR.1          = Mixing ratio [g] [g/g?]
     EE.1          = Vapor Pressure [hPa]
     EW.1          = Saturated Vapor Pressure [hPa]
    '''
    name_map = {
        'latitude': 'LatGPS.3',
        'longitude': 'LonGPS.3',
        'altitude': 'AltGPS.3',
        'pressure_altitude': 'AltPaADDU.1',
        'true_heading': 'THdgI-GPS.1',
        'track': 'TRK.1',
        'vert_vel_DPJ': 'WSZ_DPJ',
        'temperature': 'TA.1',
        'dewpoint_temperature': 'TD.1',
        'virtual_temperature': 'TVIRT.1',
        'potential_temp': 'THETA.1',
        'equiv_potential_temp': 'THETAE.1',
        'virtual_potential_temp': 'THETAV.1',
        'wind_spd': 'WS.1',
        'wind_dir': 'WD.1',
        'relative_humidity': 'HUM_REL.1',
        'specific_humidity': 'HUM_SPEC.1',
        'mixing_ratio': 'MR.1',
        'vapor_pressure': 'EE.1',
        'sat_vapor_pressure': 'EW.1', }
    return name_map


def _uwka_name_map():
    '''Map UWyo King Air variables to AWOT.'''
    name_map = {
        'time': 'time',
        # Aircraft Position
        'longitude': 'LONC',
        'latitude': 'LATC',
        'altitude': 'ztrue',
        'pressure_altitude': 'PALT',
        'tas': 'tas',
        'ias': 'aias',
        'true_heading': 'AVthead',
        'pitch': 'AVpitch',
        'roll_angle': 'AVroll',
        # Atmospheric State
        'pressure': 'pmb',
        'temperature': 'trf',
        'dewpoint_temperature': 'tdplicor',
        'thetad': 'thetad',
        'thetae': 'thetae',
        'relative_humidity': 'rh',
        'mixing_ratio': 'mr',
        'lwc': 'lwc100',
        'turb': 'turb',
        # Radiometric
        'irtop': 'irtc',
        'irbottom': 'irbc',
        'swtop': 'swt',
        'swbottom': 'swb',
        # Wind derivations
        'Uwind': 'AVuwind',
        'Vwind': 'AVvwind',
        'Wwind': 'AVwwind',
        'longitudinal_wind': 'AVux',
        'latitudinal_wind': 'AVvy',
        'wind_dir': 'AVwdir',
        'wind_spd': 'AVwmag',
        # Licor Concentrations
        'co2_conc': 'co21s',
        'h2o_conc': 'h2o1s',
        # Aerosol
        'pcasp_num': 'AS200_OBR',
        'pcasp_conc': 'CS200_OBR',
        'pcasp_mean_diam': 'DBARP_OBR',
        'pcasp_surf_area_conc': 'PSFCP_OBR',
        'pcasp_vol_conc': 'PVOLP_OBR',
        # Cloud Physics
        'conc_cpc': 'cpc_conc',
        # Miscellaneous
        'topo': 'topo',
    }
    return name_map


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


def _latmos_name_map():
    '''
    Map out names used in SAFIRE/LATMOS Falcon data to AWOT.

    This map was developed based upon data from the DYNAMO
    2011 project.
    '''
    name_map = {
        "time": "time",
        "latitude": "latitude : from GPS (degree)",
        "longitude": "longitude : from GPS (degree)",
        "altitude": "altitude : from GPS (meter)",
        "true_heading": "platform_orientation : from INS (degree)",
        "temperature":
        "air_temperature : from deiced Rosemount sensor (Celsius)",
        "dewpoint_temperature": "dew_point_temperature : from 1011B "
                                "top dew-point hygrometer (Celsius)",
        "wind_spd": "wind_speed : Attitudes and speed wrt ground from INS, "
                    "air angles from radome, air speed from pitot (m/s)",
        "wind_dir": "wind_from_direction : Attitudes and speed wrt ground "
                    "from INS, air angles from radome, air speed from "
                    "pitot (degree)",
        "relative_humidity": "relative_humidity : from Aerodata sensor (%)",
        "mixing_ratio": "humidity_mixing_ratio : from Aerodata sensor "
                        "(gram/kg)",
        "pressure": "air_pressure : from front sensor, corrected for the "
                    "so-called static defect (hPa)",
        "roll": "platform_roll_angle : from INS (degree)",
        "pitch": "platform_pitch_angle : from INS (degree)",
        "aircraft_air_speed": "platform_speed_wrt_air : from pitot (m/s)",
        "platform_ground_speed": "platform_speed_wrt_ground : from GPS (m/s)",
        "platform_ground_speed2": "platform_speed_wrt_ground : from INS (kt)",
        "aircraft_vert_accel": "platform_acceleration_along_vertical_axis : "
                               "from INS (meter second-2)",
        "altitude2": "altitude : from INS (meter)",
        "mixing_ratio2": "humidity_mixing_ratio : from top dew-point "
                         "hygrometer (GE 1011B) (gram/kg)",
        "platform_course": "platform_course : from INS (degree)",
        "platform_upward_ground_speed": "upward_platform_speed_wrt_ground : "
                                        "from INS (m/s)",
        "platform_upward_ground_speed": "upward_platform_speed_wrt_ground : "
                                        "from GPS (m/s)",
        "attack_angle": "angle_of_attack : from sensor on the boom "
                        "(degree)",
        "sideslip_angle": "angle_of_sideslip : from sensor on the "
                          "boom (degree)",
        "Uwind": "eastward_wind : Attitudes and speed wrt ground from "
                 "INS, air angles from radome, air speed from pitot (m/s)",
        "Vwind": "northward_wind : Attitudes and speed wrt ground from "
                 "INS, air angles from radome, air speed from pitot (m/s)",
        "air_vertical_velocity": "upward_air_velocity : Attitudes and "
                                 "speed wrt ground from INS, air angles "
                                 "from radome, air speed from pitot (m/s)",
        "wind_dir2": "wind_from_direction : Attitudes and speed wrt "
                     "ground from INS, air angles from radome, air speed "
                     "from pitot (degree)",
        "wind_spd2": "wind_speed : Attitudes and speed wrt ground from "
                     "INS, air angles from radome, air speed from pitot (m/s)",
    }
    return name_map


def _und_citation_name_map():
    '''
    Map out names used in UND Citation data to AWOT.

    This map was developed based upon data from the
    MC3E 2011 project.
    '''
    name_map = {
        'time': 'time',
        'temperature': "Air Temperature Corrected for Dynamic Heating "
                       "(Based first on the main temperatue/pitot "
                       "instrument and secondarly based on the backup"
                       " temperature/pitot instrument) [degC]",
        'mach_number': "Mach Number (Based first on the main "
                       "temperatue/pitot instrument and secondarly "
                       "based on the backup temperature/pitot instrument)",
        'ias': "Indicated Air Speed (Based first on the main pitot "
               "instrument and secondarly based on the backup pitot "
               "instrument) [m/s]",
        'tas': "True Air Speed (Based first on the main temperatue/pitot "
               "instrument and secondarly based on the backup "
               "temperature/pitot instrument) [m/s]",
        'pressure_altitude': "Pressure Altitude [m]",
        'theta': "Potential Temperature (Based first on the main "
                 "temperatue/pitot instrument and secondarly based "
                 "on the backup temperature/pitot instrument) [degK]",
        'pitot_pressure': "Pitot Pressure from Wing Probe [hPa] "
                          "{Calibration:  slope = 60.797701 "
                          "offset = -150.48263}",
        'cabin_pressure': "Aircraft Cabin Pressure [millibar]",
        'static_pressure': "Static Pressure [hPa] {Calibration:  "
                           "slope = 207.08000 offset = -0.71000000}",
        'dewpoint_temperature1': "Dewpoint Temperature from EG&G Probe "
                                 "[degC] {Calibration:  slope = 20.000000 "
                                 "offset = -70.000000}",
        'mixing_ratio': "Mixing Ratio by weight from the Laser "
                        "Hygrometer [ppmw]",
        'dewpoint_temperature2': "Dew Point Temperature from the Laser "
                                 "Hygrometer [degrees Celsius]",
        'frost_point_temperature': "Frost Point Temperature from the "
                                   "Laser Hygrometer [degrees Celsius]",
        'roll': "Aircraft Roll Angle from the Applanix POS System [degrees]",
        'pitch': "Aircraft Pitch Angle from the Applanix POS System [degrees]",
        'true_heading': "Aircraft Heading Angle [degrees]; 0-360 with 0 "
                        "being North",
        'aircraft_vert_accel': "Aircraft Z-direction (Vertical) Acceleration "
                               "for the Applanix POS system [m/s^2]",
        'latitude': "Aircraft Latitude from the Applanix POS System [degrees]",
        'longitude': "Aircraft Longitude from the Applanix POS system "
                     "[degrees]",
        'altitude': "Aircraft Altitude from the Applanix POS system [m]",
        'aircraft_speed': "Aircraft Speed from the Applanix POS system [m/s]",
        'track': "Aircraft Track Angle [degrees]; 0-360 with 0 "
                 "being North",
        'attack_angle': "Alpha (Attack) Angle [degrees] "
                        "{Calibration:  slope = 0.070310700 "
                        "offset = 0.33308870}",
        'sideslip_angle': "Beta (Sideslip) Angle [degrees] "
                          "{Calibration:  slope = -0.073846080 "
                          "offset = -1.6202790}",
        'aicrcaft_vert_vel': "Vertical Velocity of the aircraft based on "
                             "the change in position over a 2 seceond "
                             "interval [m/s]",
        'Wwind': "Z (Vertical) Component of the Wind Speed [m/s]",
        'wind_spd': "Horizontal Wind Speed [m/s]",
        'wind_dir': "Horizontal Wind Direction [degrees]; True Direction "
                    "From Which it Blows",
        'turb': "Turbulence parameter (Eddy Dissipation Rate) based on Wing "
                "Pitot pressure [cm^2/3*s^-1]",
        'lwc1': "Liquid Water Content based on King Probe measurement "
                "adjusted (cloud threshold =  5.1 [#/cm^3], cloud interval "
                "= 30.0 [s] and adjustment slope = 0.500) for the baseline "
                "offset [g/m^3]",
        'twc': "Total Water Content based on the Nevzorov Probe measurement",
        'lwc2': "Liquid Water Content based on the Nevzorov Probe measurement "
                "with correction for  residual ice (beta = 0.110000) [g/m^3]",
        'Conc_CDP': "Number Concentration of Droplets Based on the Cloud "
                    "Droplet Probe [#/cc]",
        'lwc3': "Liquid Water Content Based on the Cloud Droplet Probe "
                "[g/m^3]",
        'Dmean_CDP': "Cloud Droplet Probe's Mean Droplet Diameter [um]",
        'Dvol_CDP': "Cloud Droplet Probe's Mean Droplet Volume Diameter [um]",
        'Deff_CDP': "Cloud Droplet Probe's Effective Droplet Radius [um]",
        'Conc_2DC': "Number concentration of droplets based on the 2-DC "
                    "Probe measurements [#/cm^3]",
        'Dmean_2DC': "Mean droplet diameter based on the 2-DC Probe "
                     "measurements [um]",
        'Dvol_2DC': "Mean droplet volume diameter based on the 2-DC "
                    "Probe measurements [um]",
        'Deff_2DC': "Effective droplet radius based on the 2-DC "
                    "Probe measurements [um]",
        'mso_frequency': "The current Sensor (MSO) frequency from the "
                         "Icing Detector [Hz]",
        'Conc_CPC': "Total Concentration from CPC [#/cm^3]",
        'yymmdd': "Date Stamp based on data file name "
                  "(Example: 941119 is 19 November 1994) [stamp]",

    }
    return name_map
