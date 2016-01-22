import numpy as np
from .thermocalcs import ThermoCalcs
from .shearcalcs import ShearCalcs




def add_thermo_calcs(data):

    '''
    method to calculate thermodynamic parameters from dropsonde data

    Inputs:

    # Dictionary of sounding data.
    # Uses calculations from thermocalcs.py
    # ============

    output:

    # multiple arrays containing thermodynamic information
    '''

    
    '''
        method to calculate thermodynamic parameters from dropsonde data
        
        Inputs:
        
        Dictionary of sounding data.
        Uses calculations from thermocalcs.py
        
        '''
    
    tC = ThermoCalcs()
    
    T = data['temperature']['data'][:]
    Td = data['dewpoint']['data'][:]
    p = data['presssure']['data'][:]
    RH = data['relative_humidity']['data'][:]
    u = data['u_component']['data'][:]
    v = data['v_component']['data'][:]
    h = data['Height']['data'][:]

    LCLT = round((
    tC._LCL_temperature(h, T+273.15, Td+273.15)-273.15), 2)
    LCLP = round((
    tC._LCL_presssure(h, p, T+273.15, Td+273.15)), 0)
    LCLZ = round((
    tC._LCL_Height(h, p, T+273.15, Td+273.15)), 0)
    THETA = tC._PTk2_Theta(p, T+273.15)
    MIXR = tC._RH_2_MixR(RH, p, T+273.15)
    THETAE = tC._Tk_RH_MixR_2_ThetaE(
    p, T+273.15, RH, MIXR/1000.)
    ESAT = tC._esat(T+273.15)

    thermoCalcData = dict(
    thermoCalcData['lclt'] = _build_dict(LCLT, 'K', 'Lifting Condensation Level Temperature', 'LCL Temp')
    thermoCalcData['lclp'] = build_dict(LCLP, 'hPa', 'Lifting Condensation Level Pressure', 'LCL Pressure')
    thermoCalcData['lclz'] = build_dict(LCLZ, 'm', 'Lifting Condensation Level Height', 'LCL Height')
    thermoCalcData['theta'] = build_dict(THETA, 'K', 'Potential Temperature', 'Potential Temp')
    thermoCalcData['mixr'] = build_dict(MIXR, 'g/m^3', 'Mixing Ratio', 'Mixing Ratio')
    thermoCalcData['thetae'] = build_dict(THETAE, 'K', 'Equivalent Potential Temperature', 'Equiv Potential Temp')
    thermoCalcData['esat'] = build_dict(ESAT, 'Saturation Vapor Pressure', 'Saturation Vapor Pressure')


    return thermoCalcData



def add_shear_calcs(data):

    '''
    method to calculate thermodynamic parameters from dropsonde data

    Inputs:

    # Dictionary of sounding data.
    # Uses calculations from thermocalcs.py
    # ============

    output:

    # multiple arrays containing thermodynamic information

    '''

    T = data['temperature']
    Td = data['dewpoint']
    p = data['presssure']
    RH = data['relative_humidity']
    uwind = data['u_component']
    vwind = data['v_component']
    height = data['Height']

    mask = h.mask
    uwind = u[~mask]
    vwind = v[~mask]
    height = h[~mask]

    SHEAR1KM = sC._VertShear_Sfc_to_1km(h, u, v)
    SHEAR3KM = sC._VertShear_Sfc_to_3km(h, u, v)
    SHEAR6KM = sC._VertShear_Sfc_to_6km(h, u, v)
    BULKSHEAR1km = round(sC._bulkshear_sfc_1km(h, u, v), 2)
    BULKSHEAR3km = round(sC._bulkshear_sfc_3km(h, u, v), 2)
    BULKSHEAR6km = round(sC._bulkshear_sfc_6km(h, u, v), 2)


    shearCalcData = dict()
    shearCalcData['SHEAR1KM'] = SHEAR1KM
    shearCalcData['SHEAR3KM'] = SHEAR3KM
    shearCalcData['SHEAR6KM'] = SHEAR6KM
    shearCalcData['BULKSHEAR1km'] = BULKSHEAR1km
    shearCalcData['BULKSHEAR1km'] = BULKSHEAR3km
    shearCalcData['BULKSHEAR1km'] = BULKSHEAR6km

    return shearCalcData


#def dry_lift(data):
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



