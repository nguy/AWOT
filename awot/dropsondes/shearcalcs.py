import numpy as np

class ShearCalcs:
    '''A class for shear calculations.
    All calculations can be found in Lucas et al. (2000); Journal of 
    Atmospheric Sciences.
    '''
    def __init__(self):
        self.u_6km = []
        self.v_6km = []
            
    def _VertShear_Sfc_to_6km(self, Height, U, V):
        '''
        Calculate vertical shear [s^-1] using the 
        Weismann and Klemp (1982, MWR) methodology.
        This module calculates from the surface to 3 km AGL.  

        Parameters::
        ----------
        Height : float
            Altitude of observations [m]
        U : float
            Zonal Wind speed [m/s]
        V : float
            Meridional Wind speed [m/s]
        '''
        index6km = 6000.0
        sfcindex = 0
    
        U6km = np.interp(index6km, Height, U)
        V6km = np.interp(index6km, Height, V)
   
        U= U.compressed()
        V= V.compressed()
        Height = Height.compressed()
    
        Hsfc = Height[sfcindex]
        Usfc = U[sfcindex]
        Vsfc = V[sfcindex]
         
        shear = np.sqrt(((U6km - Usfc)**2 + (V6km - Vsfc)**2)/(index6km))
        
        return shear
        
    def _VertShear_Sfc_to_3km(self, Height, U, V):
        '''
        Calculate vertical shear [s^-1] using the 
        Weismann and Klemp (1982, MWR) methodology.
        This module calculates from the surface to 1 km AGL.  

        Parameters::
        ----------
        Height : float
            Altitude of observations [m]
        U : float
            Zonal Wind speed [m/s]
        V : float
            Meridional Wind speed [m/s]
        '''

        index3km = 3000.0
        sfcindex = 0

        U3km = np.interp(index3km, Height, U)
        V3km = np.interp(index3km, Height, V)

        U= U.compressed()
        V= V.compressed()
        Height = Height.compressed()

        Hsfc = Height[sfcindex]
        Usfc = U[sfcindex]
        Vsfc = V[sfcindex]

        shear = np.sqrt(((U3km - Usfc)**2 + (V3km - Vsfc)**2)/(index3km))
        
    
        return shear
            
        
    def _VertShear_Sfc_to_1km(self, Height, U, V):
    
        index1km = 1000.0
        sfcindex = 0
        
        U1km = np.interp(index1km, Height, U)
        V1km = np.interp(index1km, Height, V)
        
        
        print('test at 800 mb recorded 30 knot')
        
        u800 = np.interp(index1km, Height, U)
        v800 = np.interp(index1km, Height, U)


        U= U.compressed()
        V= V.compressed()
        Height = Height.compressed()

        Hsfc = Height[sfcindex]
        Usfc = U[sfcindex]
        Vsfc = V[sfcindex]
        
        '''
        print(U)
        print(V)
        print(U1km)
        print(V1km)
        '''
        
        shear = np.sqrt(((U1km - Usfc)**2 + (V1km - Vsfc)**2)/(index1km))
        
        return shear
            
    def _bulkshear_sfc_1km(self, Height, U, V):
        index1km = 1000.0
        sfcindex =0

        U1km = np.interp(index1km, Height, U)
        V1km = np.interp(index1km, Height, V)

        U = U.compressed()
        V = V.compressed()
        Height = Height.compressed()

        HSFC = Height[sfcindex]
        USFC = U[sfcindex]
        VSFC = V[sfcindex]

        bulk_shear = np.sqrt((U1km-USFC)**2 + (V1km-VSFC)**2)
        
        print (np.sqrt((U1km)**2 + (V1km)**2))

        return bulk_shear    
        
    def _bulkshear_sfc_3km(self, Height, U, V):
        index3km = 3000.0
        sfcindex =0

        U3km = np.interp(index3km, Height, U)
        V3km = np.interp(index3km, Height, V)

        U = U.compressed()
        V = V.compressed()
        Height = Height.compressed()

        HSFC = Height[sfcindex]
        USFC = U[sfcindex]
        VSFC = V[sfcindex]

        bulk_shear = np.sqrt((U3km - USFC)**2 + (V3km - VSFC)**2)

        return bulk_shear
        
    def _bulkshear_sfc_6km(self, Height, U, V):
        index6km = 6000.0
        sfcindex = 0

        U6km = np.interp(index6km, Height, U)
        V6km = np.interp(index6km, Height, V)

        U = U.compressed()
        V = V.compressed()
        Height = Height.compressed()

        HSFC = Height[sfcindex]
        USFC = U[sfcindex]
        VSFC = V[sfcindex]
        
        print (np.sqrt((U6km)**2 + (V6km)**2))

        bulk_shear = np.sqrt((U6km - USFC)**2+(V6km - VSFC)**2)

        return bulk_shear    
    