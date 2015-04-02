import numpy as np
class ThermoCalcs:
	
	def __init__(self):

		#######################
		#  thermocalcs
		#######################
		# This is a grouping of modules various thermodynamic calculations.
		#
		# Some routines are based upon Atmospheric Physics routines in IDL
		# written by Dominik Brunner 
		# (http://www.iac.ethz.ch/staff/dominik/idltools/idl_atmosphys.html)
		#-------------------------------------------------------------------
		# Define various constants that may be used for calculations
		# May want to use an already established library?
		self.R = 8.3145         # Universal gas constant (J deg-1 kmol-1)
		self.k = 1.381E-23      # Botzmann constant (J deg-1 molecule-1)
		self.SIGM = 5.6696E-8   # Stefan-Boltzmann constant (W m-2 deg-4)
		self.GRAV = 9.80665        # Acceleration due to gravity (m s-2)
		self.E0 = 611.0         # Saturation vapor pressure at 0 deg c (Pa)
		self.MD = 28.966        # Molecular weight of dry air
		self.RD = 287.04        # Gas constant of dry air (=R*1000/MD)
		self.CPD = 1004.67      # Specific heat of dry air at const. p (J deg-1 kg-1)
		self.CVD = 717.63       # Specific heat of dry air at const. volume (J deg-1 kg-1)
		self.MW = 18.016        # Molecular weight of water
		self.RV = 461.40        # Gas constant for water vapor (J deg-1 kg-1)
		self.CPV = 1865.1       # Specific heat of water vapor at const. p (J deg-1 kg-1)
		self.CVV = 1403.2       # Specific heat of water vapor at const. volume (J deg-1 kg-1)
		self.SLP = 1013.25      # Sea level pressure (hPa)
		self.LV = 2.501E6       # Latent heat of vaporization (J kg-1) [*-1 for vaporization]
		self.LD = 2.834E6       # Latent heat of deposition (J kg-1) [*-1 for sublimation]
		self.LF = 3.34E5        # Latent heat of fusion (J kg-1) [*-1 for melting]
		self.KAPPA = self.RD/self.CPD     # Ratio of Gas constant to specific heat
		self.T0 = 273.15          # Freezing temperature in degrees K

		#===============================================================



	def _LCL_Temperature(self, Height, TempK, TdewK):
		'''
		Calculate the temperature at the lifting condensation level [K], 
		using averages of the lowest 500 m of sounding.  

		Defined by the University of Wyoming Sounding website 
		(http://weather.uwyo.edu/upperair/indices.html) as
		LCLT	= [1 / ( 1 / ( TdewK - 56 ) + LN ( TempK / TdewK ) / 800 )] + 56  

		Parameters::
		----------
		Pressure : float
			Pressure [hPa]
		Height : float
			Altitude of observations [m]
		TempK : float
			Temperature [K]
		Tdew : float
			Dewpoint temperature [K]
		'''
	
		# Make sure Height is a numpy array
		Height = np.array(Height)

		# Find the boolean subset of the lowest 500 m
		In500 = np.where(np.logical_and(Height > Height.min(),
										Height <= Height.min() + 500.))

		# Compute averages for the layer
		TempK_mean = np.mean(TempK[In500])
		TdewK_mean = np.mean(TdewK[In500])

		LCLT = (1. / (1. / (TdewK_mean - 56.) + 
			   (np.log(TempK_mean / TdewK_mean) / 800.))) + 56.
		return LCLT

	###############

	def _LCL_Pressure(self, Height, Pressure, TempK, TdewK):
		'''
		Calculate the pressure at the lifting condensation levl [hPa], 
		using averages of the lowest 500 m of sounding.  

		Defined by the University of Wyoming Sounding website 
		(http://weather.uwyo.edu/upperair/indices.html) as
		LCLP	= PRES * ( LCLT / TempK ) ** ( 1 / KAPPA )  

		Parameters::
		----------
		Pressure : float
			Pressure [hPa]
		Height : float
			Altitude of observations [m]
		TempK : float
			Temperature [K]
		Tdew : float
			Dewpoint temperature [K]
		'''
		LCLT = self._LCL_Temperature(Height, TempK, TdewK)

		# Make sure Height is a numpy array
		Height = np.array(Height)

		# Find the boolean subset of the lowest 500 m
		In500 = np.where(np.logical_and(Height > Height.min(),
										Height <= Height.min() + 500.))

		# Compute averages for the layer
		TempK_mean = np.mean(TempK[In500])
		TdewK_mean = np.mean(TdewK[In500])
		Pres_mean = np.mean(Pressure[In500])

		LCLP = Pres_mean * (LCLT / TempK_mean)**self.KAPPA
		return LCLP

	###############

	def _LCL_Height(self,Height, Pressure, TempK, TdewK):
		'''
		Calculate the height of the lifting condensation level [m],
		using the pressure level calculated by LCLPressure.
		The Hypsometric equation is applied to find the height at this level.

		H = (Rd*Temp/g) * (LN(p0 / LCLP))
		where R is gas constant, g is gravity, p0 is SLP 

		Parameters::
		----------
		Pressure : float
			Pressure in hPa
		Height : float
			Altitude of observations [m]
		TempK : float
			Temperature [K]
		Tdew : float
			Dewpoint temperature [K]
		'''
		LCLT = self._LCL_Temperature(Height, TempK, TdewK)
		LCLP = self._LCL_Pressure(Height, Pressure, TempK, TdewK)

		# Make sure Height is a numpy array
		Height = np.array(Height)

		# Find the boolean subset of the lowest 500 m
		In500 = np.where(np.logical_and(Height > Height.min(),
										Height <= Height.min() + 500.))

		# Compute averages for the layer
		Pres_mean=np.mean(Pressure[In500])

		LCLH=(self.RD * LCLT / self.GRAV)*(np.log(self.SLP / LCLP))
		return LCLH
		
	def _esat(self, TempK, Opt='variable'):
		'''
		Calculate saturation vapor pressure [hPa].
	
		Equations from Stull 2000 "Meteorology for Scientists and Engineers"

		Parameters::
		----------
		TempK : float
			Temperature [deg K]
		Opt : str
			Calculation option, possibilites are:
				'var' or 'variable' to use temperature dependence
				'water'
				'ice'
		'''
		# Choose the latent heat based upon temperature
		if (Opt.lower() == 'water'):
			L = self.LV
		elif (Opt.lower() == 'ice'):
			L = self.LD
		elif (Opt.lower() == 'var') or (Opt.lower() == 'variable'):#
			L = np.empty_like(TempK)
			L[TempK >= -268.15] = self.LV
			L[TempK < -268.15] = self.LD
  
		# Compute sat vapor pressure, divide by 100 to convert to hPa
		esv = self.E0 * np.exp(L / self.RV * ((1. / self.T0) - (1. / TempK))) / 100.
		return esv	
		
	def _PTk2_Theta(self, Pressure, TempK):
		'''
		Calculate potential temperature [K]
		using pressure and temperature.

		Parameters::
		----------
		Pressure : float
			Pressure [hPa]
		TempK : float
			Temperature [deg K]
		'''

		Theta=TempK * ((1000. / Pressure)**(self.RD / self.CPD))
		return Theta 	


	def _RH_2_MixR(self, RH, Pressure, TempK, Opt='variable'):
		'''
		Calculate Mixing ratio [g/kg] (water vapor/dry air).  (D. Brunner)
	
		MixR = MassWater/MassDry = (Mw*e)/(Md*(p-e)) = Mw/Md * ((e*1000.)/(p-e))
		   Now because RH = e/esat*100, solving for e gives
						e = (RH/100.) * esat
		   Note: 1000. converts units to g/kg instead of g/g

		Parameters::
		----------
		RH : float
			Relative Humidity [%]
		Pressure
			Pressure [hPa]
		TempK
			Temperature [deg K]
		Opt : str
			Calculation option, possibilites are:
				'var' or 'variable' to use temperature dependence
				'water'
				'ice'
		'''
		e = (RH / 100.) * self._esat(TempK, Opt)
	
		# Compute mixing ratio, multiply by 1000 to convert to g/kg
		MixR = (self.MW / self.MD) * (e / (Pressure - e)) * 1000.
		return MixR
		
		
	def _Tk_RH_MixR_2_ThetaE(self,Pressure, TempK, RH, MixR, Opt='var', Calc=None):
		'''
		Calculate Equivalent Potential Temperature [K].
	
		AMS Glossary, Durran and Klemp 1982 modified.
		First calculate the dry air partial pressure,
		followed by the dry potential temperature.
		A discussion is provided in the above manuscript.

		Parameters::
		----------
		Pressure : float
			Pressure [hPa]
		TempK : float
			Temperature [deg K]
		RH : float
			Relative Humidity [%]
		MixR : float
			Mixing ratio [g/g]
		Opt : str
			Calculation option, possibilites are:
				'var' or 'variable' to use temperature dependence
				'water'
				'ice'
		Calc : int
			Definition to use for calculation, 
				1 = AMS Glossary/Paluch (1979)
				2 = Bolton (1980)
		'''
		# Choose the AMS glossary, Paluch 1979 as default
		if Calc is None:
			Calc = 1
		
		if (Calc == 1): 
			# Compute dry air partial pressure
			pDry = Pressure - self._esat(TempK, Opt) * RH / 100.
			# Compute dry potential temperature
			Term1 = TempK * (1000. / pDry)**(self.RD / self.CPD)
			Term3 = (self.LV * MixR) / (self.CPD * TempK)
			Term2 = RH**(-1 * self.RV * MixR / self.CPD)
		
			ThetaE = Term1 * Term2 * np.exp(Term3)
  
		if (Calc == 2): 
			Term1 = TempK * (1000. / Pressure)**((self.RD / self.CPD) * (1 - 0.28E-3 * MixR))
			Term2 = (1. / ((1. / (TempK - 55.))-(np.log(RH / 100.) / 2840.))) + 55.
			Term3 = np.exp((3.376 / Term2 - 0.00254) * (MixR * (1 + 0.81E-3 * MixR)))
	
			ThetaE=Term1 * Term3

		return ThetaE
		
	def _dewpoint_to_RH(self,T,Td):
	
		'''
		Calculate The RH given a dewpoint temperature
		
		Parameters::
		----------
	
		TempK : float
			Temperature [deg K]
		Td : float 
			Dewpoint temperature [deg k]
			
			
		Output::
		-------		
		RH : float
			Relative Humidity [%]
		
		'''
	
		esat = self._esat(T)
		
		
		e = (self.E0/100.)*np.exp((self.LV/self.RV)*((1.0/self.T0)-(1.0/Td)))
		
		RH = (e/esat)*100
		
		return RH
		
	###############

	def plot_dryadiabats(self):
		t0 = np.linspace(200,430,17)
		press = np.linspace(100,1000.)
	
		for temp in t0:
			theta = temp*(press/1000.)**(2./7.)
			
			self.ax1.semilogy((theta-273.15),press,color = '#7F4B10',linewidth = 0.3)
			
	def dry_lift(self,T,P,T_LCL,P_LCL):
		
		p =[]
		t = []
		
		t_lcl = T_LCL
		p_lcl = P_LCL
		
		p =P[np.arange(0,4)]
		
		t = T[np.arange(0,4)]
		
		t_avg = np.average(t)
		p_avg = np.average(p)
		
		print(t_avg,p_avg)
		
		parcel_press = np.linspace(p_avg,p_lcl,10)
		parcel_temp = np.linspace(t_avg,t_lcl,10)
		
		print(parcel_temp)
		print(parcel_press)
		
		print((parcel_press/1000.)**(self.RD/self.CPD))
		
		for temp in parcel_temp:
			
			theta = (temp+273)*(parcel_press/1000.)**(self.RD/self.CPD)
			
			
			
			#self.ax1.semilogy((theta-273.15),press,color = 'b-- ',linewidth = 0.3)
		
		
		
		
		#self.ax1.semilogy(t_avg,p_avg,'r^',ms = 10)
		
		
		return(theta -273.15,parcel_press)
		
		
		
		
		
		
			
		
			
			 			
'''
	def GammaW(tempk,pres,e=None):
		"""Function to calculate the moist adiabatic lapse rate (deg C/Pa) based
		on the temperature, pressure, and rh of the environment.

		INPUTS:
		tempk (K)
		pres (Pa)
		RH (%)

		RETURNS:
		GammaW: The moist adiabatic lapse rate (Dec C/Pa)
		"""

		tempc=tempk-273.15
		es=self._esat(tempk)
		ws=MixRatio(es,pres)

		if e is None:
		# assume saturated
		e=es

		w=MixRatio(e,pres)

		tempv=VirtualTempFromMixR(tempk,w)
		latent=Latentc(tempc)

		A=1.0+latent*ws/(Rs_da*tempk)
		B=1.0+Epsilon*latent*latent*ws/(Cp_da*Rs_da*tempk*tempk)
		Rho=pres/(Rs_da*tempv)
		Gamma=(A/B)/(Cp_da*Rho)
		return Gamma		
'''	