from matplotlib.ticker import ScalarFormatter, MultipleLocator
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
from thermocalcs import ThermoCalcs
from shearcalcs import ShearCalcs
from skew import SkewXTick

class DropSondes(object):

	'''
	Description:

	class to analyze dropsonde and sounding data
	#==========================
	requirements:

	thermodynaic calculations

	windshear calculations
	
	skew x axis class
	#==========================
	'''
	
	def __init__(self):

		#============================================#
		#initiate a new instance of thermocalcs class#
		#============================================#

		self.tC = ThermoCalcs()

		#===========================================#
		#initiate a new instance of shearcalcs class#
		#===========================================#

		self.sC = ShearCalcs()
		
		#Global height arrays for wind shear

		self.u_3km = []
		self.v_3km = []
		self.u_6km = [] 
		self.v_6km = []
		self.u_9km = []
		self.v_9km = []
		self.u_12km = []
		self.v_12km = []



	def get_sounding_data(self,filePath):

		'''
		method to retrieve standard sounding data from noaa radiosonde and U of Wyoming repository.

		Inputs:

		file path
		#============

		output:

		#Dictionary with variable for use in plotting

		'''

		data = filePath
		fp = open(data,'r')
		lines = fp.readlines()
		header = lines[1]
		fp.close()
		fp = open(data,'r')
		sounding_data = fp
		p,h,T,Td,RH,MIXR,wd,ws = np.genfromtxt(sounding_data, skip_header = 8, usecols= range(0,8) , dtype=float, delimiter = None, autostrip = True, missing_values='-9999.00', unpack=True, usemask=True)

		u = -ws*np.sin(np.radians(wd))
		v = -ws*np.cos(np.radians(wd))
		#RH = self.tC._dewpoint_to_RH(T+273.15,Td+273.15)

		mask = T.mask
		T=T[~mask]
		TD=Td[~mask]
		P=p[~mask]
		H = h[~mask]
		RH = RH[~mask]

		mask = u.mask
		U = u[~mask]
		V = v[~mask]

		data = dict()
		data['Header'] = header
		data['Temperature'] = T
		data['Dewpoint'] = Td
		data['Pressure'] = p
		data['Relative Humidity'] = RH
		data['U Component'] = U
		data['V Component'] = V
		data['Height'] = h
		data['Type'] = 'radioSonde'

		fp.close()





		#self.type = 'radioSonde'

		return data



	def get_dropsonde_data(self,filePath):

		'''
		method to retrieve dropsonde data

		Inputs:

		file path
		#============

		output:

		#Dictionary with variable for use in plotting

		'''

		data = filePath
		fp = open(data,'r')
		lines = fp.readlines()
		header = lines[1]
		fp.close()
		fp = open(data,'r')
		sounding_data = fp

		p,T,Td,RH,u,v,h = np.genfromtxt(sounding_data, skip_header= 4, usecols = (1,2,3,4,5,6,14), dtype = float, missing_values = '9999.0', unpack = True, usemask = True)

		#mask incoming T and dewpoint data

		T= np.ma.masked_greater_equal(T,999.0)
		Td = np.ma.masked_greater_equal(Td,999.0)
		RH = np.ma.masked_greater_equal(RH,999.0)
		h = np.ma.masked_greater_equal(h,99999.0)
		
		mask = T.mask
		T=T[~mask]
		TD=Td[~mask]
		P=p[~mask]
		H = h[~mask]
		RH = RH[~mask]
		
		mask = u.mask
		U = u[~mask]
		V = v[~mask]

		data = dict()
		data['Header'] = header
		data['Temperature'] = T
		data['Dewpoint'] = TD
		data['Pressure'] = P
		data['Relative Humidity'] = RH
		data['U Component'] = U
		data['V Component'] = V
		data['Height'] = H
		data['Type'] = 'dropsonde'

		return data
		
		

	def plot_skewtlogp(self,data,**kwargs):
	
		'''
		method to plot sounding or dropsonde data

		Inputs:
		
		dictionary of sounding data
		
		keyword arguments min, max, titles, and label.
	
		#============

		output:

		#Image

		'''
	
		T = data['Temperature']
		TD = data['Dewpoint']
		P = data['Pressure']
		
		
		
		
		self.fig = plt.figure(figsize=(10, 8))
		self.ax1 = self.fig.add_axes([0.05, 0.1, 0.6, 0.8],projection = 'skewx')
		plt.grid(True)
		
		self.plot_dryadiabats()
			
			
		self.ax1.semilogy(T, P, 'r-',linewidth =1.5)
		self.ax1.semilogy(TD, P, 'g-',linewidth =1.5)
		
			
		# Plot the data using normal plotting functions, in this case using
		# log scaling in Y, as dictated by the typical meteorological plot
		# Disables the log-formatting that comes with semilogy
		#set bounds for data on Skew T chart

		self.ax1.yaxis.set_major_formatter(ScalarFormatter())
		self.ax1.set_yticks(np.linspace(100,1000,10))
		self.ax1.set_ylim(1050,100)
		self.ax1.set_ylabel('Pressure mb')
		self.ax1.set_xlabel('Temperature C')
		self.ax1.xaxis.set_major_locator(MultipleLocator(10))
		self.ax1.set_title('Skew T log P')
		
		self.y_min = kwargs.get('y_min')
		self.y_max = kwargs.get('y_max')
		
		self.x_min = kwargs.get('x_min')
		self.x_max = kwargs.get('x_max')
		
		#data label, title and font sizes value assignments default values above .
		
		y_label = kwargs.get('y_label')
		x_label = kwargs.get('x_label')
		title = kwargs.get('title')
		font_size = kwargs.get('font_size')
		
		for item in ([self.ax1.xaxis.label,self.ax1.yaxis.label] + self.ax1.get_xticklabels()+self.ax1.get_yticklabels()):
			item.set_fontsize(font_size)

		#set the labels and titles defaults to 'none'

		self.ax1.set_ylabel(y_label)
		self.ax1.set_xlabel(x_label)
		self.ax1.set_title(title)
		
		
		
		self.ax1.set_ylim(self.y_max,self.y_min)
		self.ax1.set_xlim(self.x_min,self.x_max)
		

		
	def plot_hodograph(self,data):
	
		'''
		method to plot wind data on a hodograph

		Inputs:
		
		dictionary containing wind data. 

		file path
		#============

		output:
		
		#Image

		'''

		self.ax2 =self.fig.add_axes([.05, 0.6, 0.25, 0.3])
		
		#create axis and invert masks and assign values for U,V, and h coordinates
		
		U = data['U Component']
		V = data['V Component']
		H = data['Height']

		#=====================================================#
		# loop to place hodograph data into proer height array 
		#0-3km
		#3-6km
		#6-9km
		#9-12km
		#======#
		
		if data['Type']== 'radioSonde':

			for unew, vnew, hnew in zip(U,V,H):
				if hnew <= 3001.0:
					self.u_3km.append(unew)
					self.v_3km.append(vnew)
				elif hnew <= 6001.0:
					self.u_6km.append(unew)
					self.v_6km.append(vnew)
				elif hnew <= 9000.0:
					self.u_9km.append(unew)
					self.v_9km.append(vnew)
				elif hnew <= 12000.0:
					self.u_12km.append(unew)
					self.v_12km.append(vnew)    


		#print(len(self.u_3km))
		#print(len(self.v_3km))
		
		#=============================================================================#
		#take the first point of the next height level and attach it
		# to the last point of the previous line segment to connect hodograph segments
		#should only be done for soundings and non dense data
		#type determines how the data is placed in the array.
		#Soundings require endpoints be joined so that the hodograph is continuous.  
		#=============================================================================#
		
		if data['Type'] == 'radioSonde':

			self.u_3km.append(self.u_6km[0])
			self.u_6km.append(self.u_9km[0])
			self.u_9km.append(self.u_12km[0])

			self.v_3km.append(self.v_6km[0])
			self.v_6km.append(self.v_9km[0])
			self.v_9km.append(self.v_12km[0])

		#==============================#
		#Mask the arrays U and V coords
		#==============================#

		self.v_3km = np.ma.asarray(self.v_3km)
		self.v_6km = np.ma.asarray(self.v_6km)
		self.v_9km = np.ma.asarray(self.v_9km)
		self.v_12km = np.ma.asarray(self.v_12km)

		self.u_3km = np.ma.asarray(self.u_3km)
		self.u_6km = np.ma.asarray(self.u_6km)
		self.u_9km = np.ma.asarray(self.u_9km)
		self.u_12km = np.ma.asarray(self.u_12km)
		
		#plotting of the hodograph information
		#Different colors are sued to represent different height levels
		#r = 0-3km 
		#y = 3-6km
		#g = 6-9km
		#b = 9-12km
		
		self.ax2.plot(self.u_3km,self.v_3km,'r-',linewidth =3)  
		self.ax2.plot(self.u_6km,self.v_6km,'y-',linewidth =3,)
		self.ax2.plot(self.u_9km,self.v_9km,'g-',linewidth =3)          
		self.ax2.plot(self.u_12km,self.v_12km,'b-',linewidth =3)
	
		#set the default limits on the hodo axes
		# remove spines and set the axes to be invisible. 
		
		self.ax2.spines['left'].set_position('zero')
		self.ax2.spines['right'].set_color('none')
		self.ax2.spines['bottom'].set_position('zero')
		self.ax2.spines['top'].set_color('none')
		self.ax2.xaxis.set_ticks_position('bottom')
		self.ax2.yaxis.set_ticks_position('left')
		self.ax2.set_ylim(-50,50)
		self.ax2.set_xlim(-50,50)

		#Draw hodograph range rings
		#10 20 30 40 m/s circles

		circ1 = patches.Circle( (0, 0), 10,fc='white')
		circ2 = patches.Circle( (0, 0), 20,fc='white')
		circ3 = patches.Circle( (0, 0), 30,fc='white')
		circ4 = patches.Circle( (0, 0), 40,fc='white')

		#Add circles to plot as patches#
		#descedning order to make sure they layer ontop of each other.

		self.ax2.add_patch(circ4)
		self.ax2.add_patch(circ3)
		self.ax2.add_patch(circ2)
		self.ax2.add_patch(circ1)
	
	
	
	def plot_aux_graph(self,x_value,y_value,**kwargs):
		
		'''
		method to plot an auxiliary graph of user defined data

		Inputs:

		# sounding data (User specified)
		#kwargs to adjust plot dimensions scales label, and title.
		#============

		output:

		#Image

		'''
		
		#define axes fro figure.
		#set grid and plot user specified data
	
		self.ax3 = self.fig.add_axes([.7,0.7,.29,.24])
		plt.grid(True)
		self.ax3.plot(x_value,y_value)
		
		#Max and min graph bonds value assignment. Defaults to autoscale.
		
		y_min = kwargs.get('y_min')
		y_max = kwargs.get('y_max')
		
		x_min = kwargs.get('x_min')
		x_max = kwargs.get('x_max')
		
		#data label, title and font size value assignments .
		
		y_lable = kwargs.get('y_lable')
		x_lable = kwargs.get('x_lable')
		title = kwargs.get('title')
		font_size = kwargs.get('font_size')
		
		#keyword arguments for minimum and maximum values
		
		self.ax3.set_ylim(y_min,y_max)
		self.ax3.set_xlim(x_min,x_max)
		
		#set label and tick fontsizes to user defined value. Defaults to 10 pt. 
		
		for item in ([self.ax3.xaxis.label,self.ax3.yaxis.label] + self.ax3.get_xticklabels()+self.ax3.get_yticklabels()):
			item.set_fontsize(font_size)

		#set the labels and title from keyword arguments. Defaults to 'none'

		self.ax3.set_ylabel(y_lable)
		self.ax3.set_xlabel(x_lable)
		self.ax3.set_title(title)
		
		
		
	def generate_parameter_list(self):
	
		'''
		method to generate a list of parameters calculated from the sounding data

		Inputs:

		#None
		#============

		output:

		#Blank axis for plotting
		'''
		
		#Set axes position and set both axes invisible
		
		self.ax4 = self.fig.add_axes([.65,0.1,.35,.54])
		self.ax4.xaxis.set_visible(False)
		self.ax4.yaxis.set_visible(False)
		
		
	
	def run_thermo_calcs(self,data):
		
		'''
		method to calculate thermodynamic parameters from dropsonde data

		Inputs:

		#Dictionary of sounding data.
		#Uses calculations from thermocalcs.py 
		#============

		output:

		#multiple arrays containing thermodynamic information

		'''
	
		T = data['Temperature']
		Td = data['Dewpoint']
		p = data['Pressure']
		RH = data['Relative Humidity']
		u = data['U Component']
		v = data['V Component'] 
		h = data['Height']

		self.LCLT = round((self.tC._LCL_Temperature(h,T+273.15,Td+273.15)-273.15),2)
		self.LCLP = round((self.tC._LCL_Pressure(h, p, T+273.15,Td+273.15)),0)
		self.LCLZ = round((self.tC._LCL_Height(h,p,T+273.15,Td+273.15)),0)
		self.THETA = self.tC._PTk2_Theta(p, T+273.15)
		self.MIXR = self.tC._RH_2_MixR(RH, p, T+273.15)
		self.THETAE = self.tC._Tk_RH_MixR_2_ThetaE(p, T+273.15, RH, self.MIXR/1000.)
		self.ESAT = self.tC._esat(T+273.15)	

	def plot_thermo_calcs(self):

		'''
		method to plot the thermodynamic parameters on the parameter list. 

		Inputs:

		#None
		#============

		output:

		#Image
		'''


		#plot the parameters on the list generated. 
		self.ax4.text(.01,.01,'LCL Pressure: '+ str(self.LCLP)+(' hPa'))
		self.ax4.text(.01,.04,'LCL Temp: '+ str(self.LCLT)+' c')
		self.ax4.text(.01,.07,'LCL Height: '+ str(self.LCLZ)+ ' m')


	def plot_dryadiabats(self,**kwargs):

		'''
		method to plot the dry adibats. Used in the plotskewtlogp method. 

		Inputs:

		#kwargs (does not function)
		#============

		output:

		#Image

		'''
		#test = self.shear1km

		#temperature array and pressure array

		t0 = np.linspace(200,430,17)
		press = np.linspace(100,1000.)

		#retrieve desired line style from user kwargs

		line_style = kwargs.get('line_style')
		if line_style == None:
			line_style = '-' 	
	
		#loop to calculate using poissions equation the dry adiabats for the skew t.	

		for temp in t0:
			theta = temp*(press/1000.)**(2./7.)
	
			#plot the dry adiabats given a specified color
	
			self.ax1.semilogy((theta-273.15),press,line_style,color = '#7F4B10',linewidth = 0.5)
			


	def plot_wind_barbs(self,data,**kwargs):

		P = data['Pressure']
		U = data['U Component']
		V = data['V Component'] 
		H = data['Height']

	




		#Copy y axis to plot wind barbs
		#Set x lim for windpbarbs
		#adjust location of barbs from x =0
		#plot every 30th windbarb 
		#set axis to invisible
		self.ax1_copy = self.ax1.twiny()
		self.ax1_copy.set_xlim(self.x_min,self.x_max)
		self.ax1_copy.set_ylim(self.y_min,self.y_max)
		x_const = np.zeros(P.shape) +(self.x_max - 2)
		self.ax1_copy.xaxis.set_visible(False)
		

		if data['Type'] == 'radioSonde':
		
			mask = U.mask
			U = U[~mask]
			V = V[~mask]
			P = P[~mask]
			
			self.ax1_copy.barbs(x_const[::3],P[::3], U[::3], V[::3])
		
		else:
		
		
			self.ax1_copy.barbs(x_const[::40],P[::40], U[::40], V[::40])
	
	
	def run_shear_calcs(self,data):

		'''
		method to calculate thermodynamic parameters from dropsonde data

		Inputs:

		#Dictionary of sounding data.
		#Uses calculations from thermocalcs.py 
		#============

		output:

		#multiple arrays containing thermodynamic information

		'''

		T = data['Temperature']
		Td = data['Dewpoint']
		p = data['Pressure']
		RH = data['Relative Humidity']
		u = data['U Component']
		v = data['V Component'] 
		h = data['Height']	
		
		mask = h.mask
		u = u[~mask]
		v = v[~mask]
		h = h[~mask]
		
		

		self.SHEAR1KM = self.sC._VertShear_Sfc_to_1km(h,u,v)
		self.SHEAR3KM = self.sC._VertShear_Sfc_to_3km(h,u,v)
		self.SHEAR6KM = self.sC._VertShear_Sfc_to_6km(h,u,v)
		self.BULKSHEAR1km = round(self.sC._bulkshear_sfc_1km(h,u,v),2)
		self.BULKSHEAR3km = round(self.sC._bulkshear_sfc_3km(h,u,v),2)
		self.BULKSHEAR6km = round(self.sC._bulkshear_sfc_6km(h,u,v),2)
		#self.ax4.text(.01,.1,'0-1 km shear: '+ str(self.SHEAR1KM)+(' 1/s'))




	def plot_shear_calcs(self):

		'''
		method to plot the thermodynamic parameters on the parameter list. 

		Inputs:

		#None
		#============

		output:

		#Image
		'''
		
		#bitch = self.SHEAR1KM


		#plot the parameters on the list generated. 
		self.ax4.text(.01,.1,'0-1 km shear: '+ str(self.SHEAR1KM)+(' 1/s'))
		self.ax4.text(.01,.14,'0-3 km shear: '+ str(self.SHEAR3KM)+' 1/s')
		self.ax4.text(.01,.18,'0-6 km shear: '+ str(self.SHEAR6KM)+ ' 1/s')
		self.ax4.text(.01,.22,'0-1km Bulk Shear: '+str(self.BULKSHEAR1km)+ ' m/s')
		self.ax4.text(.01,.26,'0-3km Bulk Shear: '+str(self.BULKSHEAR3km)+ ' m/s')
		self.ax4.text(.01,.3,'0-6km Bulk Shear: '+str(self.BULKSHEAR6km)+ ' m/s')

	def dry_lift(self,data):

		T = data['Temperature']
		Td = data['Dewpoint']
		p = data['Pressure']
		RH = data['Relative Humidity']
		u = data['U Component']
		v = data['V Component'] 
		h = data['Height']


		t_parcel,p_parcel = self.tC.dry_lift(T,p,self.LCLT,self.LCLP)
		
		#print(t_parcel,p_parcel)
		
		self.ax1.semilogy(t_parcel,p_parcel,'k--',ms =1)
		
