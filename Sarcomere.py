#Sarcomere class to simulate sliding filament theory
#June 13, 2018

#import
import csv
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.axis as ax
import matplotlib.path as mpath
import scipy.stats as stats

#radius = 1.0
#Lm = 10.0

#Set a class to read csv files 
# Data file must only have lengths of sarcomeres and Tension in N/cm2
class Sarcomere_Data:
	def __init__(self, filename=None):
		self.headers=[]
		self.header2col={}
		self.L_thick = None
		self.L_thin = None
		self.L_bare = None
		self.L_bt = None
		self.L_tt= None

		#Tensions of Sliding Filament Theory
		self.D_tension = None
		self.P_tension = None
		self.A_tension = None

		self.sft_tension = None

		#Data file information
		self.data=None

		#Converstion from kg/cm2 to N/cm2 if needed
		self.TN = None

		#Converstion to cm from um
		self.cmLS = None
		self.umLs = None

		#Conversion to Newtons
		self.Force = None

		#Combined data to 2d array
		self.sarc = None
		self.sarcT=None
		self.Hmark = None

		#Partitioned Data with Force
		self.DescF = None
		self.PlatF = None
		self.AscF = None

		#Work calculation parameters
		self.ST_work = None

		#Energy Density
		self.ED = None

		#crucial points for sliding filament model
		self.marker = None

		#if filename is not None, call self.read(filename)
		if filename !=None:
			self.read(filename)

	def get_headers(self):
		return self.headers

	def get_num_dimensions(self):
		return len(self.headers)

	def get_num_points(self):
		return self.data.size

	def read(self, filename):

		#Reading the file
		with open(filename, 'rU') as fp:
			csv_reader= csv.reader(fp)

			#making a list of the headers
			header= next(csv_reader)
			#print(header)

			#making a list of types
			typelist=next(csv_reader)
			#print(typelist)
			n_count=0
			for i in range(len(typelist)):
				if typelist[i]== 'numeric':
					self.headers.append(header[i])
					self.header2col[header[i]]=n_count
					n_count+=1

			print(self.header2col)

		#update to store a non-numeric column separately as the string it was read
		#numpy matrix of the file data
			x=list(csv_reader)
			self.data=np.matrix(np.zeros((len(x), len(self.headers))))
			for ridx in range(len(x)):
				for cidx in range(len(self.headers)):
					self.data[ridx, cidx]=x[ridx][cidx]
			print(self.data)

	def colHeader(self, headers):
		mat=np.matrix(np.zeros((self.data.shape[0], len(headers)) ))
		for i in range(len(headers)):
			mat[:,i]=self.data[:, self.header2col[headers[i]]]
		return mat

	#method to calculate the length of the sarcomere for the simulation
	def calcL(self, L_thick, L_thin, L_bare, Zdisc):
		self.L_thick = L_thick
		self.L_thin = L_thin
		self.L_bare = L_bare
		self.Z_disc = Zdisc
		self.L_tz = L_thin + Zdisc

		print('\n', "These are the Lengths of the filaments: ")
		print("Myosin(um) = ", self.L_thick)
		print("Actin (um) = ", self.L_thin)
		print("Actin with Z disc (um) = ", self.L_tz)
		print("Bare zone (um) = ", self.L_bare)

		self.L_tt = self.L_thick + self.L_thin + self.Z_disc
		self.L_bt = self.L_bare + self.L_thin + self.Z_disc
		print('\n', "These are the sum of lengths of filaments: ")
		print("Sum of Thick and Thin Filaments = ", self.L_tt)
		print("Sum of Bare zone and Thin Filaments = ", self.L_bt)

		return self.L_thick, self.L_thin, self.Z_disc, self.L_bare, self.L_tz, self.L_tt, self.L_bt

	#Method to calculate the tension of Sliding Filament Theory for R2 value
	def calcT(self, length):

		#Slope Calculation
		sub = (self.L_tt )-(self.L_bt)
		Desc_slope= -1 / sub
		#print("Descending Limb Slope = ", Desc_slope, '\n')

		Asc_slope= 1/sub
		#print("Ascending Limb Slope= ", Asc_slope, '\n')

		# b calculation for the Tensions
		Desc_b = (Desc_slope*(-1))*self.L_tt
		#print("B for Descending Limb = ", Desc_b, '\n')

		Asc_b= 1 - (Asc_slope)*self.L_tz		#self.L_bt
		#print("B for Ascending Limb = ", Asc_b, '\n')

		#Calculating the Tensions based on provided length
		if length >= self.L_tt:
			#print('\n', "Tension greater than thick + thin filaments = ", 0)
			return 0

		elif length >= self.L_bt:
			self.D_tension= (Desc_slope*length) + Desc_b
			#print('\n', "Descending Limb %T0= ", self.D_tension, '\n')
			return self.D_tension

		elif length >= self.L_tz and length <= self.L_bt:
			#print('\n', "Plateau %T0= 1.00", '\n')
			self.P_tension = np.amax(self.data[:,1])
			return 1.0

		elif length <= self.L_tz:
			self.A_tension= (Asc_slope*length)+ Asc_b

			if self.A_tension <= 0.0:
				self.A_tension = 0.0

			#print('\n',"Ascending Limb %T0 = ", self.A_tension, '\n')
			return self.A_tension

		else:
			return 0

	# Method to get the 4 crucial marker points 
	def getMark(self, medTen):
		max_LT=np.hstack((self.L_tt, self.calcT(self.L_tt) ))
		#print('\n', "Marker point 1: ", max_LT)

		Rplat_LT= np.hstack((self.L_bt, self.calcTen(self.calcT(self.L_bt), medTen) ))
		#print("Marker point 2: ", self.calcT(self.L_bt)


		Lplat_LT = np.hstack((self.L_tz, self.calcTen(self.calcT(self.L_tz), medTen) )) 
		#print("Marker point 3: ", Lplat_LT)

		sub = (self.L_tt )-(self.L_bt)
		#print(sub)
		asc_slope = (Rplat_LT[1]-max_LT[1])/sub
		asc_b = 1 - (asc_slope)*self.L_tz

		x = Lplat_LT[0] - (Lplat_LT[1]/asc_slope)

		min_LT = np.hstack((x, 0.0 ))
		#print("Marker point 4: ", min_LT)

		self.marker = np.vstack((max_LT, Rplat_LT, Lplat_LT, min_LT))
		print('\n', "Marker Points: ", '\n', self.marker)

		return self.marker

	#method to convert the empirical data's Tension in newtons
	# tension parameter = array of normalized tension
	# medT parameter = median max Tension from Plateau [ kg/cm2 ]
	def calcTen(self, tension, medT):
		#PARAMETERS = 
		Ncm = 9.81
		#print('\n',"Newtons per kg =",Ncm)

		self.medT = medT
		medTen = self.medT * Ncm
		#print('\n', "Median Max Tension in Newtons = ", medTen )
		self.tension = tension

		#undo the normalization of the data
		self.TN = self.tension * medTen

		#print('\n', "Tension (N/cm2) = ", '\n', self.TN)

		return self.TN

	#method to convert the lengths from micrometers to cm
	def calcCM(self):
		self.umLs = self.data[:,0]
		#print('\n', "Sarcomere Lengths in um = ", '\n', self.umLs)

		self.cmLs = self.umLs / 10000
		#print('\n', "Sarcomere Lengths in cm = ", '\n', self.cmLs)

		return self.cmLs, self.umLs

	def calcR(self, medT):
		fig = plt.figure()
		#Partitioning the empirical data:
		isDesc = []
		isPlat= []
		isAsc= []
		for idx in self.data[:,0]:
			if idx >= (self.L_bt):
				isDesc.append(idx)
				self.isDesc= np.concatenate(isDesc)

			elif idx >= (self.L_tz) and idx <= (self.L_bt):
				isPlat.append(idx)
				self.isPlat= np.concatenate(isPlat)

			elif idx <= self.L_tz:
				isAsc.append(idx)
				self.isAsc= np.concatenate(isAsc)

		#print('\n', "Desc = ", '\n', self.isDesc)
		#print('\n', "Plat = ", '\n', self.isPlat)
		#print('\n', "Asc = ", '\n', self.isAsc)

		#To partition the actual data
		Desc=[]
		Plat=[]
		Asc=[]
		for mdx in range(self.isDesc.shape[0]):
			Desc.append(self.data[mdx])
			self.Desc= np.vstack(Desc)
		print('\n', "Descending limb coordinates", '\n', self.Desc)

		for jdx in range(self.isPlat.shape[0]):
		 	Plat.append(self.data[mdx+jdx+1])
		 	self.Plat=np.vstack(Plat)
		print('\n', "Plateau coordinates", '\n', self.Plat)

		for kdx in range(self.isAsc.shape[0]):
		 	Asc.append(self.data[mdx+jdx+kdx+2])
		 	self.Asc=np.vstack(Asc)
		print('\n', "Ascending limb coordinates", '\n', self.Asc)

		#Calculate the R_squared of Descending Limb
		#D_mean = 
		D_tension=[]
		Des_ten=[]
		Desc_SSres=0.0
		for i in range(self.isDesc.shape[0]):
		 	f = self.calcTen(self.calcT(self.Desc[i,0]), medT)
		 	D_tension.append(f)
		 	self.D_ten = np.vstack(D_tension)
		 	Desc_SSres+= (self.Desc[i,1]-f)**2
		 	#plt.plot(self.Desc[i,0], f, 'ro' )
		D_mean = np.mean(self.data[:,1])
		Desc_SStot = np.sum(np.square(self.Desc[:,1] - D_mean))
		#print("This is the array of f values: ", '\n', self.D_ten)
		#print("Descending Limb SSres = ", Desc_SSres)
		#print("This is the mean of Descending Limb Tension = ", D_mean)
		#print("This is the SStot for Descending Limb = ", Desc_SStot)

		Desc_r2 = 1 - (Desc_SSres / Desc_SStot)
		print('\n', "R-squared value for Descending Limb = ", Desc_r2)

		# #Calculate the R_squared of Plateau
		P_tension=[]
		P_SSres=0.0
		for j in range(self.isPlat.shape[0]):
		 	fi = self.calcTen(self.calcT(self.Plat[j,0]), medT)
		 	P_tension.append(fi)
		 	self.P_ten = np.vstack(P_tension)
		 	P_SSres += (self.Plat[j,1]-fi)**2
		 	#plt.plot(self.Plat[j,0], fi, 'yo' )
		P_mean = np.mean(self.data[:,1])
		P_SStot = np.sum(np.square(self.Plat[:,1] - P_mean))
		#print('\n', "This is the array of f values: ", '\n', self.P_ten)
		#print("SSres = ", P_SSres)
		#print("mean of Plateau tension = ", P_mean)
		#print("SStot = ", P_SStot)

		Plat_r2 = 1 - (P_SSres / P_SStot)
		print("R-squared value for Plateau = ", Plat_r2)

		#Calculate the R_squared of Ascending Limb
		A_tension = []
		A_SSres=0.0

		for n in range(self.isAsc.shape[0]):
		 	fj = self.calcTen(self.calcT(self.Asc[n,0]), medT)
		 	A_tension.append(fj)
		 	self.A_ten = np.vstack(A_tension)
		 	A_SSres += (self.Asc[n,1]-fj)**2
		 	#plt.plot(self.Asc[n,0], fj, 'go' )
		A_mean = np.mean(self.data[:,1])
		A_SStot= np.sum(np.square(self.Asc[:,1] - A_mean))
		#print('\n', "This is the array of f values: ", '\n', self.A_ten)
		#print("SSres = ", A_SSres)
		#print("mean of Ascending Limb = ", A_mean)
		#print("SStot= ", A_SStot)

		Asc_r2= 1-(A_SSres/A_SStot)
		print("R-squared value for Ascending Limb = ", Asc_r2)

		#Total R-squared calculation
		y_bar= np.mean(self.data[:,1])
		SStot= np.sum(np.square(self.data[:,1]-y_bar))
		SSres=0.0
		for sdx in range(self.data.shape[0]):
		 	fk=self.calcTen(self.calcT(self.data[sdx,0]), medT)
		 	SSres += (self.data[sdx,1]-fk)**2
		 	#plt.plot(self.data[sdx,0], fk, 'ko' )
		self.R2= 1-(SSres/SStot)
		print("Overall R-squared = ", self.R2)

		#plt.show()

		# #print("Lm = ", Lm)

		#Combine theoretical tensions of lengths
		self.sft_tension = np.vstack((self.D_ten, self.P_ten, self.A_ten))
		print('\n', "Sliding Filament Theoretical Tension Values= ", '\n', self.sft_tension)

		return self.R2, Desc_r2, Plat_r2, Asc_r2, self.sft_tension

	#method to calculate the Force of Sarcomere
	# takes in radius of cylinder in cm 
	# returns Force ( N )
	def sarcF(self, tension, radius):
		self.radius = radius
		print('\n', "The radius of this model = ", self.radius)
		tension = tension

		#Cross Sectional Area calculations
		self.CSA = (self.radius**2) * np.pi
		#print('\n', "Cross Section Area (cm2) = ", self.CSA)

		#Force in Newtons ( N )
		self.Force = self.CSA * tension
		print('\n', "Force (N) = ", '\n', self.Force)

		return self.Force

	#method to combine cmLs and sarcF 
	# should make a 2D array of [cm , N] and another array of [ cm , N/cm2]
	def Homsarc(self):

		#numpy array of [ cm , N ]
		self.sarc = np.hstack((self.cmLs, self.Force))
		print('\n', "Combined Data for Homogeneous Muscle Sarcomere: ", '\n', self.sarc)

		#numpy array of [ cm, N/cm2 ]
		self.sarcT = np.hstack((self.cmLs, self.data[:,1]))
		print('\n', "Sarcomere N/cm2 Data: ", '\n', self.sarcT)

		return self.sarc, self.sarcT

	#method to partition the homogeneous sarcomere data
	# data parameter = partitioning the column array of either Force or Tension array from Homsarc
	def partHS(self, data):
		#Parameters:
		x = 10000
		self.Lb = self.L_bare / x
		self.Ltz= (self.L_thin+self.Z_disc) / x
		self.Lthick = self.L_thick / x
		sarc = data

		isDesc=[]
		isPlat=[]
		isAsc=[]
		for idx in sarc[:,0]:
			if idx >= (self.Lb + self.Ltz):
				isDesc.append(idx)
				self.isDesc = np.concatenate(isDesc)

			elif idx >= self.Ltz and idx <= (self.Lb + self.Ltz):
				isPlat.append(idx)
				self.isPlat = np.concatenate(isPlat)

			elif idx <= self.Ltz:
				isAsc.append(idx)
				self.isAsc = np.concatenate(isAsc)

		#print("Testing Descending Limb = " , self.isDesc)

		#To partitioin the actual data
		Desc=[]
		Plat=[]
		Asc=[]
		for mdx in range(self.isDesc.shape[0]):
			Desc.append(sarc[mdx])
			self.DescF = np.vstack(Desc)
		print('\n', "Descending limb coordinates: ", '\n', self.DescF)

		for jdx in range(self.isPlat.shape[0]):
			Plat.append(sarc[mdx+jdx+1])
			self.PlatF = np.vstack(Plat)
		print('\n', "Plateau coordinates: ", '\n', self.PlatF)

		for kdx in range(self.isAsc.shape[0]):
			Asc.append(sarc[mdx+jdx+kdx+2])
			self.AscF = np.vstack(Asc)
		print('\n', "Ascending limb coordinates: ", '\n', self.AscF)
		return self.DescF, self.PlatF, self.AscF

	#method to get the number of sarcomeres in series
	def numSarc(self, Lm):
		#average of markers in plateau
		A = self.marker[1,0]
		B = self.marker[2,0]

		### Middle of the Plateau for L0 ###
		L0 = ((A + B) / 2) / 10000
		self.Lm = Lm  

		self.num_sarcomere = self.Lm / L0
		self.halfSarc = self.num_sarcomere * 2

		print('\n', "Number of Whole Sarcomeres = ", self.num_sarcomere)
		print("Number of Half Sarcomeres = ", self.halfSarc)

		return self.halfSarc, self.num_sarcomere

	#method to calculate the work for the sliding filament theory --> kept in um
	def calcWork(self, radius):
		#print(self.marker)
		#don't use markers for the work! Use actual points from sliding filament theory
		#
		#5.5 J (N*m)

		#Units in [cm, N]
	
		max_Ld = self.marker[0,0] / 10000
		min_Ld = self.marker[1,0] / 10000
		DL = max_Ld - min_Ld
		max_DT = self.sarcF(self.marker[1,1], radius)

		max_Lp = self.marker[1,0] / 10000
		min_Lp = self.marker[2,0] / 10000
		PL = max_Lp - min_Lp
		max_PT = self.sarcF(self.marker[2,1], radius)

		max_La = self.marker[2,0] / 10000
		min_La = self.marker[3,0] / 10000
		AL = max_La - min_La
		max_AT = self.sarcF(self.marker[2,1], radius)

		#print('\n', "Min of Descending Limb = ", min_Ld, '\n', "Max of Descending Limb = ", max_Ld)
		#print('\n', "Min of Plateau = ", min_Lp, '\n', "Max of Plateau = ", max_Lp)
		#print('\n', "Min of Ascending Limb = ", min_La, '\n', "Max of Ascending Limb = ", max_La)

		self.Desc_work = (0.5 * ((max_DT*DL) / 100)) * self.num_sarcomere
		print('\n', "Work of Descending Limb = ", self.Desc_work)

		self.Plat_work = ((max_PT * PL) / 100 ) * self.num_sarcomere
		print("Work of Plateau region = ", self.Plat_work)

		self.Asc_work = (0.5 * ((max_AT * AL) / 100 )) * self.num_sarcomere
		print("Work of Ascending Limb = ", self.Asc_work)

		self.ST_work = (self.Desc_work + self.Plat_work + self.Asc_work)
		print('\n',"Calculating the Work for the Sliding Filament Theory:" , self.ST_work)

		return self.ST_work

	#method to calculate the energy density
	def calcED(self, radius, Lm, work):
		self.radius = radius
		self.Lm = Lm
		self.work = work

		#calculate the volume
		self.muscle_vol = (self.radius**2)*self.Lm*(np.pi)
		#print('\n', "Volume of the homogeneous muscle = ", self.muscle_vol)

		#calculate the energy density of Sliding Filament Theory
		self.ED = self.work / self.muscle_vol
		print("Energy Density = ", self.ED)

		return self.ED

	#method to plot the overall Force and total sum of sarcomere length
	def HsPlot(self):
		Ls = self.sarc[:,0] 
		#print(Ls)
		F = self.sarc[:,1]
		#print(F)

		#x = (self.marker[:,0] / 10000)		#converting the marker to cm from um
		#print('\n', x)
		#markTen = self.calcTen(self.marker[:,1], 2.6)
		#b = self.calcTen(self.marker[:,1], 2.6)
		#y = self.sarcF(sarc.data[:,1], 1.0)
		#print('\n', y)

		plt.plot(Ls, F, '^')

		plt.xlabel('Length (cm)', fontsize= 16)
		plt.ylabel('Force (N)', fontsize = 16)
		plt.title('Sarcomere Length and Force Curve', fontsize=20)

		plt.tick_params(axis='both', which='major', labelsize=14)

		#A = self.marker[:,0] / 10000
		#print(A)
		#B = y
		#print(B)
		#plt.plot(x, y, 'rs')
		#plt.plot(A,B, 'r--')
		plt.show()

	#method to plot the sarcomere and Tension curve
	def sarcPlot(self):
		L_sarc = self.data[:,0]
		#print('\n', L_sarc)

		Ten = self.data[:,1]
		#print('\n', "Tension in N/cm2: ", '\n', Ten)

		x = self.marker[:,0]
		print('\n', x)

		y = self.marker[:,1]

		#plt.plot(L_sarc, Ten, 'o')
		plt.xlabel('Length (um)', fontsize = 16)
		plt.ylabel('Tension (N/cm2)', fontsize = 16)
		plt.title('Isometric Tetani and Length Curve', fontsize = 20)

		plt.tick_params(axis = 'both', which = 'major', labelsize = 14)

		plt.plot(x,y, 'mp', markersize = 10)
		plt.plot(x,y, 'm--')
		plt.show()

	#method to just plot the empirical data
	def plot(self):
		L=self.data[:,0]
		#print(L)
		T=self.data[:,1]
		#print(T)

		#x = self.marker[:,0]
		#print(x)
		#y= self.marker[:,1]
		#print(y)

		#AE = self.marker [0:5]
		#print(AE)
		#A=AE[:,0]
		#E=AE[:,1]

		#Plotting points of Gordon's data
		plt.plot(L,T, 'o')
		plt.xlabel('Length (um)', fontsize=16)
		plt.ylabel('Tension (%T0)', fontsize=16)
		plt.title('Zachar Isometric Tetani Length and Tension Curve', fontsize=18)

		#ticks
		plt.tick_params(axis='both', which='major', labelsize=14)

		#Plotting crucial points of sliding filament model
		#plt.plot(x, y, 'yo', markersize = 10)
		#plt.plot(A,E, 'y--')
		plt.show()

sarc=Sarcomere_Data("Isometric Tetani.csv")
## Input thick, thin, bare, Zdisc lengths in this order!
sarc.calcL(1.6, 2.0, 0.2, 0.05)

### Plug in the length of interest ###
sarc.calcT(1.5)
sarc.getMark(2.6)

### Plug in median length of sarcomere ###
sarc.calcR(2.6)
sarc.calcCM()

### Plug in tension and radius of muscle model ###
sarc.sarcF(sarc.data[:,1],1.0)

sarc.Homsarc()
sarc.partHS(sarc.sarc)
#sarc.plot()

### Plug in the length of the muscle model ###
sarc.numSarc(10.0)

#sarc.HsPlot()
sarc.sarcPlot()

### Work calculation ###
sarc.calcWork(1.0)

### Plug in: radius, length of muscle model, and work value ### 
sarc.calcED(1.0, 10.0, sarc.ST_work)