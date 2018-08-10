### Sarcomere class that simulates just by plugging numbers in ###
# July 29, 2018 #

### import ###
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm

### Set a class to simulate the numbers ###
class Sarcomere_Simulation:
	def __init__(self, thick=None, thin=None, bare=None, zdisc=None, T0=None, radius=None, LM=None):

		print('\n', "INPUT VALUES ACCORDINGLY", '\n')

		### Inputting values in the terminal ###

		if thick == None:
			thick= float(input("Thick filament/Myosin (um) = "))
		else:
			print("List of Thick Filament/Myosin (um) = ", thick)

		if thin == None:
			thin = float(input("Thin Filament/Actin (um) = "))
		else:
			print("List of Thin Filament/Actin (um) = ", thin)

		if bare ==None:
			bare = float(input("Bare Zone (um) = "))
		else:
			print("Bare zone (um) = ", bare)

		if zdisc == None:
			zdisc = float(input("Z Disc (um) = "))
		else:
			print("Z Disc (um) = ", zdisc)

		if T0 == None:
			T0 = float(input("Median Max Tension (N/cm2) = "))
		else: 
			print("Median Max Tension (N/cm2) = ", T0)

		if radius == None:
			radius = float(input("Radius (cm) = "))
		else: 
			print("Radius of Muscle Model (cm) = ", radius)

		if LM == None:
			LM = float(input("Muscle Length (cm) = "))
		else:
			print("Length of Muscle Model (cm) = ", LM)

		### Parameters for Sarcomere Lengths and Median Max Tension (T0) ###
		self.L_thick = thick 	# Thick Filament/Myosin Length
		self.L_thin = thin 		# Thin Filament/Actin Length
		self.L_bare = bare 		# Bare Zone Length
		self.L_zdisc = zdisc 	# Z Disc Length
		self.T0 = T0 			# Median Max Tension 
		self.radius = radius 	# Radius of the Muscle Model
		self.L_muscle = LM 		# Lenght of the Muscle Model

		#sum of filament lengths
		self.L_tt = None 		# Thick and Thin Filaments and Z disc
		self.L_bt = None 		# Thin Filament, Bare zone and Z disc
		self.L_tz = None 		# Thin Filament and Z disc

		#slope and y-intercept values
		self.Desc_slope = None 	# Slope of the Descending Limb
		self.Desc_y_int = None 	# Y-intercept of the Descending Limb
		self.Asc_slope = None 	# Slope of the Ascending Limb
		self.Asc_y_int = None 	# Y-intercept of the Ascending Limb

		# Theoretical Tension values
		self.Desc_Tension = None 	# Theoretical Tension of Descending Limb
		self.Asc_Tension = None 	# Theoretical Tension of Ascending Limb

		#Characteristic Points
		self.marker = None 

		#Force calculation
		self.Force = None

		#Number of sarcomeres
		self.num_sarcomere = None 	# Number of Whole Sarcomeres in series
		self.halfSarc = None 		# Number of half sarcomeres in series

		#Work calculation
		self.Desc_work = None 	# Calculated work of Descending Limb
		self.Plat_work = None 	# Calculated work of Plateau region
		self.Asc_work = None 	# Calculated work of Ascending Limb
		self.Work = None 		# Calculated Work Total

		#Energy Density
		self.ED = None 			# Calculated Energy Density

	def calcLength(self):
		## Printing out the individual lengths ##
		print('\n', "These are the Filament Lengths: ")
		print("Myosin (um) = ", self.L_thick)
		print("Actin (um) = ", self.L_thin)
		print("Bare Zone (um) = ", self.L_bare)
		print("Z Disc (um) = ", self.L_zdisc)

		print('\n', "Median Max Tension = ", self.T0, " N/cm2 ")

		print('\n', "----------------------------------------", '\n')

		print("The sum of the filament lengths: ")
		
		# Length of thick, thin and z disc combined #
		self.L_tt = self.L_thick + self.L_thin + self.L_zdisc
		print("Sum of Thick, Thin and Zdisc = ", self.L_tt, "um")

		#Length of Bare, Thin and Z disc combined #
		self.L_bt = self.L_thin + self.L_bare + self.L_zdisc
		print("Sum of Thin Filament, Bare Zone and Zdisc = ", self.L_bt,
			"um")

		#Length of Thin Filament and Z disc combined #
		self.L_tz = self.L_thin + self.L_zdisc
		print("Sum of Thin Filament and Zdisc = ", self.L_tz, "um")

		return self.L_tt, self.L_bt, self.L_tz

	### Calculate the slope and y-intercept of the three parts 
	### of the length-tension curve ###
	def calcLinear(self):

		print('\n', "----------------------------------------", '\n')

		xsub = self.L_tt - self.L_bt

		# Slope of the Descending Limb
		self.Desc_slope = -self.T0 / xsub
		print("Slope of the Descending Limb = ", self.Desc_slope)

		# Y-intercept of the Descending Limb
		self.Desc_y_int = (self.Desc_slope*(-1))*self.L_tt
		print("Y-intercept of the Descending Limb = ", self.Desc_y_int)

		# Slope of the Ascending Limb
		self.Asc_slope = self.T0 / xsub
		print('\n', "Slope of the Ascending Limb = ", self.Asc_slope)

		#Y-intercept of the Ascending Limb
		self.Asc_y_int = self.T0 - (self.Asc_slope)*(self.L_tz)
		print("Y-intercept of the Ascending Limb = ", self.Asc_y_int)

		print('\n', "----------------------------------------", '\n')

		return self.Desc_slope, self.Asc_slope, self.Desc_y_int, self.Asc_y_int

	### Calculating the Theoretical Tension in Newtons ###
	def calcT(self, length):

		if length >= self.L_tt:
			print("Tension greater than thick + thin filaments = ", 0, "N/cm2")
			return 0

		elif length > self.L_bt:
			self.Desc_Tension = (self.Desc_slope*length) + self.Desc_y_int
			print("Descending Limb %T0= ", self.Desc_Tension, "N/cm2")
			return self.Desc_Tension

		elif length >= self.L_tz  and length <= self.L_bt:
			print("Plateau %T0 = ", self.T0, "N/cm2")
			return self.T0

		elif length <= self.L_tz:
			self.Asc_Tension = (self.Asc_slope * length) + self.Asc_y_int

			if self.Asc_Tension <= 0.0:
				self.Asc_Tension = 0.0

			print("Ascending Limb %T0 = ", self.Asc_Tension, "N/cm2")
			return self.Asc_Tension

		else:
			return 0

		print('\n', "----------------------------------------", '\n')

	### Calculating the Characteristic Points and returning 
	### Numpy Array of it. Units are as follows --> [ um, N ]
	def getMark(self):

		# Getting the maximum characteristic point
		max_LT = np.hstack((self.L_tt, self.calcT(self.L_tt)))
		#print(max_LT)

		#Getting the Right side of the Plateau of characteristic point
		Rplat_LT = np.hstack((self.L_bt, self.calcT(self.L_bt)))
		#print(Rplat_LT)

		#Getting the Left side of the Plateau of characteristic point
		Lplat_LT = np.hstack((self.L_tz, self.calcT(self.L_tz)))
		#print(Lplat_LT)

		#Getting the minimum characteristic point
		x = Lplat_LT[0] - (Lplat_LT[1]/self.Asc_slope)
		min_LT = np.hstack((x, 0.0))
		#print(min_LT)

		self.marker = np.vstack((max_LT, Rplat_LT, Lplat_LT, min_LT))
		print('\n', "Characteristic Points: ", '\n', self.marker)

		print('\n', "----------------------------------------", '\n')

		return self.marker

	### Calculating Force from a given radius and tension values
	def calcForce(self, tension):
		print("The radius of the muscle = ", self.radius, "cm")
		print("Length of the muscle = ", self.L_muscle, "cm")

		self.tension = tension

		# Cross Sectional Area (CSA) calculations 
		self.CSA = (self.radius**2) * np.pi
		print("Cross Sectional Area = ", self.CSA, "cm2")

		#Calculate Force (Newtons)
		self.Force = self.CSA * self.tension
		print("Force of the Muscle = ", self.Force, "N")

		print('\n', "----------------------------------------", '\n')

		return self.Force

	### Calculating the number of Sarcomeres in series ###
	def numSarc(self):
		# Average of markers in plateau
		A = self.marker[1,0]
		B = self.marker[2,0]

		### Middle of the Plateau for L0 ### 
		L0 = (( A + B ) / 2) / 10000

		self.num_sarcomere = self.L_muscle / L0
		self.halfSarc = self.num_sarcomere * 2

		print("Number of Whole Sarcomeres = ", self.num_sarcomere)
		print("Number of Half Sarcomeres = ", self.halfSarc)

		print('\n', "----------------------------------------", '\n')

		return self.halfSarc, self.num_sarcomere

	### Calculating the Work of the Muscle ###
	def calcWork(self):

		#Descending Limb Region ( in cm )
		max_LD = self.marker[0,0] / 10000
		min_LD = self.marker[1,0] / 10000
		DL = max_LD - min_LD

		#Descending Limb T0
		Desc_T0 = self.calcForce(self.marker[1,1])

		#Plateau Region ( in cm )
		max_LP = self.marker[1,0] / 10000
		min_LP = self.marker[2,0] / 10000
		PL = max_LP - min_LP 

		#Plateau T0
		Plat_T0 = self.calcForce(self.marker[2,1])

		#Ascending Limb Region ( in cm)
		max_LA = self.marker[2,0] / 10000
		min_LA = self.marker[3,0] / 10000
		AL = max_LA - min_LA

		#Ascending Limb T0
		Asc_T0 = self.calcForce(self.marker[2,1])

		print("Calculating the Work of the Muscle: ", '\n')

		### Calculating Work for Limbs Individually ###
		self.Desc_work = (0.5 * ((Desc_T0*DL) / 100)) * self.num_sarcomere
		print("Descending Limb Work = ", self.Desc_work, "J")
		self.Plat_work = ((Plat_T0 * PL) / 100 ) * self.num_sarcomere
		print("Plateau Work = ", self.Plat_work, "J")
		self.Asc_work = (0.5 * ((Asc_T0*DL) / 100)) * self.num_sarcomere
		print("Ascending Limb Work = ", self.Asc_work, "J")

		self.Work = self.Desc_work + self.Plat_work + self.Asc_work
		print('\n', "Calculating the Work for the Sliding Filament Theory = ", 
			self.Work, "J")

		print('\n', "----------------------------------------", '\n')

		return self.Desc_work, self.Plat_work, self.Asc_work, self.Work

	### Calculating the Energy Density of the Muscle ###
	def calcED(self):

		# Calculate the Volume
		muscle_volume = (self.radius**2)* self.L_muscle * (np.pi)
		#print("Volume of the Muscle = ", muscle_volume, "cm3")

		# Calculate the Energy Density of the Muscle
		self.ED = self.Work / muscle_volume
		print("Energy Density of the Muscle = ", self.ED, "J/cm3")

		print('\n', "----------------------------------------", '\n')

		return self.ED

	### Plotting the Characteristic Points ### 
	def plotMark(self):
		x = self.marker[:,0]
		y = self.marker[:,1]

		plt.xlabel('Length (um)', fontsize = 16)
		plt.ylabel('Tension (N/cm2)', fontsize = 16)
		plt.title('Isometric Tetani and Length Curve', fontsize = 20)

		plt.tick_params(axis = 'both', which = 'major', labelsize = 14)

		plt.plot(x,y, 'mp', markersize = 10)
		plt.plot(x,y, 'm--')
		plt.show()

	def main(self):
		print('\n', "----------------------------------------", '\n')
		print("CALCULATING...")
		self.calcLength()
		self.calcLinear()
		self.getMark()
		self.numSarc()
		self.calcWork()
		self.calcED()

### class for the loop Simulation ###
class Sarcomere_Loop:

	def __init__(self, thick=None, thin=None, bare=None, zdisc=None, slope=None, radius=None, LM=None):

		print('\n', "INPUT VALUES ACCORDINGLY", '\n')

		if thick == None:
			thick = float(input("Thick Filament/Myosin Length (um) = "))
		else:
			print("Thick Filament/Myosin Length (um) = ", thick)

		if thin == None:
			thin = float(input("Thin Filament/Actin Length (um) = "))
		else:
			print("Thin Filament/Actin Length (um) = ", thin)

		if bare == None:
			bare = float(input("Bare zone Length (um) = "))
		else:
			print("Bare zone (um) = ", bare)

		if zdisc == None:
			zdisc = float(input("Z Disc length (um) = "))
		else:
			print("Z Disc (um) = ", zdisc)

		if slope == None:
			slope = float(input("Slope (N/(cm2*um) = "))
		else:
			print("Slope (N/(cm2*um) = ", slope)

		if radius == None:
			radius = float(input("Radius of Muscle Model (cm) = "))
		else:
			print("Radius of Muscle Model (cm)= ", radius)

		if LM == None:
			LM = float(input("Length of Muscle Model (cm)= "))
		else:
			print("Length of Muscle Model (cm) = ", LM)

		### Parameters for Sarcomere Lengths and Median Max Tension (T0) ###
		self.L_thick = thick 	# Thick Filament/Myosin Length
		self.L_thin = thin 		# Thin Filament/Actin Length
		self.L_bare = bare 		# Bare Zone Length
		self.L_zdisc = zdisc 	# Z Disc Length
		self.slope = slope 		# Slope
		self.radius = radius 	# Radius of the Muscle Model
		self.L_muscle = LM 		# Lenght of the Muscle Model

		#sum of filament lengths
		self.L_tt = None 		# Thick and Thin Filaments and Z disc
		self.L_bt = None 		# Thin Filament, Bare zone and Z disc
		self.L_tz = None 		# Thin Filament and Z disc

		#slope and y-intercept values
		self.Desc_slope = None 	# Slope of the Descending Limb
		self.Desc_y_int = None 	# Y-intercept of the Descending Limb
		self.Asc_slope = None 	# Slope of the Ascending Limb
		self.Asc_y_int = None 	# Y-intercept of the Ascending Limb

		# Theoretical Tension values
		self.Desc_Tension = None 	# Theoretical Tension of Descending Limb
		self.T0 = None				# Theoretical Max Tension
		self.Asc_Tension = None 	# Theoretical Tension of Ascending Limb

		#Characteristic Points
		self.marker = None 

		#Force calculation
		self.Force = None
		self.F0 = None				# Peak Force Value

		#Number of sarcomeres
		self.num_sarcomere = None 	# Number of Whole Sarcomeres in series
		self.halfSarc = None 		# Number of half sarcomeres in series

		#Work calculation
		self.Desc_work = None 	# Calculated work of Descending Limb
		self.Plat_work = None 	# Calculated work of Plateau region
		self.Asc_work = None 	# Calculated work of Ascending Limb
		self.Work = None 		# Calculated Work Total

		#Energy Density
		self.ED = None 			# Calculated Energy Density

	### Calculate the Length of the Different Filament Lengths ###
	def calcLength(self):
		# Printing the individual Lengths #
		print('\n', "These are the Filament Lengths: ")
		print("Myosin (um) = ", self.L_thick)
		print("Actin (um) = ", self.L_thin)
		print("Bare Zone (um) = ", self.L_bare)
		print("Z disc (um) = ", self.L_zdisc)

		print('\n', "Slope = ", self.slope, " N/(cm2*um) ")

		print('\n', "----------------------------------------", '\n')

		print("The sum of the filament lengths: ")
		
		# Length of thick, thin and z disc combined #
		self.L_tt = self.L_thick + self.L_thin + self.L_zdisc
		print("Sum of Thick, Thin and Zdisc = ", self.L_tt, "um")

		#Length of Bare, Thin and Z disc combined #
		self.L_bt = self.L_thin + self.L_bare + self.L_zdisc
		print("Sum of Thin Filament, Bare Zone and Zdisc = ", self.L_bt,
			"um")

		#Length of Thin Filament and Z disc combined #
		self.L_tz = self.L_thin + self.L_zdisc
		print("Sum of Thin Filament and Zdisc = ", self.L_tz, "um")

		#If statement for L3, L2 and L1 
		self.L3 = self.L_tt

		if self.L_thin < (self.L_thick - self.L_bare):
			self.L2 = self.L_thick + self.L_zdisc
			self.L1 = self.L2

		else:
			self.L2 = self.L_bt	# The right side of the plateau
			self.L1 = self.L_tz

		return self.L_tt, self.L_bt, self.L_tz, self.L1, self.L2, self.L3

	### Calculate the Slope and Y-intercept of the Descending and Ascending Limb ###
	def calcLinear(self):
		print('\n', "----------------------------------------", '\n')

		# Slope of the Descending Limb
		self.Desc_slope = self.slope
		print("Slope of the Descending Limb = ", self.Desc_slope)

		# Y-intercept of the Descending Limb
		self.Desc_y_int = (self.Desc_slope*(-1))*self.L_tt
		print("Y-intercept of the Descending Limb = ", self.Desc_y_int)

		# Slope of the Ascending Limb
		self.Asc_slope = (-1)*self.slope
		print('\n', "Slope of the Ascending Limb = ", self.Asc_slope)

		#Max Tension Value
		self.T0 = (self.Desc_slope*(-1)) * (self.L_tt-self.L_bt)

		#Y-intercept of the Ascending Limb
		self.Asc_y_int = self.T0 - (self.Asc_slope)*(self.L_tz)
		print("Y-intercept of the Ascending Limb = ", self.Asc_y_int)

		print('\n', "----------------------------------------", '\n')

		return self.Desc_slope, self.Asc_slope, self.Desc_y_int, self.Asc_y_int

	### Calculating the Theoretical Tension in Newtons ###
	def calcT(self, length):

		if length >= self.L_tt:
			print("Tension greater than thick + thin filaments = ", 0, "N/cm2")
			return 0

		elif length > self.L_bt:
			self.Desc_Tension = (self.Desc_slope*length) + self.Desc_y_int
			print("Descending Limb = ", self.Desc_Tension, "N/cm2")
			return self.Desc_Tension

		elif length >= self.L_tz  and length <= self.L_bt:
			self.T0 = (self.Desc_slope*-1)*(self.L_tt - self.L_bt)
			print("Plateau = ", self.T0, "N/cm2")
			return self.T0

		elif length <= self.L_tz:
			self.Asc_Tension = (self.Asc_slope * length) + self.Asc_y_int

			if self.Asc_Tension <= 0.0:
				self.Asc_Tension = 0.0

			print("Ascending Limb = ", self.Asc_Tension, "N/cm2")
			return self.Asc_Tension

		else:
			return 0

		print('\n', "----------------------------------------", '\n')

	### Calculating the Characteristic Points and returning 
	### Numpy Array of it. Units are as follows --> [ um, N ]
	def getMark(self):

		x = self.L3 - self.L2
		print("The difference of two lengths = ", x, "um")

		# calculate the Force (Newtons)
		self.T0 = np.absolute(self.slope * x)

		# Getting the maximum characteristic point
		max_LT = np.hstack((self.L3, self.calcT(self.L3)))
		#print(max_LT)

		#Getting the Right side of the Plateau of characteristic point
		Rplat_LT = np.hstack((self.L2, self.T0))#self.calcT(self.L2)))
		#print(Rplat_LT)

		#Getting the Left side of the Plateau of characteristic point
		Lplat_LT = np.hstack((self.L1, self.T0))#self.calcT(self.L1)))
		#print(Lplat_LT)

		#Getting the minimum characteristic point
		x = Lplat_LT[0] - (Lplat_LT[1]/self.Asc_slope)
		min_LT = np.hstack((x, 0.0))
		#print(min_LT)

		self.marker = np.vstack((max_LT, Rplat_LT, Lplat_LT, min_LT))
		print('\n', "Characteristic Points: ", '\n', self.marker)

		print('\n', "----------------------------------------", '\n')

		return self.marker

	### Calculating Force from a given radius and tension values
	def calcForce(self, tension):
		print("The radius of the muscle = ", self.radius, "cm")
		print("Length of the muscle = ", self.L_muscle, "cm")

		self.tension = tension

		# Cross Sectional Area (CSA) calculations 
		self.CSA = (self.radius**2) * np.pi
		print("Cross Sectional Area = ", self.CSA, "cm2")

		#Calculate Force (Newtons)
		self.Force = self.CSA * self.tension
		print("Force of the Muscle = ", self.Force, "N")

		print('\n', "----------------------------------------", '\n')

		return self.Force

	### Calculating the peak Force ###
	def calcF0(self):
		print("The radius of the muscle = ", self.radius, "cm")
		print("Length of the muscle = ", self.L_muscle, "cm")

		#Cross Sectional Area (CSA) calculations
		self.CSA = (self.radius**2) * np.pi
		print("Cross Sectioal Area = ", self.CSA, "cm2")

		# Getting the difference of the two filament lengths
		# using the following lengths
		
		print("The max length of the sarcomere = ", self.L3, "um")
		print("The length at the right side of the plateau = ", self.L2, "um")

		x = self.L3 - self.L2
		print("The difference of two lengths = ", x, "um")

		# calculate the Force (Newtons)
		self.F0 = np.absolute(self.slope * self.CSA * x)
		print("Force of the Muscle = ", self.Force, "N")

		print('\n', "----------------------------------------", '\n')

		return self.F0

	### Calculating the number of Sarcomeres in series ###
	def numSarc(self):

		### Middle of the Plateau for L0 ### 
		L0 = (( self.L1 + self.L2 ) / 2) / 10000

		self.num_sarcomere = self.L_muscle / L0
		self.halfSarc = self.num_sarcomere * 2

		print("Number of Whole Sarcomeres = ", self.num_sarcomere)
		print("Number of Half Sarcomeres = ", self.halfSarc)

		print('\n', "----------------------------------------", '\n')

		return self.halfSarc, self.num_sarcomere

	### Calculating the Work of the Muscle ###
	def calcWork(self):

		#Descending Limb Region ( in cm )
		max_LD = self.L3 / 10000
		min_LD = self.L2/ 10000
		DL = max_LD - min_LD

		#Descending Limb T0
		Desc_T0 = self.calcForce(self.marker[1,1])

		#Plateau Region ( in cm )
		max_LP = self.L2 / 10000
		min_LP = self.L1 / 10000
		PL = max_LP - min_LP 

		#Plateau T0
		Plat_T0 = self.calcForce(self.marker[2,1])

		#Ascending Limb Region ( in cm)
		max_LA = self.L1 / 10000
		min_LA = self.marker[3,0] / 10000
		AL = max_LA - min_LA

		#Ascending Limb T0
		Asc_T0 = self.calcForce(self.marker[2,1])

		print("Calculating the Work of the Muscle: ", '\n')

		### Calculating Work for Limbs Individually ###
		self.Desc_work = (0.5 * ((Desc_T0*DL) / 100)) * self.num_sarcomere
		print("Descending Limb Work = ", self.Desc_work, "J")
		self.Plat_work = ((Plat_T0 * PL) / 100 ) * self.num_sarcomere
		print("Plateau Work = ", self.Plat_work, "J")
		self.Asc_work = (0.5 * ((Asc_T0*DL) / 100)) * self.num_sarcomere
		print("Ascending Limb Work = ", self.Asc_work, "J")

		self.Work = self.Desc_work + self.Plat_work + self.Asc_work
		print('\n', "Calculating the Work for the Sliding Filament Theory = ", 
			self.Work, "J")

		print('\n', "----------------------------------------", '\n')

		return self.Desc_work, self.Plat_work, self.Asc_work, self.Work

	### Calculating the Energy Density of the Muscle ###
	def calcED(self):

		# Calculate the Volume
		muscle_volume = (self.radius**2)* self.L_muscle * (np.pi)
		#print("Volume of the Muscle = ", muscle_volume, "cm3")

		# Calculate the Energy Density of the Muscle
		self.ED = self.Work / muscle_volume
		print("Energy Density of the Muscle = ", self.ED, "J/cm3")

		print('\n', "----------------------------------------", '\n')

		return self.ED

	### Plotting the Characteristic Points ### 
	def plotMark(self):
		x = self.marker[:,0]
		y = self.marker[:,1]

		plt.xlabel('Length (um)', fontsize = 16)
		plt.ylabel('Tension (N/cm2)', fontsize = 16)
		plt.title('Isometric Tetani and Length Curve', fontsize = 20)

		plt.tick_params(axis = 'both', which = 'major', labelsize = 14)

		plt.plot(x,y, 'mp', markersize = 10)
		plt.plot(x,y, 'm--')
		plt.show()

	def MAIN(self):
		print('\n', "----------------------------------------", '\n')
		print("CALCULATING...")
		self.calcLength()
		self.calcLinear()
		self.getMark()
		self.calcF0()
		self.numSarc()
		self.calcWork()
		self.calcED()

	### The code to use for the nested for loop! ###
	###			input: thick filament, thin filament lengths ###
	def main_loop(self, Lthick=None, Lthin=None):
		print('\n', "----------------------------------------", '\n')

		if Lthick == None:													# If not inputted when called,
			Lthick = float(input("Thick Filament/Myosin Length (um) = ")) 	# 
		else: 
			print("Thick Filament/Myosin Length (um) = ", Lthick)

		if Lthin == None:
			Lthin = float(input("Thin Filament/Actin (um) = "))
		else:
			print("Thin Filament/Actin (um) = ", Lthin)

		self.calcLength()
		self.getMark()
		self.calcF0()		# Make sure to copy and paste the slope from the main() function!
		self.numSarc()
		#self.calcWork()
		#self.calcED()

if __name__ == "__main__":
	### Input the filament lengths in following order...
			### thick, thin, bare, zdisc, T0, radius, muscle length ###
	sarc_ls = []		# list of [thick, thin] filament lengths in um
	num_sarc = []		# list of shortening velocity coefficients through number of sarcomeres
	sarc_work = []		# list of work values (J)
	sarc_ed = []		# list of energy density values (J/cm3)
	sarc_force = []		# list of force values (N)
	
	Lthick = [1.0, 2.0, 3.0]
	Lthin = [1.0, 2.0, 3.0]

	Thick = []			# 2D array of Thick filament lengths

	Thin = [] 			# 2D array of Thin filament lengths

	row=0
	col=0

	### Inputting the Organism and Muscle of Interest ###
	organism = str(input("Organism and Muscle of Interest: "))

	sarcomere = Sarcomere_Simulation()
	sarcomere.main()

	plt.xlabel('Length (um)', fontsize = 16)
	plt.ylabel('Force (N)', fontsize = 16)
	plt.title('Isometric Tetani and Length Curve', fontsize = 20)

	plt.tick_params(axis = 'both', which = 'major', labelsize = 14)

	### Nested For Loop for the different muscle lengths ###
	for m in Lthick:
		row += 1
		for a in Lthin:
			col += 1 
			sarc = Sarcomere_Loop(m, a, 0.2, 0.05, sarcomere.Desc_slope, 1.0, 10.0)
			sarc.MAIN()
			sarc_ls.append([m, a])
			x = sarc.marker[:,0]
			y = sarc.marker[:,1]
			plt.plot(x,y, 'p', markersize = 10)
			plt.plot(x,y, '--')
			sarc_force.append([sarc.F0])
			sarc_work.append([sarc.Work])
			num_sarc.append([sarc.num_sarcomere])
			sarc_ed.append([sarc.ED])
			Thick.append([m])
			Thin.append([a])
	F0_cont = np.column_stack([Thick, Thin, sarc_force])
	Work_cont = np.column_stack([Thick, Thin, sarc_work])
	ED_cont = np.column_stack([sarc_ed, sarc_ed, sarc_ed])
	X = np.column_stack([Thick, Thick, Thick])
	Y = np.column_stack([Thin, Thin, Thin])

	plt.show()

	print('\n', "List of Sarcomere Lengths: ", *sarc_ls, sep = "\n")
	print('\n', "List of Force done: ", *sarc_force, sep = "\n")
	print('\n', "Shortening Velocity Coefficients: ", *num_sarc, sep="\n")
	print('\n', "List of Work done through Loop: ",*sarc_work, sep = "\n")
	print('\n', "List of Energy Densities: ", *sarc_ed, sep="\n")


	## The section to help plot for the Force Contour plot ##

	fig = plt.figure()
	ax = fig.add_subplot(111, projection = '3d')
	Force_cont = ax.contour(X, Y, F0_cont) #Contour plot of the force
	#plt.clabel(Force_cont, fontsize = 12)
	plt.title('Force Contour Plot', fontsize = 20)
	plt.xlabel('Thick Filament (um)' , fontsize = 16)
	plt.ylabel('Thin Filament (um)', fontsize = 16)
	plt.show()

	## The section to help plot for the Work Contour plot ##

	fig1 = plt.figure()
	ax1 = fig1.add_subplot(111, projection= '3d')
	W_cont = ax1.contour(X, Y, Work_cont)
	#plt.clabel(Work_cont, fontsize = 12)
	plt.title("Work Contour Plot", fontsize = 20)
	plt.xlabel('Thick Filament (um)', fontsize = 14)
	plt.ylabel('Thin Filament (um)', fontsize = 14)
	plt.show()


			

			




