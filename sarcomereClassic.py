# sarcomereClassic.py : A sliding filament model of a single sarcomere's isometric
# length-tension relationship. Visualizes the length-tension curve with 
# matplotlib.
#
# To generate response surfaces, run: 	python3 sarcomereClassic.py response
# To simulate a single sarcomere model: python3 sarcomereClassic.py 
#
# Caitrin Eaton
# 13 July 2018

import random # only used in responseSurface() for probability of plotting an isometric curve

import sys
import numpy as np
import matplotlib.pyplot as plt

class SlidingFilament:
	
	# Default sarcomere model characteristics drawn from Heejoon's 
	# median values for vertebrates.
	lthick = 1.6			# thick filament length, um
	lthin = 2.2			# thin filament length, um
	lz = 0.08			# z disc thickness, um
	lbare = 0.155			# bare zone thickness, um
	t0 = None			# peak isometric tension, N/cm**2
	slope = -17.1922		# slope of the isometric length-tension curve's segments, Ncm**-2/um_overlap
	csv = None			# empirical data file
		
	def __init__( self, lthick=None, lthin=None, lz=None, lbare=None, slope=None, t0=None, csv=None):
		''' CAUTION: Can pass in either peak tension t0 or slope of the descending limb.
					 If both are present (not None), slope will be given priority, and 
					 a theoretical self.t0 will be calculated based on maximal theoretical
					 crossbridge-populated overlap and the given slope.'''
		# Constants
		if lthick != None:
			self.lthick = lthick
		if lthin != None:
			self.lthin = lthin
		if lz != None:
			self.lz = lz
		if lbare != None:
			self.lbare = lbare
		if slope != None:
			self.slope = np.array([-slope, 0, slope])
			self.t0 = None
		elif t0 != None:
			self.t0 = t0
			self.slope = np.zeros((3,))
		else:
			self.t0 = None
			self.slope = np.zeros((3,))
		if csv != None:
			self.csv = csv
			data = np.genfromtxt(self.csv, delimiter=',', autostrip=True, skip_header=1)
			self.L_empirical = data[:,0]
			self.Ta_empirical = data[:,1]
			
		self.calcCharacteristicPoints()	# Vertices of the sliding filament model's graph
		self.calcWork()
		self.calcRMSE()		# R-squared of model relative to empirical data
		self.calcRsquared()	# Root mean squared error of model relative to empirical data
			
	def calcCharacteristicPoints( self ):
		'''Calculate this sarcomere geometry's characteristic points'''
		self.L = np.zeros((4,))								# Sarcomere length, um
		self.Ta = np.zeros((4,))							# Active tension, N/cm**2
		
		# Characteristic lengths
		lcb = min(self.lthin, self.lthick - self.lbare)		# length of maximum theoretical crossbridge-populated overlap
		self.L[3] = self.lthick + self.lthin + self.lz		# x-intercept of the descending limb: too long, no overlap
		self.L[2] = self.L[3] - lcb							# right-hand side of the bare zone: longest length maintaining max overlap
		self.L[1] = max(self.lthick, self.lthin) + self.lz	# left-hand side of the bare zone: shortest length maintaining max overlap
		self.L[0] = self.L[1] - lcb							# x-intercept of the ascending limb: too short, compression negates crossbridge forces. 
															# 		Note that In the classic sliding filament model, we have no mechanism by which 
															# 		to predict ascending limb slope, so we mirror the descendling limb.
		
		# Limb slopes and peak tension
		if (np.min(self.slope) == 0) and (self.t0 != None):	
			# Isometric length-tension curve slopes, um*N/cm**2
			self.slope[2] = self.t0/(self.L[2]-self.L[3])	# descending limb
			self.slope[1] = 0.0								# plateau
			self.slope[0] = -self.slope[2]					# ascending limb
		elif (self.t0 == None) and (np.min(self.slope) != 0):
			# Peak active isometric tension (N/cm**2)
			self.t0 = -self.slope[2] * lcb
		else:
			# No slope or tension to go off of, so normalize tension
			self.t0 = 1.0
			self.slope[2] = self.t0/(self.L[2]-self.L[3])	# descending limb
			self.slope[1] = 0.0								# plateau
			self.slope[0] = -self.slope[2]					# ascending limb
		
		# Characteristic tensions
		self.Ta[3] = 0.0
		self.Ta[2] = self.t0
		self.Ta[1] = self.t0
		self.Ta[0] = 0.0
		
		# "Rest length" can be chosen as the right-hand side of the plateau
		# or the middle of the plateau, as long as we're consistent throughout.
		# This length gets used as the mystical muscle cylinder's sarcomere
		# packing length in the full muscle model.
		self.l0 = self.L[2]								# right-hand side of the plateau
		#self.l0 = np.mean(self.L[1:3])					# middle of the plateau
		
		# If we've recalculated characteristic points, then we must
		# recalculate work.
		self.calcWork()
	
	def calcWork( self, lmax=None, lmin=None ):
		'''Calculate the max theoretical specific work (in J/cm**2) done by the 
		sarcomere as it shortens across its entire force-producing 
		range.'''
		self.wDesc = 0.5*self.t0*(self.L[3]-self.L[2])*10**-6
		self.wPlat = self.t0*(self.L[2]-self.L[1])*10**-6
		self.wAsc  = 0.5*self.t0*(self.L[1]-self.L[0])*10**-6
		self.w = self.wDesc + self.wPlat + self.wAsc
		return self.w
	
	def calcTa( self, l ):
		'''Calculate the active isometric tension (in N/cm**2) at length
		l (in um).'''
		if type(l)==float or type(l)==int:
			# single length value
			if l > self.L[3]:
				# Too long, no overlap
				return 0
			elif l > self.L[2]:
				# Descending limb
				return self.t0 + (l - self.L[2])*self.slope[2]
			elif l >= self.L[1]:
				# Plateau
				return self.t0
			elif l > self.L[0]:
				# Ascending limb
				return (l-self.L[0])*self.slope[0]
			# Too short, no active force production
			return 0
		else:
			# ndarray of length values
			t = np.zeros(l.shape)
			# Descending limb
			desc = (self.L[2] < l) & (l < self.L[3] )
			t[ desc ] = self.t0 + (l[desc] - self.L[2])*self.slope[2]
			# Plateau
			plat = (self.L[1] <= l) & (l <= self.L[2] )
			t[ plat ] = self.t0
			# Ascending limb
			asc = (self.L[0] < l) & (l < self.L[1])
			t[ asc ] = (l[ asc ]-self.L[0])*self.slope[0]
			return t
	
	def calcRMSE( self ):
		'''Calculate the root mean squared error of the model
		relative to the empirical data ("truth").'''
		if self.csv == None:
			# No empirical data available
			self.rmse = None
			return self.rmse
		
		tmodel = self.calcTa( self.L_empirical )
		self.rmse = np.sqrt(np.mean(np.square(self.Ta_empirical - tmodel)))
		
		# Calculate RMSE for each limb:
		
		# Descending limb
		desc = (self.L[2] < self.L_empirical) & (self.L_empirical < self.L[3] )
		ttrue = self.Ta_empirical[ desc ]
		tmodel = self.calcTa( self.L_empirical[ desc ] )
		self.rmseDesc = np.sqrt(np.mean(np.square(ttrue - tmodel)))
		
		# Plateau
		plat = (self.L[1] <= self.L_empirical) & (self.L_empirical <= self.L[2] )
		ttrue = self.Ta_empirical[ plat ]
		tmodel = self.calcTa( self.L_empirical[ plat ] )
		self.rmsePlat = np.sqrt(np.mean(np.square(ttrue - tmodel)))
		
		# Ascending limb
		asc = (self.L[0] < self.L_empirical) & (self.L_empirical < self.L[1])
		ttrue = self.Ta_empirical[ asc ]
		tmodel = self.calcTa( self.L_empirical[ asc ] )
		self.rmseAsc = np.sqrt(np.mean(np.square(ttrue - tmodel)))
		
		return (self.rmse, self.rmseAsc, self.rmsePlat, self.rmseDesc)
	
	def calcRsquared( self, plot=True ):
		'''Calculate the R**2 correlation coefficient of the model
		relative to the empirical data ("truth").'''
		if self.csv == None:
			# No empirical data available
			self.rsq = None
			return self.rsq
		
		# Overall R-squared of the model relative to the entire dataset:
		tmean = np.mean(self.Ta_empirical)
		tmodel = self.calcTa( self.L_empirical )
		sstot = np.sum(np.square(self.Ta_empirical - tmean))
		ssres = np.sum(np.square(self.Ta_empirical - tmodel))
		self.rsq = 1 - ssres/sstot
		if plot:
			# Superimpose the residuals on a length-tension curve.
			# Also helps us verify that the tensions we're calculating
			# lie along the model's line segments
			fig = plt.figure()	
			plt.figure(fig.number)
			plt.xlabel("Length (um)")
			plt.ylabel("Tension (N/cm**2)")
			plt.title("Residuals")	
			plt.plot(self.L, self.Ta, 'd-k', markerfacecolor='k')	
			plt.plot(self.L_empirical, self.Ta_empirical, 'ob', markerfacecolor='w')	
			for i in range(len(self.L_empirical)):
				plt.plot([self.L_empirical[i], self.L_empirical[i]], [tmodel[i], self.Ta_empirical[i]], ':r')	
				plt.plot(self.L_empirical[i], tmodel[i], '.r')
				plt.plot(self.L_empirical[i], tmodel[i], '.r')
		
		# Calculate R-squared for each limb:
		
		# Descending limb
		desc = (self.L[2] < self.L_empirical) & (self.L_empirical < self.L[3] )
		l = self.L_empirical[ desc ]
		ttrue = self.Ta_empirical[ desc ]
		#tmean = np.mean( ttrue )
		tmodel = self.calcTa( l )
		sstot = np.sum(np.square(ttrue - tmean))
		ssres = np.sum(np.square(ttrue - tmodel))
		self.rsqDesc = 1 - ssres/sstot
		
		# Plateau
		plat = (self.L[1] <= self.L_empirical) & (self.L_empirical <= self.L[2] )
		l = self.L_empirical[ plat ]
		ttrue = self.Ta_empirical[ plat ]
		#tmean = np.mean( ttrue )
		tmodel = self.calcTa( l )
		sstot = np.sum(np.square(ttrue - tmean))
		ssres = np.sum(np.square(ttrue - tmodel))
		self.rsqPlat = 1 - ssres/sstot
		
		# Ascending limb
		asc = (self.L[0] < self.L_empirical) & (self.L_empirical < self.L[1])
		l = self.L_empirical[ asc ]
		ttrue = self.Ta_empirical[ asc ]
		#tmean = np.mean( ttrue )
		tmodel = self.calcTa( l )
		sstot = np.sum(np.square(ttrue - tmean))
		ssres = np.sum(np.square(ttrue - tmodel))
		self.rsqAsc= 1 - ssres/sstot
		
		return (self.rsq, self.rsqAsc, self.rsqPlat, self.rsqDesc)
					
	# Accessors
		
	def get_lthick( self ):
		'''Return the thick filament length, in um'''
		return self.lthick
		
	def get_lthin( self ):
		'''Return the thin filament length, in um'''
		return self.lthin
		
	def get_lz( self ):
		'''Return the z disc width, in um'''
		return self.lz
		
	def get_lbare( self ):
		'''Return the bare zone width, in um'''
		return self.lbare
		
	def get_t0( self ):
		'''Return the peak isometric tension, in N/cm**2'''
		return self.t0
		
	def get_slope( self ):
		'''Return the slopes of the ascending limb, plateau, and 
		descending limb'''
		return self.slopes.copy()
		
	def get_csv( self ):
		'''Return the name of the empirical data file'''
		return self.csv
		
	# Utilities
	
	def print( self ):
		'''Describe model characteristics with string formatting.'''
		stateStr  = "\n------------------------------------------------"
		stateStr += "\n\t\tsarcomere geometry"
		stateStr += "\n------------------------------------------------"
		stateStr += "\nL thick\t{0:<5.3f} um".format(self.lthick)
		stateStr += "\nL thin\t{0:<5.3f} um".format(self.lthin)
		stateStr += "\nL lz\t{0:<5.3f} um".format(self.lz)
		stateStr += "\nL bare\t{0:<5.3f} um".format(self.lbare)
		stateStr += "\nL0\t{0:<5.3f} um".format(self.l0)
		stateStr += "\nT0\t{0:<5.3f} N/cm**2".format(self.t0)
		stateStr += "\n\ncharacteristic length-tension pairs"
		stateStr += "\n\t{0} um".format(self.L.T)
		stateStr += "\n\t{0} N/cm**2".format(self.Ta.T)
		stateStr += "\n\nWork\t{0:<10.9f} J/cm**2".format(self.w)
		stateStr += "\n\n\tAscending  {0:<10.9f} J/cm**2".format(self.wAsc)
		stateStr += "\n\tPlateau    {0:<10.9f} J/cm**2".format(self.wPlat)
		stateStr += "\n\tDescending {0:<10.9f} J/cm**2".format(self.wDesc)
		if self.rmse != None:
			stateStr += "\n\nIsometric curve RMSE {0:<10.9f}".format(self.rmse)
			stateStr += "\n\n\tAscending  {0:<10.9f}".format(self.rmseAsc)
			stateStr += "\n\tPlateau    {0:<10.9f}".format(self.rmsePlat)
			stateStr += "\n\tDescending {0:<10.9f}".format(self.rmseDesc)
		if self.rsq != None:
			stateStr += "\n\nIsometric curve R**2 {0:<10.9f}".format(self.rsq)
			stateStr += "\n\n\tAscending  {0:<10.9f}".format(self.rsqAsc)
			stateStr += "\n\tPlateau    {0:<10.9f}".format(self.rsqPlat)
			stateStr += "\n\tDescending {0:<10.9f}".format(self.rsqDesc)
		stateStr += "\n"
		print( stateStr )
		
	def isStable( self ):
		'''Decide whether the model is stable (True) or whether it has 
		gone irrevocably awry (False). Failure is defined as a tension
		in excess of 1.8*T0 (McMahon 1974).'''
		if self.t > 1.8*self.t0:
			# McMahon thinks the muscle will yield at 1.8*T0
			print("FAILURE: tension in excess of 1.8*T0")
			return False
		# Lookin' good!
		return True
		
	def plot( self, fig=None ):
		'''Display the sliding filament model's curve. If a figure
		window is passed in through fig, plot in that window. Otherwise,
		create a new window. Return the figure handle, either way.'''
		if fig == None:
			fig = plt.figure()
			
		plt.figure(fig.number)
		plt.xlabel("Length (um)")
		plt.ylabel("Tension (N/cm**2)")
		plt.title("Classic sliding filament model")	
		plt.plot(self.L, self.Ta, 'dk', markerfacecolor='r')	
		
		if self.csv != None:
			# If empirical data is normalized, un-normalize it
			# (Keep in mind that sometimes re-digitized normalized measurements stray above 1.0)
			plt.plot(self.L_empirical, self.Ta_empirical, "ob", markerfacecolor='w')
				
		plt.plot(self.L, self.Ta, '-k')
		
		return fig
	
	
def responseSurfaces():
	'''Generating response surfaces by simulating sliding filament 
	models for a range of thick and thin filament combinations.'''
	
	# Gordon's 1966 frog semitendinosus descending limb slope
	# slopeDescending = -18.2	# N/cm**2 per micrometer overlap
	
	# Median vertebrate descending limb slope (Heejoon's reported value)
	slopeDescending = -17.1922
	
	# Thick and thin filament lengths to simulate, in um
	thick = np.linspace(0.1, 20.1, num=201)
	thin = np.linspace(0.1, 20.1, num=201)
	
	# 2D arrays to aggregate performance metrics
	nThick = max(thick.shape)
	nThin = max(thin.shape)
	tension = np.zeros( (nThin, nThick) )
	work = np.zeros( (nThin, nThick) )
	restLength = np.zeros( (nThin, nThick) )
	X, Y = np.meshgrid(thick, thin)
	
	# Figure in which to aggregate sliding filament models
	fig = plt.figure()
	plt.xlabel("Length (um)")
	plt.ylabel("Tension (N/cm**2)")
	plt.title("Classic sliding filament model")	
	
	
	for m in range(nThick):
		for a in range(nThin):
			
			# Classic sliding filament model for this new combination of filament lengths
			s = SlidingFilament( lthick=thick[m], lthin=thin[a], slope=slopeDescending )
			
			# Add it to the plot of aggregated sliding filament models
			if random.random() < 0.0005:
				plt.plot(s.L, s.Ta, 'd-', markerfacecolor='w')
				plt.pause(0.01)
			
			# Peak active isometric tension
			tension[a][m] = s.t0
			
			# Specific work (J/cm**2)
			work[a][m] = s.w
			
			# Rest length correlates inversely with shortening velocity
			restLength[a][m] = s.l0
	
	# Peak tension response surface
	fig = plt.figure()
	plt.xlabel("Thick filament (um)")
	plt.ylabel("Thin filament (um)")	
	plt.title("Classic model's predicted\npeak tension (N/cm**2)")
	contourT = plt.contourf(X, Y, tension)
	plt.colorbar(contourT)
	
	# Specific work response surface
	fig = plt.figure()
	plt.xlabel("Thick filament (um)")
	plt.ylabel("Thin filament (um)")	
	plt.title("Classic model's predicted\nspecific work (J/cm**2)")
	contourW = plt.contourf(X, Y, work)
	plt.colorbar(contourW)
	
	# Rest length response surface
	fig = plt.figure()
	plt.xlabel("Thick filament (um)")
	plt.ylabel("Thin filament (um)")	
	plt.title("Classic model's predicted\nrest length (um)")
	contourL = plt.contourf(X, Y, restLength)
	plt.colorbar(contourL)
	
	plt.show()
	
	

def main( argv ):
	'''Test the simulation by creating a single sliding filament model
	and analyzing its characteristics.'''
	
	# Check to see if we should run the responseSurfaces() function
	if argv[1].strip() == "response":
		responseSurfaces()
		return
	
	# Default sarcomere characteristics
	lthick 	= SlidingFilament.lthick
	lthin	= SlidingFilament.lthin
	lz		= SlidingFilament.lz
	lbare	= SlidingFilament.lbare
	slope	= SlidingFilament.slope
	t0		= SlidingFilament.t0
	csv		= SlidingFilament.csv
	
	# I can never remember what all these are, so print them every time
	print("\nUSAGE: python3 sarcomereClassic.py [optional: -flag <value> pairs]\n")
	print("FLAGS:")
	print("\t-m <thick filament length, um>, default: {} um".format(SlidingFilament.lthick))
	print("\t-a <thin filament length, um>, default: {} um".format(SlidingFilament.lthin))
	print("\t-z <z disc width, um>, default: {} um".format(SlidingFilament.lz))
	print("\t-b <bare zone width, um>, default: {} um".format(SlidingFilament.lbare))
	print("\t-s <descending limb slope, Ncm**-2/um>, default: {}".format(SlidingFilament.slope))
	print("\t-t <peak isometric tension, N/cm**2>, default: {}".format(SlidingFilament.t0))
	print("\t-f <empirical data file, *.csv>, default: {}".format(SlidingFilament.csv))
			
	# Check commandline arguments for flags. Use flags to update model 
	# characteristics before simulation.
	for i in range(1,len(argv)-1):
		flag = argv[i].strip()
		flagParam = argv[i+1].strip()
		if flag == "-m":
			lthick = float(flagParam)
		elif flag == "-a":
			lthin = float(flagParam)
		elif flag == "-z":
			lz = float(flagParam)
		elif flag == "-b":
			lbare = float(flagParam)
		elif flag == "-s":
			slope = float(flagParam)
			t0 = None
		elif flag == "-t":
			t0 = float(flagParam)
			slope = None
		elif flag == "-f":
			csv = flagParam
	
	# Construct new SlidingFilament model
	sf = SlidingFilament(lthick, lthin, lz, lbare, slope, t0, csv)
	
	#sf.calcCharacteristicPoints()
	
	# Print and plot sarcomere geometry info
	sf.print()
	sf.plot()
	plt.show()
	
	
if __name__=="__main__":
	main( sys.argv )
		
