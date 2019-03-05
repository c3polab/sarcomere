# muscleAnchored.py : A muscle model composed of homogeneous sarcomeres.
# The sliding filament template for isometric contractions is "anchored" here with
# ascending limb slopes observed in physiological studies of vertebrate muscle.
# The force of a contracting muscle is modeled at a range of different lengths. 
# Myofilament geometries are constant (homogeneous) across all sarcomeres in the 
# muscle. Visualizes the length-tension curve with matplotlib.
#
#
# To model a single muscle's isometric length-tension curve, run with:  
# >> python3 muscleAnchored.py
# Optional command line argument [-flag <value>] pairs:
#	-m <myosin (thick) filament length, in micrometers>
#	-a <actin (thin) filament length, in micrometers>
#	-z <z-disk width, in micrometers>
#	-b <bare zone width, in micrometers>
#	-r <radius of the relaxed muscle cylinder, in cm>
#	-l <length of the relaxed muscle cylinder, in cm>
#	-f <filename containing animal data that should be superimposed on the model>
#
#
# To generate response surfaces illustrating performance metrics (peak force,
# velocity, energy density) as a function of myofilament lengths, run with:
# >> python3 muscleAnchored.py response
#
#
# Caitrin Eaton
# 13 July 2018


import random # only used in test() for probability of plotting

import sys
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import SFAnchored as anchored


class Homogeneous:
	
	# Default sarcomere geometry
	lthick = 1.6	# thick filament length, um
	lthin = 2.0		# thin filament length, um
	lz = 0.1		# z disc thickness, um
	lbare = 0.15	# bare zone thickness, um
	csv = None		# empirical data file
	rCyl = 1.0		# radius of the muscle cylinder, cm
	lCyl = 10.0		# length of the muscle cylinder, cm
	
	# Median vertebrate isometric curve segment slopes (Heejoon's reported values)
	slopeVertebrate = [67.2000, 6.1935, 0, 0, -17.3913]
		
	def __init__( self, lthick=None, lthin=None, lz=None, lbare=None, slope=None, t0=None, csv=None, rCyl=None, lCyl=None, verbose=False):
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
			self.slope = slope.copy()
			self.t0 = None
		elif t0 != None:
			self.t0 = t0
			self.slope = None
		else:
			self.t0 = None
			self.slope = self.slopeVertebrate.copy()
		if csv != None:
			self.csv = csv
		if rCyl != None:
			self.rCyl = rCyl
		if lCyl != None:
			self.lCyl = lCyl
		self.verbose = verbose
		self.init_sarcomere_anchored()
		self.calcWork()
		self.calcSpeed()
		
	def init_sarcomere_anchored( self ):
		'''Initialize the muscle's sarcomere using the anchored sliding filament template'''
		self.sarcomere = anchored.SlidingFilament( self.lthick, self.lthin, self.lz, self.lbare) #, self.slope, self.csv)
		if self.verbose:
			self.sarcomere.print()
		self.calcWork()
		self.calcSpeed()
					
	def calcWork( self ):
		'''Calculate the theoretical maximum work this muscle could
		produce as a multiple of the sarcomeres' max work.'''
		
		self.csa = np.pi * self.rCyl**2		# cross sectional area of the muscle, cm**2
		l0 = self.sarcomere.l0 * 10**-4		# optimal sarcomere length, cm
		self.Ns = self.lCyl / l0			# number of sarcomeres along the muscle's length
		self.Nhs = self.Ns*2.0				# velocity coefficient: number of half sarcomeres along the muscles length
		
		self.L = self.sarcomere.L * 10**-4 * self.Ns	# length of muscle cylinder at sarcomeres' characteristic lengths
		self.F = self.sarcomere.Ta * self.csa			# force of muscle cylinder at sarcomeres' characteristic lengths
		
		# Work can be calculated as the area under the limbs' 3 line segments
		# Descending limb:
		self.f0 = self.sarcomere.t0*self.csa
		self.wDesc = self.sarcomere.wDesc*self.csa *self.Ns
		self.wPlat = self.sarcomere.wPlat*self.csa *self.Ns
		self.wAsc3  = self.sarcomere.wAsc3*self.csa *self.Ns
		self.wAsc2  = self.sarcomere.wAsc2*self.csa *self.Ns
		self.wAsc1  = self.sarcomere.wAsc1*self.csa *self.Ns
		self.w = self.wDesc + self.wPlat + self.wAsc3 + self.wAsc2 + self.wAsc1;
		
		# Energy density!
		self.energyDensity = self.w / (self.csa * self.lCyl) # J / cm**3
		
		return self.w
					
	def calcNhs( self ):
		'''Calculate the number of half sarcomeres in series along the 
		length of the muscle cylinder.'''
		l0 = self.sarcomere.l0 / 2.0 * 1e-4	# convert um to cm
		self.Nhs = self.lCyl / l0	# number of sarcomeres along the muscle's length
		return self.Nhs
		
	def calcSpeed( self, vCrossbridge=1.0 ):
		'''Calculate the theoretical maximum speed of shortening as the 
		product of the number of half sarcomeres in series along the 
		length of the muscle cylinder and the speed of crossbridge 
		cycling, vCrossbridge (in um/s).'''
		self.v = self.calcNhs() * vCrossbridge * 1e-4
		return self.v
		
	# Mutators
		
	def set_lthick( self, lthick ):
		'''Set thick filament length to new lthick, in um'''
		self.lthick = lthick
		
	def set_lthin( self, lthin ):
		'''Set thin filament length to new lthin, in um'''
		self.lthin = lthin
		
	def set_lz( self, lz ):
		'''Set z disc width to new lz, in um'''
		self.lz = lz
		
	def set_lbare( self, lbare ):
		'''Set bare zone width to new lbare, in um'''
		self.lbare = lbare
	
	def set_t0( self, t0 ):
		'''Set peak isometric tension to new t0, in N/cm**2'''
		self.t0 = t0
	
	def set_rCyl( self, rCyl ):
		'''Set radius of the muscle cylinder to new rCyl, in cm'''
		self.rCyl = rCyl
	
	def set_lCyl( self, lCyl ):
		'''Set length of the muscle cylinder to new lCyl, in cm'''
		self.lCyl = lCyl
	
	def set_csv( self, csv ):
		'''Set emprical data file: a CSV in which column 0 contains
		isometric tetanus lengths (in um) and column 1 contains the
		corresponding active tensions (in N/cm**2 or %T0).'''
		self.csv = csv
			
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
		
	def get_csv( self ):
		'''Return the name of the empirical data file'''
		return self.csv
		
	# Utilities
	
	def print( self ):
		'''Describe model characteristics with string formatting.'''
		stateStr  = "\n------------------------------------------------"
		stateStr += "\n\t\tmuscle performance"
		stateStr += "\n------------------------------------------------"
		stateStr += "\nL thick\t{0:<5.3f} um".format(self.lthick)
		stateStr += "\nL thin\t{0:<5.3f} um".format(self.lthin)
		stateStr += "\nL lz\t{0:<5.3f} um".format(self.lz)
		stateStr += "\nL bare\t{0:<5.3f} um".format(self.lbare)
		stateStr += "\nL cyl\t{0:<5.3f} cm".format(self.lCyl)
		stateStr += "\nr cyl\t{0:<5.3f} cm".format(self.rCyl)
		stateStr += "\nCSA  \t{0:<5.3f} cm**2".format(np.pi * self.rCyl**2)
		stateStr += "\nF0\t{0:<5.3f} N".format(self.f0)
		stateStr += "\nNhs\t{0:<7.1f} half sarcomeres".format(self.Ns)
		stateStr += "\n\nWork\t{0:<5.3f} J".format(self.w)
		stateStr += "\n\n\tAscending  {0:<5.3f} J".format(self.wAsc3 + self.wAsc2)
		stateStr += "\n\tPlateau    {0:<5.3f} J".format(self.wPlat + self.wAsc1)
		stateStr += "\n\tDescending {0:<5.3f} J".format(self.wDesc)
		stateStr += "\n\nEnergy Density\t{0:<6.3f} J/cm**3".format(self.energyDensity)
		stateStr += "\n"
		print( stateStr )
		
	def isStable( self ):
		'''Decide whether the model is stable (True) or whether it has 
		gone irrevocably awry (False). Failure is defined as a tension
		in excess of 1.8*T0.'''
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
		plt.xlabel("Length (cm)")
		plt.ylabel("Force (N)")
		plt.title("Muscle model\nr = {0:3.2f} cm, L = {1:3.2f} cm".format(self.rCyl, self.lCyl))	
		
		if self.csv != None:
			# If empirical data is normalized, un-normalize it
			# (Keep in mind that sometimes re-digitized normalized measurements stray above 1.0)
			l = self.sarcomere.L_empirical * 10**-6 * self.Ns
			f = self.sarcomere.Ta_empirical * self.csa
			plt.plot(l, f, "ob", markerfacecolor='w')
				
		px = self.sarcomere.L  * 10**-6 * self.Ns
		py = self.sarcomere.Ta * self.csa
		plt.plot(px, py, '-dk', markerfacecolor='r')
		
		plt.show()
		
		return fig



	
def superimposeAnimalData( thick, thin, slope, intercept, fig=None, quadratic=False, poly2=None ):
	'''Superimpose animal data on the response surface and measure
	root mean squared error and R**2 relative the identified trend.'''	
		
	data = np.genfromtxt("animal_data.csv", delimiter=',', autostrip=True, skip_header=1)
	
	# Keep only data that fits within in the bounds of simulated sarcomere dimensions
	data = data[data[:,0] <= max(thick)]
	data = data[data[:,0] >= min(thick)]
	data = data[data[:,1] <= max(thin)]
	data = data[data[:,1] >= min(thin)]
	
	animalThick = data[:,0]
	animalThin  = data[:,1]
	
	# Animal data relative to model's prediction of ultrastructure for optimal performance
	if( not quadratic ):
		modelThin   = animalThick*slope + intercept	# where our regression expects to see animals' thin filaments
	else:
		modelThin 	= poly2[0]*np.square(animalThick) + poly2[1]*animalThick + poly2[0]
	rmseModel   = (np.sum(np.mean(np.square(animalThin - modelThin))))**0.5
	animalMean  = np.mean(animalThin)
	ss_tot      = np.sum(np.square(animalThin - animalMean))
	ss_res      = np.sum(np.square(animalThin - modelThin))
	rsqModel    = 1 - ss_res/ss_tot
	
	# Trend of the animal data in isolation (quadratic shape identified through offline analysis)
	poly2 = np.polyfit( animalThick, animalThin, 2 )
	poly2Thin	= poly2[0]*np.square(animalThick) + poly2[1]*animalThick + poly2[0]
	ss_res      = np.sum(np.square(animalThin - poly2Thin))
	rsqPoly2    = 1 - ss_res/ss_tot
	poly2Thin	= poly2[0]*np.square(thick) + poly2[1]*thick + poly2[0]
	poly2Thin	= poly2Thin[poly2Thin <= np.max(thin)]
	poly2Thick	= thick[0:max(poly2Thin.shape)]
	if( fig!= None ):
		plt.plot( animalThick, animalThin, 'or', label="animal data, R**2={0:5.3f} wrt peak performance".format(rsqModel) )
		plt.plot( poly2Thick, poly2Thin, '-r', label="thin={0:5.3f}*thick**2 + {1:5.3f}*thick + {2:5.3f}, R**2={3:5.3f}".format(poly2[0], poly2[1], poly2[2], rsqPoly2))
	
	return (rsqModel, rmseModel)
	
def trendLine( thick, thin, surface, fig=None, quadratic=False ):
	'''Find the thin filament length that optimizes performance for each
	thick filament length, and perform a regression to characterize the 
	trend.'''	
	
	# Find the thin filament length, a[i], that maximize performance for
	# each thick filament length, thick[i].
	nThick = max(thick.shape)
	nThin = max(thin.shape)
	i = np.argmax(surface, axis=0)
	a = thin[i]
	m = thick[a < max(thin)]
	a = a[a < max(thin)]
	m = m[ min(thin) < a ]
	a = a[ min(thin) < a ]
	
	if max(a.shape) < 1:
		i = np.argmax(np.flipud(surface), axis=0)
		a = thin[nThin-1-i]
		m = thick[a < max(thin)]
		a = a[a < max(thin)]
		m = m[ min(thin) < a ]
		a = a[ min(thin) < a ]
	
	if( not quadratic ):
		# Perform a least squares regression on optimal (thick, a) pairs
		slope, intercept, r_value, p_value, std_err = stats.linregress(m, a)
		rmse = (np.sum(np.mean(np.square(a - slope*m+intercept))))**0.5
		rsqLinear   = r_value**2
			
		# Plot the regression
		if( fig != None ):
			plt.plot( m, a, '.c', markersize=1, label="peak performance" )
			plt.plot( m, slope*m+intercept, 'c', label="thin=thick*{0:5.3f}+{1:5.3f}, R**2={2:5.3f}".format(slope, intercept, rsqLinear))
			
		return (slope, intercept, rmse)
	else:
		# Performing 2nd order polynomial regressions on optimal (thick, a) 
		# pairs resulted in weaker correlations across all response surfaces
		# in classic model
		poly2 = np.polyfit( m, a, 2 )
		poly2Thin	= poly2[0]*np.square(m) + poly2[1]*m + poly2[0]
		rmse = (np.sum(np.mean(np.square(a - poly2Thin))))**0.5
		poly2Thin	= poly2Thin[poly2Thin <= np.max(thin)]
		poly2Thick	= m[0:np.max(poly2Thin.shape)]
		# Plot the regression
		if( fig != None ):
			plt.plot( m, a, '.c', markersize=1, label="peak performance" )
			plt.plot( poly2Thick, poly2Thin, 'c', label="thin={0:5.3f}*thick**2 + {1:5.3f}*thick + {2:5.3f}".format(poly2[0], poly2[1], poly2[2]))
	
		return (poly2, rmse)
	

def responseSurfaces():
	'''Generating response surfaces by simulating sliding filament 
	models for a range of thick and thin filament combinations.'''
	
	# Gordon's 1966 frog semitendinosus descending limb slope
	# slopeGordon1966 = [ 53.544309, 10.198916, 0, 0, -18.21235]	# N/cm**2 per micrometer overlap
	
	# Median vertebrate descending limb slope (Heejoon's reported value)
	slopeVertebrate = [67.2000, 6.1935, 0, 0, -17.3913]
	
	# Thick and thin filament lengths to simulate, in um
	# Cross-phyla range
	thick = np.linspace(0.1, 20.1, num=201)
	thin = np.linspace(0.1, 20.1, num=201)
	# Vertebrate range
	#thick = np.linspace(1.25, 2.0, num=176*2)
	#thin = np.linspace(1.5, 3.0, num=151*2)
	
	# 2D arrays to aggregate performance metrics
	nThick = max(thick.shape)
	nThin = max(thin.shape)
	force = np.zeros( (nThin, nThick) )
	work = np.zeros( (nThin, nThick) )
	ed = np.zeros( (nThin, nThick) )
	nhs = np.zeros( (nThin, nThick) )
	X, Y = np.meshgrid(thick, thin)
	
	# Figure in which to aggregate muscle models
	fig = plt.figure()
	plt.xlabel("Length (cm)")
	plt.ylabel("Force (N)")
	plt.title("Anchored muscle model")	
	
	
	for m in range(nThick):
		for a in range(nThin):
			
			# Generate a sliding filament model for this new combination of filament lengths
			muscle = Homogeneous( lthick=thick[m], lthin=thin[a], slope=slopeVertebrate )
			
			# Add it to the plot of aggregated sliding filament models
			if random.random() < 0.00025:
				plt.plot(muscle.L, muscle.F, 'd-', markerfacecolor='w')
				plt.pause(0.01)
			
			# Peak active isometric force
			force[a][m] = muscle.f0
			
			# Theoretical max work (J)
			work[a][m] = muscle.w
			
			# Energy density = theoretical max work / volume
			ed[a][m] = muscle.energyDensity
			
			# Number of half sarcomeres in series correlates with shortening velocity
			nhs[a][m] = muscle.Nhs
	
	# Peak tension response surface
	fig = plt.figure()
	plt.xlabel("Thick filament (um)")
	plt.ylabel("Thin filament (um)")	
	plt.title("Anchored model's predicted\npeak force (N)")
	contourF = plt.contourf(X, Y, force)
	plt.colorbar(contourF)
	#slope, intercept, rmse = trendLine( thick, thin, force, fig )
	slope, intercept, r_value, p_value, std_err = stats.linregress(thick, thick+Homogeneous.lz)
	rmse = (np.sum(np.mean(np.square(thick+Homogeneous.lz - slope*thick+intercept))))**0.5
	rsq   = r_value**2
	plt.plot( thick, thick+Homogeneous.lz, '.c', markersize=1, label="peak performance" )
	plt.plot( thick, slope*thick+intercept, 'c', label="thin=thick*{0:5.3f}+{1:5.3f}, R**2={2:5.3f}".format(slope, intercept, rsq))
	rsq, rmse = superimposeAnimalData( thick, thin, slope, intercept, fig )
	plt.xlim((np.min(thick), np.max(thick)))
	plt.ylim((np.min(thin), np.max(thin)))
	plt.legend()
	
	# Specific work response surface
	fig = plt.figure()
	plt.xlabel("Thick filament (um)")
	plt.ylabel("Thin filament (um)")	
	plt.title("Anchored model's predicted\ntheoretical max work (J)")
	contourW = plt.contourf(X, Y, work)
	plt.colorbar(contourW)
	slope, intercept, rmse = trendLine( thick, thin, work, fig )
	rsq, rmse = superimposeAnimalData( thick, thin, slope, intercept, fig )
	plt.xlim((np.min(thick), np.max(thick)))
	plt.ylim((np.min(thin), np.max(thin)))
	plt.legend()
	
	# Specific work response surface
	fig = plt.figure()
	plt.xlabel("Thick filament (um)")
	plt.ylabel("Thin filament (um)")	
	plt.title("Anchored model's predicted\nenergy density (J/cm**3)")
	contourE = plt.contourf(X, Y, ed)
	plt.colorbar(contourE)
	slope, intercept, rmse = trendLine( thick, thin, ed, fig )
	rsq, rmse = superimposeAnimalData( thick, thin, slope, intercept, fig )
	plt.xlim((np.min(thick), np.max(thick)))
	plt.ylim((np.min(thin), np.max(thin)))
	plt.legend()
	
	# Number of half sarcomeres in series (NHS) response surface
	# provides a heuristic for maximal velocity
	fig = plt.figure()
	plt.xlabel("Thick filament (um)")
	plt.ylabel("Thin filament (um)")	
	plt.title("Anchored model's predicted\nnumber of half sarcomeres in series")
	contourN = plt.contourf(X, Y, nhs)
	plt.colorbar(contourN)
	slope, intercept, rmse = trendLine( thick, thin, nhs, fig )
	rsq, rmse = superimposeAnimalData( thick, thin, slope, intercept, fig )
	plt.xlim((np.min(thick), np.max(thick)))
	plt.ylim((np.min(thin), np.max(thin)))
	plt.legend()
	
	# Plot composite of normalized surfaces
	fig = plt.figure()
	plt.xlabel("Thick filament (um)")
	plt.ylabel("Thin filament (um)")	
	plt.title("Anchored model's normalized\nperformance metrics superimposed")
	composite = (force/force.max() + 2*nhs/nhs.max() + ed/ed.max())/4
	contourComp = plt.contourf(X, Y, composite )
	plt.colorbar( contourComp )
	poly2, rmse = trendLine( thick, thin, composite, fig, True )
	rsq, rmse = superimposeAnimalData( thick, thin, slope, intercept, fig, True, poly2 )
	plt.xlim((np.min(thick), np.max(thick)))
	plt.ylim((np.min(thin), np.max(thin)))
	plt.legend()
	
	
	plt.show()

def main( argv ):
	'''Test the simulation by creating a single sliding filament model
	and analyzing its characteristics.'''
	
	# Check to see if we should run the responseSurfaces() function
	if len(argv) > 1 and argv[1].strip() == "response":
		responseSurfaces()
		return
		
	# Construct new SlidingFilament model with default characteristics
	muscle = Homogeneous()
	
	print("\nUSAGE: python3 muscleAnchored.py [optional: -flag <value> pairs]\n")
	print("FLAGS:")
	print("\t-m <thick filament length, um>, default: {} um".format(Homogeneous.lthick))
	print("\t-a <thin filament length, um>, default: {} um".format(Homogeneous.lthin))
	print("\t-z <z disc width, um>, default: {} um".format(Homogeneous.lz))
	print("\t-b <bare zone width, um>, default: {} um".format(Homogeneous.lbare))
	print("\t-r <radius of muscle cylinder>, default: {} cm".format(Homogeneous.rCyl))
	print("\t-l <length of muscle cylinder>, default: {} cm".format(Homogeneous.lCyl))
	print("\t-f <empirical data file, *.csv>, default: {}".format(Homogeneous.csv))

	# Check commandline arguments for flags. Use flags to update model 
	# characteristics before simulation. 
	flagMap = { "-m":muscle.set_lthick,
				"-a":muscle.set_lthin,
				"-z":muscle.set_lz,
				"-b":muscle.set_lbare,
				"-t":muscle.set_t0,
				"-r":muscle.set_rCyl,
				"-l":muscle.set_lCyl,
				"-f":muscle.set_csv }
	for i in range(1,len(argv)-1):
		flag = argv[i]
		if flag in flagMap.keys():
			flagFunc = flagMap[ flag ]
			if flag != "-f":
				flagParam = float(argv[i+1])
			else:
				flagParam = argv[i+1]
			flagFunc( flagParam )
	
	muscle.init_sarcomere_anchored()
	muscle.calcWork()
	muscle.calcSpeed()
	
	# Print and plot sarcomere geometry info and muscle performance
	muscle.sarcomere.print()
	muscle.print()
	muscle.sarcomere.plot()
	muscle.plot()
	plt.show()
	
	
if __name__=="__main__":
	main( sys.argv )
		

