import numpy as np
import matplotlib.pyplot as plt

# Read in dataset of published sarcomere geometries
data = np.genfromtxt("geometries.csv", delimiter=',', autostrip=True, skip_header=1)
animalThick = data[:,0]
animalThin  = data[:,1]

# Scatter plot dataset of published sarcomere geometries
fig = plt.figure()
plt.xlabel("Thick Filament (um)")
plt.ylabel("Thin Filament (um)")
plt.title("Known Sarcomere Geometries")	
plt.plot( animalThick, animalThin, 'or' )

# Trend of the animal data in isolation (quadratic shape identified through offline analysis)
poly2 = np.polyfit( animalThick, animalThin, 2 )
poly2Thin	= poly2[0]*np.square(animalThick) + poly2[1]*animalThick + poly2[0]
ss_res      = np.sum(np.square(animalThin - poly2Thin))
animalMean  = np.mean(animalThin)
ss_tot      = np.sum(np.square(animalThin - animalMean))
rsqPoly2    = 1 - ss_res/ss_tot
thick = np.linspace(0.1, 10.1, num=801)
thin = np.linspace(0.1, 20.1, num=801)
poly2Thin	= poly2[0]*np.square(thick) + poly2[1]*thick + poly2[0]
poly2Thin	= poly2Thin[poly2Thin <= np.max(thin)]
poly2Thick	= thick[0:max(poly2Thin.shape)]
plt.plot( poly2Thick, poly2Thin, '-r', label="thin={0:5.3f}*thick**2 + {1:5.3f}*thick + {2:5.3f}, R**2={3:5.3f}".format(poly2[0], poly2[1], poly2[2], rsqPoly2))

# Display the figure window
plt.show()
