# plot from .dat files produced by C++ integration

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

# constants
T = 1.334824922283800106e+03	# solar temp at core
amu = 931.49410242e6			# atomic mass unit in eV


plt.style.use("style.txt")	# import plot style

# setup plot
start = 1e0
stop = 1e7
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set( xlim=(start, stop))#, ylim=(-0.01, 1.01) )
#ax2.set( ylim = (1e0, 1e14) )
#ax2.set( xlim=(0,3e3))#, ylim=(1e-20,1e-8))

def gaussian(x, w, phi, T, m):
	s2 = w*w*T/m	# sigma squared
	g = phi * np.exp(-(x-w)*(x-w) / s2*s2) / np.sqrt(2*np.pi*s2)
	print(g)
	print("s2 = " + str(s2))
	return g


# get data and plot

# T plasmon
#datT = loadtxt("data/Espectrum0.dat")	# m = 1 eV
#phiT = datT[:,1] / np.nanmax(datT[:,1])	# normalise phi
#ax2.plot(datT[:,0], phiT, ls = '--', label = 'T plasmon m = 1 eV')

#datT = loadtxt("data/Espectrum1.dat")	# m = 10 eV
#phiT = datT[:,1] / np.nanmax(datT[:,1])	# normalise phi
#ax2.plot(datT[:,0], phiT, ls = '--', label = 'T plasmon m = 10 eV')

datT = loadtxt("data/Espectrum2.dat")	# m = 1e2 eV
phiT = datT[:,1]	# normalise phi
ax2.plot(datT[:,0], phiT, ls = '--', label = 'T plasmon m = 100 eV')

#datT = loadtxt("data/Espectrum3.dat")	# m = 1e2 eV
#phiT = datT[:,1] / np.nanmax(datT[:,1])	# normalise phi
#ax2.plot(datT[:,0], phiT, ls = '--', label = 'T plasmon m = 1 keV')


# pp-chain
ax2.vlines(5.49e6, 0, 1, label = 'deuterium fusion')	# deuterium fusion
ax2.vlines(5.11e5, 0, 1, label = 'e+ e- annihilation', color='red')	# e+ e- annihilation

# 57Fe
w = 14.4e3				# 14.4 keV
m = 56.9353928 * amu	# mass of nucleus [eV]
phi = 85.4099			# taken from babyIAXO-57Fe-3.dat
x = np.linspace(w - 5, w+5 ,1000)
ax2.plot(x, gaussian(x, w, phi, T, m), label = '57Fe')

# p+D -> 3He + g
w = 5.49e6				# 14.4 keV
m = 1.5 * amu	# mass of nucleus [eV]
phi = 9.28108e-13			# taken from babyIAXO-57Fe-3.dat
wid = 28793710.044686843
x = np.linspace(w - 5, w+5 ,1000)
ax2.plot(x, gaussian(x, w, phi, T, m), label = 'pp')



# axes
ax2.set_xlabel("Dark photon energy [eV]")
ax2.set_ylabel("Flux (normalised)")
#ax2.set_ylabel("Differential flux [eV2]")
ax2.set_xscale('log')
#ax2.set_yscale('log')
ax2.tick_params(axis='both')
ax2.legend()

fig2.savefig('plots/spectrum_lines-1.jpg')

plt.show()

