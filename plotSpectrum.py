# plot from .dat files produced by C++ integration

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set( xlim=(1e0, 1e3), ylim=(-0.01, 1.01) )


# get data and plot

# T plasmon
datT = loadtxt("data/Espectrum-1.dat")	# m = 1e-1 eV
phiT = datT[:,1] / np.nanmax(datT[:,1])	# normalise phi
ax2.plot(datT[:,0], phiT, ls = '--', color=(1., 0.8, 0.), label = 'T plasmon m = 0.1 eV')

datT = loadtxt("data/Espectrum0.dat")	# m = 1 eV
phiT = datT[:,1] / np.nanmax(datT[:,1])	# normalise phi
ax2.plot(datT[:,0], phiT, ls = '--', color=(1., 0.4, 0.), label = 'T plasmon m = 1 eV')

datT = loadtxt("data/Espectrum1.dat")	# m = 10 eV
phiT = datT[:,1] / np.nanmax(datT[:,1])	# normalise phi
ax2.plot(datT[:,0], phiT, ls = '--', color=(1., 0., 0.), label = 'T plasmon m = 10 eV')


# lMixing
dat = loadtxt("data/Espectrum-lMixing-1.dat")
phi = dat[:,1] / np.nanmax(dat[:,1])
ax2.plot(dat[:,0], phi, color="cyan", label = 'L mixing')


# pureL
dat = loadtxt("data/Espectrum-pureL-1.dat")
phi = dat[:,1] / np.nanmax(dat[:,1])
ax2.plot(dat[:,0], phi, color="magenta", label = 'pure L')	


# axes
ax2.set_xlabel("Plasma frequency [eV]")
ax2.set_ylabel("Flux (normalised)")
ax2.set_xscale('log')
ax2.tick_params(axis='both')
ax2.legend()

fig2.savefig('plots/spectrum-pureL.jpg')

plt.show()

