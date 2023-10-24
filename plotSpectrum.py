# plot from .dat files produced by C++ integration

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set( xlim=(1e0, 1e4), ylim=(-0.01, 1.01) )
#ax2.set( ylim = (1e0, 1e14) )
#ax2.set( xlim=(0,3e3))#, ylim=(1e-20,1e-8))


# get data and plot

# T plasmon
datT = loadtxt("data/Espectrum-1.dat")	# m = 1 eV
t1 = np.nanmax(datT[:,1])
phiT = datT[:,1] / np.nanmax(datT[:,1])	# normalise phi
ax2.plot(datT[:,0], phiT,lw=5,ls='-',color='black', label = 'T plasmon m = 0.1 eV')

datT = loadtxt("data/Espectrum2.dat")	# m = 1e2 eV
t2 = np.nanmax(datT[:,1])
phiT = datT[:,1] / np.nanmax(datT[:,1])	# normalise phi
ax2.plot(datT[:,0], phiT,lw=2,ls='-',color='black', label = 'T plasmon m = 100 eV')

print(t1/t2)


# lMixing
dat = loadtxt("data/Espectrum-lMixing-1.dat")
phi = dat[:,1] / np.nanmax(dat[:,1])
ax2.plot(dat[:,0], phi,ls=':',color='black',label = 'L mixing')


# pureL
dat = loadtxt("data/Espectrum-pureL-1.dat")
phi = dat[:,1] / np.nanmax(dat[:,1])
ax2.plot(dat[:,0],phi,ls='--',color='black',label = 'pure L')



# axes
ax2.set_xlabel("Dark photon energy [eV]")
ax2.set_ylabel("Flux (normalised)")
#ax2.set_ylabel("Differential flux [eV2]")
ax2.set_xscale('log')
#ax2.set_yscale('log')
ax2.tick_params(axis='both')
ax2.legend()

fig2.savefig('plots/spectrum-black.jpg')

plt.show()

