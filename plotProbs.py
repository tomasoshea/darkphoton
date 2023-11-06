# comparing different contributions for babyIAXO

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(1e-4,1e4), ylim=(1e-8,1e9))
alp = 0.6
lw = 4



############## Atlas data #################

# vacuum
dat = loadtxt("data/probVac-50eV.dat")
ax2.plot(dat[:,0],dat[:,1], label='vacuum (50 eV)', color='black', ls='-', lw=lw, alpha=alp)

dat = loadtxt("data/probVac-500eV.dat")
ax2.plot(dat[:,0],dat[:,1], label='vacuum (500 eV)', color='black', ls=':', lw=lw, alpha=alp)

dat = loadtxt("data/probVac-5keV.dat")
ax2.plot(dat[:,0],dat[:,1], label='vacuum (5 keV)', color='black', ls='--', lw=lw, alpha=alp)


# gas
dat = loadtxt("data/probGas-50eV.dat")
ax2.plot(dat[:,0],dat[:,1], label='gas (50 eV)', color='red', ls='-', lw=lw, alpha=alp)

dat = loadtxt("data/probGas-500eV.dat")
ax2.plot(dat[:,0],dat[:,1], label='gas (500 eV)', color='red', ls=':', lw=lw, alpha=alp)

dat = loadtxt("data/probGas-5keV.dat")
ax2.plot(dat[:,0],dat[:,1], label='gas (5 keV)', color='red', ls='--', lw=lw, alpha=alp)


# axes
ax2.set_xlabel("Dark photon mass [eV]")
ax2.set_ylabel("Back-conversion probablity")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/backConversion.jpg')
plt.show()

