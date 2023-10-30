# Tom O'Shea 2023

# plot Gammas

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(1,10), ylim=(1e-7,1e1))

dat = loadtxt("data/GammaCAST08.dat")
ax2.plot(dat[:,0],dat[:,1],color='black',ls='--',label="CAST 0.08 mbar")
dat = loadtxt("data/GammaCAST13.dat")
ax2.plot(dat[:,0],dat[:,1],color='black',ls=':',label="CAST 13.43 mbar")
dat = loadtxt("data/GammaIAXO08.dat")
ax2.plot(dat[:,0],dat[:,1],color='red',ls='--',label="IAXO 0.08 mbar")
dat = loadtxt("data/GammaIAXO13.dat")
ax2.plot(dat[:,0],dat[:,1],color='red',ls=':',label="IAXO 13.43 mbar")

# axes
ax2.set_xlabel("Energy [keV]")
ax2.set_ylabel("Absorption Gamma [m-1]")
#ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

wtab = 10**np.arange(0,4.0002,0.001)
mtab = 10**np.arange(-3,4.1,0.1)

np.savetxt('data/wtab.dat', wtab, delimiter=',')
np.savetxt('data/mtab.dat', mtab, delimiter=',')

plt.savefig('plots/GammaCAST-IAXO.jpg')
plt.show()
