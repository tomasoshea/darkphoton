# comparing different contributions for babyIAXO

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(1e-3,1e4), ylim=(1e-15,1e-7))


# T-plasmon Atlas
dat = loadtxt("data/limits/stats-babyIAXO-Atlas-100eV.dat")
ax2.plot(dat[:,0],dat[:,1], label='old (100 eV)', color='red',ls='--', lw=2)

dat = loadtxt("data/limits/statsAtlas-babyIAXO-Atlas.dat")
ax2.plot(dat[:,0],dat[:,1], label='new (100 eV)', color='green',ls=':', lw=2)

dat = loadtxt("data/limits/stats-babyIAXO-Atlas-1eV.dat")
ax2.plot(dat[:,0],dat[:,1], label='old (1 eV)', color='red',ls='--', lw=4)

dat = loadtxt("data/limits/statsAtlas-babyIAXO-1eV.dat")
ax2.plot(dat[:,0],dat[:,1], label='new (1 eV)', color='green',ls=':', lw=4)


# axes
#ax2.set_xlabel("Dark photon mass [eV]")
#ax2.set_ylabel("Flux (normalised)")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()#loc='upper right')

plt.savefig('plots/limit-comparison.jpg')
plt.show()

