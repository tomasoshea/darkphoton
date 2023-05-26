# comparing different contributions for babyIAXO

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(1e-3,1e3), ylim=(1e-20,2e0))


# tPlasmon
dat = loadtxt("data/limits/babyIAXO-tPlasmon-flux.dat")
top = np.nanmax(dat[:,1])
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], color='red', label='T-plasmon')

# tPlasmon gas
dat = loadtxt("data/limits/babyIAXO-tPlasmon-100eV-gas.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], color='magenta', label='T-plasmon (gas)')

# lMixing
dat = loadtxt("data/limits/babyIAXO-1-lMixingRes-flux.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], color='blue', label='L-plasmon mixing')

# lMixing gas
dat = loadtxt("data/limits/babyIAXO-1-lMixingResGas-flux.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], color='cyan', label='L-plasmon mixing (gas)')

# axes
ax2.set_xlabel("Dark photon mass [eV]")
ax2.set_ylabel("Flux (normalised)")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/flux-comparison1.jpg')
plt.show()
