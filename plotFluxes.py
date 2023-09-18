# comparing different contributions for babyIAXO

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(1e-3,1e7), ylim=(1e-35,2e0))


# tPlasmon
dat = loadtxt("data/limits/babyIAXO-tPlasmon-flux.dat")
top = np.nanmax(dat[:,1])
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='T-plasmon')

#dat = loadtxt("data/limits/babyIAXO-tPlasmon-70eV.dat")
#dat[:,1] = dat[:,1] / top
#ax2.plot(dat[:,0],dat[:,1], ls=':', label='T-plasmon (70 eV)')

# tPlasmon gas
dat = loadtxt("data/limits/babyIAXO-tPlasmon-100eV-gas.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='T-plasmon (gas)')

# tPlasmon gas T=3K
#dat = loadtxt("data/limits/babyIAXO-tPlasmon-100eV-3K-gas.dat")
#dat[:,1] = dat[:,1] / top
#ax2.plot(dat[:,0],dat[:,1], ls=":", label='T-plasmon (cold gas)')

# lMixing
dat = loadtxt("data/limits/babyIAXO-1-lMixingRes-flux.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='L-plasmon mixing')

# lMixing gas
dat = loadtxt("data/limits/babyIAXO-1-lMixingResGas-flux.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='L-plasmon mixing (gas)')

# pp chain
dat = loadtxt("data/limits/babyIAXO-pp-2.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='pp-chain')

# prob
#fig = plt.figure(2)	# display is 1920 x 1080 (16:9)
#ax = fig.add_axes((.1,.1,.8,.8))
#ax.set(xlim=(1e-3,1e7), ylim=(0.,1.1))
#dat = loadtxt("data/highEprob.dat")
#ax.plot(dat[:,0],dat[:,1], label='conversion prob in vacuum')
#ax.set_xscale('log')

# axes
ax2.set_xlabel("Dark photon mass [eV]")
ax2.set_ylabel("Flux (normalised)")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/flux-comparison-pp2.jpg')
plt.show()
