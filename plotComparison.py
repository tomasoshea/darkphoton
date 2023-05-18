# comparing different contributions for babyIAXO

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(1e-3,1e0), ylim=(1e-12,3e-1))
#ax2.set(xlim=(1e-80,1e0), ylim=(1e-100,3e-6))
#ax2.set(xlim=(1e-6,1e5), ylim=(1e-13,3e-4))

# tPlasmon
dat = loadtxt("data/limits/babyIAXO-tPlasmon-clean.dat")
ax2.plot(dat[:,0],dat[:,1], color='black', label='T-plasmon')

# tPlasmon (newE comparison)
#dat = loadtxt("data/limits/babyIAXO-tPlasmon-newE.dat")
#ax2.plot(dat[:,0],dat[:,1], color='black',  label='T-plasmon')

# tPlasmon gas
dat = loadtxt("data/limits/babyIAXO-tPlasmon-clean-gas.dat")
ax2.plot(dat[:,0],dat[:,1], color='green', label='T-plasmon (gas)')

# tPlasmon gas new phi
dat = loadtxt("data/limits/babyIAXO-tPlasmon-newphi-gas.dat")
ax2.plot(dat[:,0],dat[:,1], color='green', ls=':', label='T-plasmon (gas)')

# tPlasmon new 95%
dat = loadtxt("data/stats-babyIAXO.dat")
ax2.plot(dat[:,0],dat[:,1], color='red', label='T-plasmon (gas) CL95')

# tPlasmon new 95%
dat = loadtxt("data/limits/stats-babyIAXO2.dat")
ax2.plot(dat[:,0],dat[:,1], color='red', ls=':', label='T-plasmon (gas) (2) CL95')

# tPlasmon gas newE
#dat = loadtxt("data/limits/babyIAXO-tPlasmon-100eV-gas.dat")
#ax2.plot(dat[:,0],dat[:,1], color='green',  label='T-plasmon (gas)')

# L-plasmon
#dat = loadtxt("data/limits/babyIAXO-clean-lMixingRes.dat")
#ax2.plot(dat[:,0], dat[:,1], color='red', label='L-plasmon mixing')

# L-plasmon newE
#dat = loadtxt("data/limits/babyIAXO-newE-lMixingRes.dat")
#ax2.plot(dat[:,0], dat[:,1], color='red',  label='L-plasmon mixing')

# L-plasmon gas
dat = loadtxt("data/limits/babyIAXO-clean-lMixingResGas.dat")
ax2.plot(dat[:,0], dat[:,1], color='magenta', label='L-plasmon mixing (gas)')

# L-plasmon new gas newE
#dat = loadtxt("data/limits/babyIAXO-newerE-lMixingResGas.dat")
#ax2.plot(dat[:,0], dat[:,1], color='magenta', ls=':')#, label='L-plasmon mixing (gas) (newE)')

# L-plasmon new gas newE newP
#dat = loadtxt("data/limits/babyIAXO-newP-lMixingResGas.dat")
#ax2.plot(dat[:,0], dat[:,1], color='magenta',  label='L-plasmon mixing (gas)')

# pure L contribution
#dat = loadtxt("data/limits/babyIAXO-conversion-pureL.dat")
#ax2.plot(dat[:,0], dat[:,1], color='cyan', label='IAXO B-field conversion')

# axes
ax2.set_xlabel("Dark photon mass [eV]")
ax2.set_ylabel("Kinetic mixing")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/comparison-cl95.jpg')
plt.show()
