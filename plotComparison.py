# comparing different contributions for babyIAXO

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(3e-4,1e0), ylim=(1e-13,3e-6))

# original
dat = loadtxt("data/limits/babyIAXO-tPlasmon-clean.dat")
ax2.plot(dat[:,0],dat[:,1], color='black', label='T-plasmon')

# 1day gas original
dat = loadtxt("data/limits/babyIAXO-tPlasmon-clean-gas.dat")
ax2.plot(dat[:,0],dat[:,1], color='green', label='T-plasmon (gas)')

# L-plasmon new
dat = loadtxt("data/limits/babyIAXO-newB-lMixingRes.dat")
ax2.plot(dat[:,0], dat[:,1], color='red', label='L-plasmon mixing')

# L-plasmon new gas
dat = loadtxt("data/limits/babyIAXO-newB-lMixingResGas.dat")
ax2.plot(dat[:,0], dat[:,1], color='magenta', label='L-plasmon mixing (gas)')

# pure L contribution
dat = loadtxt("data/limits/babyIAXO-clean-pureL.dat")
ax2.plot(dat[:,0], dat[:,1], color='cyan', label='pure L conversion')

# axes
ax2.set_xlabel("Dark photon mass [eV]")
ax2.set_ylabel("Kinetic mixing")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/comparison-clean.jpg')
plt.show()
