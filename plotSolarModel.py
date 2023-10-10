# comparing different contributions for babyIAXO

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(0.,1.01), ylim=(1e-3,1.01))

# radius (fraction)
r = loadtxt("data/rFrac.dat")

# plasma freq
dat = loadtxt("data/wp.dat")
dat = dat / np.nanmax(dat)
ax2.plot(r,dat, label='Plasma frequency')

# temperature
dat = loadtxt("data/T.dat")
dat = dat / np.nanmax(dat)
ax2.plot(r,dat,label="Temperature")

# electron number density
dat = loadtxt("data/ne.dat")
dat = dat / np.nanmax(dat)
ax2.plot(r,dat,label="Electron number density")

# H number density
dat = loadtxt("data/nH.dat")
dat = dat / np.nanmax(dat)
ax2.plot(r,dat,label="H number density")

# 3He number density
dat = loadtxt("data/nHe3.dat")
dat = dat / np.nanmax(dat)
ax2.plot(r,dat,label="3He number density")

# 4He number density
dat = loadtxt("data/nHe4.dat")
dat = dat / np.nanmax(dat)
ax2.plot(r,dat,label="4He number density")

# 57Fe number density
dat = loadtxt("data/n57Fe.dat")
dat = dat / np.nanmax(dat)
ax2.plot(r,dat,label="57Fe number density")

# axes
ax2.set_xlabel("Solar radius")
#ax2.set_ylabel("Plasma frequency [eV]")
#ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/solarmodel-log.jpg')
plt.show()
