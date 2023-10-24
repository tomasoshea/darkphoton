# comparing different contributions for babyIAXO

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

m2eV = (1.973269804e-7)

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(0.,1.01), ylim=(1e-3,1.01))

# radius (fraction)
r = loadtxt("data/rFrac.dat")

# plasma freq
dat = loadtxt("data/wp.dat")
print("plasma freq:	{} eV".format(np.nanmax(dat)))
dat = dat / np.nanmax(dat)
ax2.plot(r,dat,ls='--',color='black',label='Plasma frequency')

# temperature
dat = loadtxt("data/T.dat")
print("temperature:	{}	eV".format(np.nanmax(dat)))
dat = dat / np.nanmax(dat)
ax2.plot(r,dat,ls=':',color='black',label="Temperature")

# H number density
dat = loadtxt("data/nH.dat")
print("nH:		{} eV3".format(np.nanmax(dat)))
print("nH:		{} m-3".format(np.nanmax(dat)/(m2eV**3)))
dat = dat / np.nanmax(dat)
#ax2.plot(r,dat,ls=(0, (3, 5, 1, 5, 1, 5)),color='black',label="H number density")

# 57Fe number density
dat = loadtxt("data/n57Fe.dat")
print("n57Fe:		{} eV3".format(np.nanmax(dat)))
print("n57Fe:		{} m-3".format(np.nanmax(dat)/(m2eV**3)))
dat = dat / np.nanmax(dat)
ax2.plot(r,dat,ls='-.',color='black',label="57Fe number density")

# B-field
dat = loadtxt("data/Bfields-R.dat")
print("B:		{} T".format(np.nanmax(dat[:,1])))
dat[:,1] = dat[:,1] / np.nanmax(dat[:,1])
ax2.plot(dat[:,0],dat[:,1],ls='-',color='black',label="Magnetic field strength")

# axes
ax2.set_xlabel("Solar radius")
#ax2.set_ylabel("Plasma frequency [eV]")
#ax2.set_xscale('log')
#ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/solarmodel-black2.jpg')
plt.show()
