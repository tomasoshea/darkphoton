# comparing different contributions for babyIAXO

# conversion factors
s2eV = (6.582119569e-16)
J2eV = (1. / 1.602176634e-19)
m2eV = (1.973269804e-7)
K2eV = (8.617333262e-5)
kg2eV = 5.609588604e35

import numpy as np
from numpy import loadtxt
from numpy import sqrt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(1e-3,1e4), ylim=(1e-15,1e-7))

# solar bounds
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Solar.txt")
ax2.plot(dat[:,0],dat[:,1], label='Stellar', color='black',ls='-', lw=4)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/RG.txt")
ax2.plot(dat[:,0],dat[:,1], color='black',ls='-', lw=4)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/HB.txt")
ax2.plot(dat[:,0],dat[:,1], color='black',ls='-', lw=4)


# XENON
col='cyan'
y2= 1e-5
alp=1
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Xenon1T.txt")
dat[:,1] = dat[:,1]*sqrt(0.3/0.45)
plt.fill_between(1e3*dat[:,0],dat[:,1],y2=y2,label='XENON',edgecolor=None,facecolor=col,zorder=0.5,alpha=alp)

dat = loadtxt("DPlimits/limit_data/DarkPhoton/Xenon1T_S1S2.txt")
dat[:,1] = dat[:,1]*sqrt(0.3/0.45)
plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.5,alpha=alp)

dat = loadtxt("DPlimits/limit_data/DarkPhoton/XENON1T_SE.txt")
dat[:,1] = dat[:,1]*sqrt(0.3/0.45)
plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.5,alpha=alp)

dat = loadtxt("DPlimits/limit_data/DarkPhoton/XENON1T_Solar_SE.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.5,alpha=alp)

dat = loadtxt("DPlimits/limit_data/DarkPhoton/XENONnT.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.5,alpha=alp)



# vacuum stats
dat = loadtxt("data/limits/statsAtlas-babyIAXO-1eV-test.dat")
ax2.plot(dat[:,0],dat[:,1], label='BabyIAXO (stats)', color='red',ls='--', lw=4)

dat = loadtxt("data/limits/statsAtlas-baselineIAXO-1eV-test.dat")
ax2.plot(dat[:,0],dat[:,1], label='IAXO (stats)', color='red',ls='-', lw=4)

dat = loadtxt("data/limits/statsAtlas-upgradedIAXO-1eV-test.dat")
ax2.plot(dat[:,0],dat[:,1], label='IAXO+ (stats)', color='red',ls=':', lw=4)


# vacuum naive
dat = loadtxt("data/limits/babyIAXO-Atlas-1eV.dat")
a = 0.6 * 1e-4 		# m2
A = 0.77			# m2
bg = 1e-7	# cm-2 s-1 keV-1
bg *= 10 * (100*m2eV)**2 * s2eV	# eV3
dat[:,1] = ((bg/dat[:,1]) * a/A) **(1/4)
ax2.plot(dat[:,0],dat[:,1], label='BabyIAXO (bg = sig)', color='green',ls='--', lw=4)

dat = loadtxt("data/limits/baselineIAXO-Atlas-1eV.dat")
a = 1.2 * 1e-4 		# m2
A = 2.3				# m2
bg = 1e-8	# cm-2 s-1 keV-1
bg *= 10 * (100*m2eV)**2 * s2eV	# eV3
dat[:,1] = ((bg/dat[:,1]) * a/A) **(1/4)
ax2.plot(dat[:,0],dat[:,1], label='IAXO (bg = sig)', color='green',ls='-', lw=4)

dat = loadtxt("data/limits/upgradedIAXO-Atlas-1eV.dat")
a = 1.2 * 1e-4 		# m2
A = 3.9				# m2
bg = 1e-9	# cm-2 s-1 keV-1
bg *= 10 * (100*m2eV)**2 * s2eV	# eV3
dat[:,1] = ((bg/dat[:,1]) * a/A) **(1/4)
ax2.plot(dat[:,0],dat[:,1], label='IAXO+ (bg = sig)', color='green',ls=':', lw=4)


# gas Atlas
#dat = loadtxt("data/limits/stats-babyIAXO-AtlasGas-30eV.dat")
#ax2.plot(dat[:,0],dat[:,1], label='BabyIAXO (gas)', color='g',ls='--', lw=3)
#
#dat = loadtxt("data/limits/stats-baselineIAXO-AtlasGas-30eV.dat")
#ax2.plot(dat[:,0],dat[:,1], label='IAXO (gas)', color='g',ls='-', lw=3)
#
#dat = loadtxt("data/limits/stats-upgradedIAXO-AtlasGas-30eV.dat")
#ax2.plot(dat[:,0],dat[:,1], label='IAXO+ (gas)', color='g',ls=':', lw=3)



# axes
#ax2.set_xlabel("Dark photon mass [eV]")
#ax2.set_ylabel("Flux (normalised)")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()#loc='upper right')

plt.savefig('plots/limit-comparison-test.jpg')
plt.show()


## gas Atlas
#dat = loadtxt("data/limits/statsAtlas-babyIAXO-30eV-gas.dat")
#ax2.plot(dat[:,0],dat[:,1], label='BabyIAXO (gas)', color='g',ls='--', lw=3)
#
#dat = loadtxt("data/limits/statsAtlas-baselineIAXO-30eV-gas.dat")
#ax2.plot(dat[:,0],dat[:,1], label='IAXO (gas)', color='g',ls='-', lw=3)
#
#dat = loadtxt("data/limits/statsAtlas-upgradedIAXO-30eV-gas.dat")
#ax2.plot(dat[:,0],dat[:,1], label='IAXO+ (gas)', color='g',ls=':', lw=3)
