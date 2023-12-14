# comparing different contributions for babyIAXO

import numpy as np
from numpy import loadtxt
from numpy import sqrt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(1e-3,7e3), ylim=(1e-15,1e-7))
y2= 1e-5
alp=1

# solar bounds
col='lime'
z=0.6
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Solar.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,label='Stellar',facecolor=col,edgecolor=None,zorder=z)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/RG.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,edgecolor=None,zorder=z)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/HB.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,edgecolor=None,zorder=z)


# XENON
col='hotpink'
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

col='black'
lw=5
dat = loadtxt("data/limits/stats-babyIAXO-Atlas-1eV.dat")
ax2.plot(dat[:,0],dat[:,1], label='BabyIAXO',color=col,ls='--',lw=lw)

dat = loadtxt("data/limits/stats-baselineIAXO-Atlas-1eV.dat")
ax2.plot(dat[:,0],dat[:,1], label='IAXO',color=col,ls='-',lw=lw)

dat = loadtxt("data/limits/stats-upgradedIAXO-Atlas-1eV.dat")
ax2.plot(dat[:,0],dat[:,1], label='IAXO+',color=col,ls=':',lw=lw)


lw=2
# gas Atlas
dat = loadtxt("data/limits/stats-babyIAXO-AtlasGas-30eV.dat")
ax2.plot(dat[:,0],dat[:,1], color=col,ls='--',lw=lw)

dat = loadtxt("data/limits/stats-baselineIAXO-AtlasGas-30eV.dat")
ax2.plot(dat[:,0],dat[:,1], color=col,ls='-',lw=lw)

dat = loadtxt("data/limits/stats-upgradedIAXO-AtlasGas-30eV.dat")
ax2.plot(dat[:,0],dat[:,1], color=col,ls=':',lw=lw)

# axes
ax2.set_xlabel("Dark photon mass [eV]")
ax2.set_ylabel("Kintetic Mixing Parameter")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()#loc='upper right')

plt.savefig('plots/limit-comparison-2.jpg')
plt.show()

