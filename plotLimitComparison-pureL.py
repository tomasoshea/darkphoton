# comparing different contributions for babyIAXO

import numpy as np
from numpy import loadtxt
from numpy import sqrt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(1e-6,1e3), ylim=(1e-15,1e-7))
y2= 1e-5
alp=1

# solar bounds
col='deepskyblue'
z=0.6
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Solar.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,label='Stellar',facecolor=col,edgecolor=None,zorder=z)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/RG.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,edgecolor=None,zorder=z)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/HB.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,edgecolor=None,zorder=z)


# XENON
col='violet'
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

# DP DM
col='blue'
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Cosmology_Arias.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,label='Dark matter',edgecolor=None,facecolor=col,zorder=0.5,alpha=alp)

# haloscopes
fs=17
projection=True
text_on=True
col='darkred'
y2 = ax2.get_ylim()[1]
zo = 0.3
# ADMX
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/ADMX.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1,lw=3,label='Haloscopes')
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/ADMX2018.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/ADMX2019_1.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/ADMX2019_2.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/ADMX2021.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/ADMX_Sidecar.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
# HAYSTAC
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/HAYSTAC.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/HAYSTAC_2020.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/HAYSTAC_2022.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
# CAPP
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/CAPP-1.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/CAPP-2.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/CAPP-3.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/CAPP-4.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/CAPP-5.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/CAPP-6.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/CAST-CAPP.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/TASEH.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/ORGAN-1a.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.21)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/QUAX.txt")
dat[0,1] = 1e0
plt.plot(dat[:,0],dat[:,1],zorder=0.2,color=col,lw=2)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/QUAX2.txt")
dat[0,1] = 1e0
plt.plot(dat[:,0],dat[:,1],zorder=0.2,color=col,lw=2)
dat = loadtxt("DPlimits/limit_data/DarkPhoton/Rescaled/QUAX3.txt")
dat[0,1] = 1e0
plt.plot(dat[:,0],dat[:,1],zorder=0.2,color=col,lw=2)

# IAXO
col='black'
lw=5
dat = loadtxt("data/limits/stats-babyIAXO-pureL-30eV-ideal-2.dat")
ax2.plot(dat[:,0],dat[:,1], label='L-DP',color=col,ls='--',lw=lw)


# axes
ax2.set_xlabel("Dark photon mass [eV]")
ax2.set_ylabel("Kintetic mixing parameter")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend(loc='upper right')

plt.savefig('plots/limit-comparison-pureL.jpg')
plt.show()

