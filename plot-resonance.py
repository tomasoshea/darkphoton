# plot from .dat files produced by C++ integration

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
#ax2.set( xlim=(0.7,0.8), ylim=(1e-10,1e-2) )
#ax2.set( xlim = (0, 1), ylim=(1e-5, 1.e3) )
x = np.arange(0, 1, 0.001)
#col = 100

m = 1e-2 # eV
log = 1

### CALCULATED ###
dat = loadtxt("RESTstuff/HiddenPhotonFlux_OShea_202309(2).dat", delimiter='\t')    # stripped flux
wp = loadtxt("data/rVwp3.dat")[:,1]    # curve of wp against r
wG = loadtxt("data/wGammaT3.dat")#, usecols=col) # omega Gamma for resonance width
datT = loadtxt("data/flux_m{}-newbins.dat".format(int(np.log10(m))), delimiter='\t')    # pre-calculated flux

totalFlux = np.zeros(len(dat[:,0]))
totalFluxT = np.zeros(len(dat[:,0]))
for i in range(len(dat[0,:]) - 1):
    flux1 =  wG[:,i] * dat[:,i] / ( (m**2 - wp**2)**2 + wG[:,i]**2 )
    flux2 = wG[:,i+1] * dat[:,i+1] / ( (m**2 - wp**2)**2 + wG[:,i+1]**2 )

    totalFlux += -(m**4) * (flux2 - flux1) / 0.1  # integral in 0.1keV steps
    totalFluxT += -(datT[:,i+1] - datT[:,i]) / 0.1

ax2.plot(x, totalFluxT, color='magenta', ls='--', alpha=0.7, label = 'm={}eV, precalculated'.format(m))
ax2.plot(x, totalFlux, color='cyan', alpha=0.7, label = 'm={}eV, calculated'.format(m))

if log:
    name = 'plots/REST_m{}_log.jpg'.format(int(np.log10(m)))    # save name
    ax2.set_yscale('log')
else:
    name = 'plots/REST_m{}.jpg'.format(int(np.log10(m)))    # save name

# axes
ax2.set_xlabel("Solar radius fraction")
ax2.set_ylabel("Flux (cm-2 s-1 keV-1)")
#ax2.set_ylabel("Flux cm-2 s-1 keV-1")
##ax2.set_xscale('log')
#ax2.tick_params(axis='both')
ax2.legend()

fig2.savefig(name)

plt.show()

