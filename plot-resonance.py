# plot from .dat files produced by C++ integration

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set( xlim=(0, 1), ylim=(1e1, 1e20) )
#ax2.set( xlim = (0, 1), ylim=(-0.01, 1.01) )
x = np.arange(0, 1, 0.01)
col = 100


### CALCULATED ###
sup = loadtxt("data/flux_0.dat", usecols=col)	# suppressed flux * m-4 * chi-2
res = loadtxt("data/flux_1.dat", usecols=col)    # resonant flux * m-4 * chi-2
unsup = loadtxt("data/flux_2.dat", usecols=col)	# unsuppressed flux * chi-2
wp = loadtxt("data/rvwp2")[:,1]    # curve of wp against r
wG = loadtxt("data/wGammaT2.dat", usecols=col) # omega Gamma for resonance width
chi = 1e-11

m = 1e2 # eV
dat1 = []
buff = 2

"""for j in range(buff):
    print(j)
    if m < wp[j]:
        item = (m**4) * (chi**2) * sup[j]
        dat1.append(item)
    elif m > wp[j]:
        item = (chi**2) * unsup[j]
        dat1.append(item)
    else:
        item = (m**4) * (chi**2) * res[j]
        dat1.append(item)
"""
for i in range( len(sup) ):

    if ( wG[i]**2 < ( m**2 - wp[i]**2 )**2 ):
        if m < wp[i]:
            item = (m**4) * (chi**2) * sup[i]
            dat1.append(item)
        elif m > wp[i]:
            item = (chi**2) * unsup[i]
            dat1.append(item)
    else:
        item = (m**4) * (chi**2) * res[i]
        dat1.append(item)

"""for i in range(buff):
    j = len(sup) - buff + i
    print(j)
    if m < wp[j]:
        item = (m**4) * (chi**2) * sup[j]
        dat1.append(item)
    elif m > wp[j]:
        item = (chi**2) * unsup[j]
        dat1.append(item)
    else:
        item = (m**4) * (chi**2) * res[j]
        dat1.append(item)
"""
#dat1 = dat1 / np.nanmax(dat1)
ax2.plot(x, dat1, color='magenta', ls='--')

"""
m = 1e2 # eV
dat1 = []
for i in range(len(sup)-1):
    if m < wp[i]:
        item = (m**4) * (chi**2) * sup[i]
        dat1.append(item)
    elif m > wp[i] and m < wp[i+1]:
        item = (m**4) * (chi**2) * res[i]
        dat1.append(item)
    elif m > wp[i]:
        item = (chi**2) * unsup[i]
        dat1.append(item)

dat1.append(0)
dat1 = dat1 / np.nanmax(dat1)
ax2.plot(x, dat1, color='red', ls='--')

m = 1e1 # eV
dat1 = []
for i in range(len(sup)-1):
    if m < wp[i]:
        item = (m**4) * (chi**2) * sup[i]
        dat1.append(item)
    elif m > wp[i] and m < wp[i+1]:
        item = (m**4) * (chi**2) * res[i]
        dat1.append(item)
    elif m > wp[i]:
        item = (chi**2) * unsup[i]
        dat1.append(item)

dat1.append(0)
dat1 = dat1 / np.nanmax(dat1)
ax2.plot(x, dat1, color='cyan', ls='--')

m = 1e-1 # eV
dat1 = []
for i in range(len(sup)-1):
    if m < wp[i]:
        item = (m**4) * (chi**2) * sup[i]
        dat1.append(item)
    elif m > wp[i] and m < wp[i+1]:
        item = (m**4) * (chi**2) * res[i]
        dat1.append(item)
    elif m > wp[i]:
        item = (chi**2) * unsup[i]
        dat1.append(item)

dat1.append(0)
#dat1 = dat1 / np.nanmax(dat1)
ax2.plot(x, dat1, color='blue', ls='--')"""



### PRESET ###

#datT = loadtxt("data/flux_m3_X-11.dat", usecols=col)	# m = 1 keV
#datT = datT / np.nanmax(datT)
#ax2.plot(x, datT, color='magenta', ls=':', label = 'm = 1 keV')

datT = loadtxt("data/flux_m2_X-11.dat", usecols=col)	# m = 100 eV
#datT = datT / np.nanmax(datT)
ax2.plot(x, datT, color='red', label = 'm = 100 eV')

datT = loadtxt("data/flux_m2_X-11-again.dat", usecols=col)	# m = 100 eV
#datT = datT / np.nanmax(datT)
ax2.plot(x, datT, color='green', ls=':')

#
#datT = loadtxt("data/flux_m1_X-11.dat", usecols=col)	# m = 10 eV
#datT = datT / np.nanmax(datT)
#ax2.plot(x, datT, color='cyan', ls=':', label = 'm = 10 eV')

#datT = loadtxt("data/flux_m0_X-11.dat", usecols=col)	# m = 1 eV
#datT = datT / np.nanmax(datT)
#ax2.plot(x, datT, color='blue', label = 'm = 1 eV, w = 500 eV')

#datT = loadtxt("data/flux_m-1_X-11.dat", usecols=col)	# m = 1 eV
#datT = datT / np.nanmax(datT)
#ax2.plot(x, datT, color='blue', ls=':', label = 'm = 0.1 eV')

#datT = loadtxt("data/flux_m3_X-11.dat", usecols=col)	# m = 1 eV
#datT = datT / np.nanmax(datT)
#ax2.plot(x, datT, color='blue', ls=':', label = 'm = 1 keV')

# axes
ax2.set_xlabel("Solar radius fraction")
ax2.set_ylabel("Flux (cm-2 s-1 keV-1)")
#ax2.set_ylabel("Flux cm-2 s-1 keV-1")
##ax2.set_xscale('log')
ax2.set_yscale('log')
#ax2.tick_params(axis='both')
ax2.legend()

fig2.savefig('plots/fluxcompar2norm.jpg')

plt.show()

