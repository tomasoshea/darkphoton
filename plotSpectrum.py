# plot from .dat files produced by C++ integration

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# constants
s2eV = (6.582119569e-16)
m2eV = (1.973269804e-7)

# setup plot
log = True				# log or linear normalised
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
if log:	ax2.set( xlim=(1e0, 1e4), ylim=(1e-12, 1.1e0) )
else: ax2.set( xlim=(1e-2, 3e4), ylim=(0, 1.01))


# get data and plot
tPlasmon = loadtxt("data/limits/output_T.dat")
wtab = 10**np.arange(0,4.0002,0.001)
mtab = 10**np.arange(-3,4.1,0.1)


# T-plasmon 100 eV
i=40
m=mtab[i]
dat = tPlasmon[i]*(m**2)
#print(m)
k = wtab>m
p = np.sqrt( (1-m*m/(wtab*wtab))*(k) )
print("m = 10 eV:	max = {}".format(np.nanmax(dat)))
max100 = np.nanmax(dat)
if log:	dat = dat/max100
else: dat = dat/np.nanmax(dat)
ax2.plot(wtab, p*dat, ls = '--', label = 'T-DP m = 10 eV')

# T-plasmon 0.1 eV
i=20
m=mtab[i]
dat = tPlasmon[i]*(m**2)
#print(m)
k = wtab>m
p = np.sqrt( (1-m*m/(wtab*wtab))*(k) )
print("m = 0.1 eV:	max = {}".format(np.nanmax(dat)))
if log:	dat = dat/max100
else: dat = dat/np.nanmax(dat)
ax2.plot(wtab, p*dat, ls='-', label = 'T-DP m = 0.1 eV')

# pure L
dat = loadtxt("data/limits/output_L.dat")
print("pure L:	max = {}".format(np.nanmax(dat[:,1])))
if log:	dat[:,1] = dat[:,1]/max100
else: dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0], dat[:,1], ls=':', label = 'L-DP')


# axes
ax2.set_xlabel("Dark photon energy [eV]")
ax2.set_ylabel(r'$\frac{1}{m^2} \frac{d \Phi}{d \omega}$ (normalised)')
#ax2.set_ylabel("Differential flux [eV2]")
ax2.set_xscale('log')
if log:
	ax2.set_yscale('log')
	ext = "-log.jpg"
else:
	ax2.set_yscale('linear')
	ext = "-normalised.jpg"
ax2.tick_params(axis='both')
ax2.legend()

fig2.savefig('plots/spectrumAtlas{}'.format(ext))

plt.show()
