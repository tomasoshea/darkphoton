# plot from .dat files produced by C++ integration

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# constants
s2eV = (6.582119569e-16)
m2eV = (1.973269804e-7)


#ax2.set( ylim = (1e0, 1e14) )
#ax2.set( xlim=(0,3e3))#, ylim=(1e-20,1e-8))
"""
atlas = np.loadtxt('data/limits/output_T.dat')

wtab = 10**np.arange(0,4.0002,0.001)
mtab = 10**np.arange(-3,4.1,0.1)

def printdata(im,th=False):
    m = mtab[im]
    w = wtab
    k = w>m
    p=1
    print(m)
    if th:
        p = np.sqrt( (1-m*m/(w*w))*(k) )
    plt.loglog(w,p * atlas[im])


# print the spectrum for the 31st mass
printdata(43)
# with threshold
printdata(43,th=True)
plt.ylabel(r'$\frac{d\Phi}{d\omega} [\rm cm^2 s^{-1} eV^{-1}]$')
plt.xlabel(r'$\omega$ [eV]')
plt.savefig('plots/test.jpg')
plt.show()



"""

# setup plot
log = True				# log or linear normalised
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
if log:	ax2.set( xlim=(1e0, 1e4), ylim=(1e-28, 4e0) )
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
ax2.plot(wtab, p*dat, ls = '--', label = 'T plasmon m = 10 eV')

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
ax2.plot(wtab, p*dat, ls='-', label = 'T plasmon m = 0.1 eV')

# lMixing
dat = loadtxt("data/Espectrum-lMixing-1.dat")
m=0.1		# [eV]
dat[:,1] = dat[:,1]*(m**2)
print("l-Mixing:	max = {}".format(np.nanmax(dat[:,1])))
if log:	dat[:,1] = dat[:,1]/max100
else: dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0], dat[:,1], ls=':', label = 'L mixing')

# pure L
dat = loadtxt("data/limits/output_L.dat")
print("pure L:	max = {}".format(np.nanmax(dat[:,1])))
if log:	dat[:,1] = dat[:,1]/max100
else: dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0], dat[:,1], ls='-.', label = 'pure L')


# axes
ax2.set_xlabel("Dark photon energy [eV]")
ax2.set_ylabel("dPhi/dw per chi2 m2 [cm-2 s-1 eV-3] ")
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









"""

max100 = 7.69581e-08

# T plasmon
dat = loadtxt("data/Espectrum0.dat")	# m = 1 eV
print("m = 1 eV:	max = {}".format(np.nanmax(dat[:,1])/max100))
if log:	dat[:,1] = dat[:,1]/max100
else: dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0], dat[:,1], ls = '-.', label = 'T plasmon m = 1 eV')

dat = loadtxt("data/Espectrum1.dat")	# m = 10 eV
print("m = 10 eV:	max = {}".format(np.nanmax(dat[:,1])/max100))
if log:	dat[:,1] = dat[:,1]/max100
else: dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0], dat[:,1], ls=(0, (3, 5, 3, 5, 1, 5)), label = 'T plasmon m = 10 eV')

dat = loadtxt("data/Espectrum2.dat")	# m = 1e2 eV
print("m = 100 eV:	max = {}".format(np.nanmax(dat[:,1])/max100))
if log:	dat[:,1] = dat[:,1]/max100
else: dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0], dat[:,1], ls=(0, (3, 5, 1, 5, 1, 5)), label = 'T plasmon m = 100 eV')

dat = loadtxt("data/Espectrum3.dat")	# m = 1e3 eV
print("m = 1 keV:	max = {}".format(np.nanmax(dat[:,1])/max100))
if log:	dat[:,1] = dat[:,1]/max100
else: dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0], dat[:,1], ls = '--', label = 'T plasmon m = 1 keV')


# lMixing
dat = loadtxt("data/Espectrum-lMixing-1.dat")
print("l-Mixing:	max = {}".format(np.nanmax(dat[:,1])/max100))
if log:	dat[:,1] = dat[:,1]/max100
else: dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0], dat[:,1], ls=':', label = 'L mixing')


# pureL
dat = loadtxt("data/Espectrum-pureL-1.dat")
print("pure-l:		max = {}".format(np.nanmax(dat[:,1])/max100))
if log:	dat[:,1] = dat[:,1]/max100
else: dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0], dat[:,1], ls='-', label = 'pure L')


# pp-chain
#ax2.vlines(5.49e6, 0, 1, label = 'deuterium fusion')	# deuterium fusion
#ax2.vlines(5.11e5, 0, 1, label = 'e+ e- annihilation', color='red')	# e+ e- annihilation
"""


