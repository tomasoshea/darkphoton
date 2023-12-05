# comparing different contributions for babyIAXO

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

s2eV = (6.582119569e-16)
m2eV = (1.973269804e-7)

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
#ax2.set(xlim=(1e-4,1e7), ylim=(1e-35,2e0))
ax2.set(xlim=(1e-3,1e3), ylim=(1e-27,2e0))
lw = 4



############## Atlas data #################
fs=30

# T-plasmon Altas (vacuum)
dat = loadtxt("data/limits/babyIAXO-Atlas-1eV.dat")
top = np.nanmax(dat[:,1])
print(top/(((100*m2eV)**2)*s2eV))		# cm-2 s-1
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='1 eV', color='black', ls='-', lw=4)

dat = loadtxt("data/limits/babyIAXO-Atlas-10eV.dat")
top = np.nanmax(dat[:,1])
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='10 eV', color='black', ls='--', lw=4)

dat = loadtxt("data/limits/babyIAXO-Atlas-100eV.dat")
top = np.nanmax(dat[:,1])
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='100 eV', color='black', ls=':', lw=4)

dat = loadtxt("data/limits/babyIAXO-Atlas-1keV.dat")
top = np.nanmax(dat[:,1])
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='1 keV', color='black', ls='-.', lw=4)


# axes
ax2.set_xlabel("Dark photon mass [eV]")
ax2.set_ylabel("Flux (normalised)")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()#loc='upper right')

plt.savefig('plots/flux-comparison-E2.jpg')
plt.show()




"""




# tPlasmon
dat = loadtxt("data/limits/babyIAXO-tPlasmon-flux.dat")
top = np.nanmax(dat[:,1])
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='T-plasmon', color='black', ls='-')


# tPlasmon gas
dat = loadtxt("data/limits/babyIAXO-tPlasmon-1e-6-gas.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='T-plasmon (gas)', color='black',ls=(5,(10,3)))


# lMixing
dat = loadtxt("data/limits/babyIAXO-1-lMixingRes-flux.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='L-plasmon mixing',ls=(0,(5,6)))

# lMixing gas
dat = loadtxt("data/limits/babyIAXO-1-lMixingResGas-flux.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='L-plasmon mixing (gas)',ls=(0,(5,2)))

# pureL B-field
dat = loadtxt("data/limits/babyIAXO-new2-pureL.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='L-DP conversion',ls=(0,(1,3)))

# pp chain
dat = loadtxt("data/limits/babyIAXO-pp-new.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='pp-chain',ls=':')

# ee annihilation
dat = loadtxt("data/limits/babyIAXO-pp-ee-new.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='e+ e- annihilation',ls='-.')

# 57Fe
dat = loadtxt("data/limits/babyIAXO-57Fe-3.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='57Fe',ls=(0, (3, 5, 1, 5, 1, 5)))




# T-plasmon old
dat = loadtxt("data/limits/babyIAXO-tPlasmon-old.dat")
dat[:,1] = dat[:,1] / top
#ax2.plot(dat[:,0],dat[:,1], label='T-plasmon old', color='red', ls='-', lw=5)

# T-plasmon old (gas)
dat = loadtxt("data/limits/babyIAXO-tPlasmon-old-100eV-gas.dat")
dat[:,1] = dat[:,1] / top
##ax2.plot(dat[:,0],dat[:,1], label='T-plasmon old (gas)', color='red',ls='-', lw=2)


dat = loadtxt("data/limits/babyIAXO-Atlas-test-gas.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='T-plasmon high p (gas)', color='red',ls='--', lw=2)

# T-plasmon old old (gas)
dat = loadtxt("data/limits/babyIAXO-tPlasmon-1e-6-gas.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='T-plasmon old-old (gas)', color='magenta', lw=2)

#dat = loadtxt("data/limits/babyIAXO-Atlas-100eV.dat")
#dat[:,1] = dat[:,1] / top
#ax2.plot(dat[:,0],dat[:,1], label='T-plasmon', color='red', ls='-', lw=5)


"""
