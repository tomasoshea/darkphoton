# comparing different contributions for babyIAXO

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
#ax2.set(xlim=(1e-4,1e7), ylim=(1e-35,2e0))
ax2.set(xlim=(1e-3,1e7), ylim=(1e-35,2e1))



############## Atlas data #################

# T-plasmon Altas (vacuum)
dat = loadtxt("data/limits/babyIAXO-Atlas-1eV.dat")
top = np.nanmax(dat[:,1])
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='T-plasmon', color='black', ls='-', lw=5)

# T-plasmon Atlas (gas)
dat = loadtxt("data/limits/babyIAXO-Atlas-30eV-gas.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='T-plasmon (gas)', color='black',ls='-', lw=2)

# L-plasmon mixing (vacuum)
dat = loadtxt("data/limits/babyIAXO-lMixing-x.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='L-plasmon', color='black', ls='--', lw=5)

# L-plasmon mixing (gas)
dat = loadtxt("data/limits/babyIAXO-lMixing-x-gas.dat")
dat[:,1] = dat[:,1] / top
for i in range(len(dat[:,0])):
		if dat[i,0] >= 0.835: dat[i,1] = 0
ax2.plot(dat[:,0],dat[:,1], label='L-plasmon (gas)', color='black',ls='--', lw=2)

# pureL B-field
dat = loadtxt("data/limits/babyIAXO-Atlas-pureL-x.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='L-DP conversion',ls=(0,(10,3)),color='black')

# pureL idealised conversion
#dat = loadtxt("data/limits/babyIAXO-Atlas-pureL-100eV-ideal.dat")
#dat[:,1] = dat[:,1] / top
#ax2.plot(dat[:,0],dat[:,1], label='L-DP conversion (idealised)',ls=(0,(1,3)))

# pp chain
dat = loadtxt("data/limits/babyIAXO-pp-new.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='pp-chain',ls=':',color='black', lw=5)

# ee annihilation
dat = loadtxt("data/limits/babyIAXO-pp-ee-new.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1], label='e+ e- annihilation',ls='-.',color='black', lw=5)

# 57Fe
dat = loadtxt("data/limits/babyIAXO-57Fe-3.dat")
dat[:,1] = dat[:,1] / top
dat[-1,1] = 0
ax2.plot(dat[:,0],dat[:,1], label='57Fe', color='black', ls=(0, (3, 5, 1, 5, 1, 5)), lw=5)

# 55Mn
#dat = loadtxt("data/limits/babyIAXO-55Mn-1.dat")
#dat[:,1] = dat[:,1] / top
#ax2.plot(dat[:,0],dat[:,1], label='55Mn')


# axes
ax2.set_xlabel("Dark photon mass [eV]")
ax2.set_ylabel("Flux (normalised)")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend(loc='upper right')

plt.savefig('plots/flux-comparison-x.jpg')
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
