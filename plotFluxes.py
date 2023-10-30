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
#ax2.set(xlim=(1e-4,1e7), ylim=(1e-35,2e1))
ax2.set(xlim=(1e-3,1e4), ylim=(1e-5,2e0))


# tPlasmon
dat = loadtxt("data/limits/babyIAXO-tPlasmon-flux.dat")
top = np.nanmax(dat[:,1])
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1]**(1/4),lw=5,color='black',label='T-plasmon')

#dat = loadtxt("data/limits/babyIAXO-tPlasmon-70eV.dat")
#dat[:,1] = dat[:,1] / top
#ax2.plot(dat[:,0],dat[:,1], ls=':', label='T-plasmon (70 eV)')

# tPlasmon gas
dat = loadtxt("data/limits/babyIAXO-tPlasmon-100eV-gas.dat")
dat[:,1] = dat[:,1] / top
ax2.plot(dat[:,0],dat[:,1]**(1/4),lw=2,color='black',label='T-plasmon (gas)')

# tPlasmon gas T=3K
#dat = loadtxt("data/limits/babyIAXO-tPlasmon-100eV-3K-gas.dat")
#dat[:,1] = dat[:,1] / top
#ax2.plot(dat[:,0],dat[:,1], ls=":", label='T-plasmon (cold gas)')

# lMixing
#dat = loadtxt("data/limits/babyIAXO-1-lMixingRes-flux.dat")
#dat[:,1] = dat[:,1] / top
#ax2.plot(dat[:,0],dat[:,1],lw=5,ls=':',color='black',label='L-plasmon mixing')
#
## lMixing gas
#dat = loadtxt("data/limits/babyIAXO-1-lMixingResGas-flux.dat")
#dat[:,1] = dat[:,1] / top
#ax2.plot(dat[:,0],dat[:,1],lw=2,ls=':',color='black',label='L-plasmon mixing (gas)')
#
## pureL idealised
#dat = loadtxt("data/limits/babyIAXO-perfect-pureL.dat")
#dat[:,1] = dat[:,1] / top
##ax2.plot(dat[:,0],dat[:,1], label='L-DP conversion (idealised)')
#
## pureL B-field
#dat = loadtxt("data/limits/babyIAXO-new2-pureL.dat")
#dat[:,1] = dat[:,1] / top
#ax2.plot(dat[:,0],dat[:,1],ls='-.',color='black',label='L-DP conversion')
#
## pp chain
#dat = loadtxt("data/limits/babyIAXO-pp-new.dat")
#dat[:,1] = dat[:,1] / top
#ax2.plot(dat[:,0],dat[:,1],lw=5,ls='--',color='black',label='pp-chain')
#
## ee annihilation
#dat = loadtxt("data/limits/babyIAXO-pp-ee-new.dat")
#dat[:,1] = dat[:,1] / top
#ax2.plot(dat[:,0],dat[:,1],lw=3,ls='--',color='black',label='e+ e- annihilation')
#
## 57Fe
#dat = loadtxt("data/limits/babyIAXO-57Fe-3.dat")
#dat[:,1] = dat[:,1] / top
#ax2.plot(dat[:,0],dat[:,1],lw=2,ls='--',color='black',label='57Fe')

# 55Mn
#dat = loadtxt("data/limits/babyIAXO-55Mn-1.dat")
#dat[:,1] = dat[:,1] / top
#ax2.plot(dat[:,0],dat[:,1], label='55Mn')

# Javier Atlas data
atlas = loadtxt("data/limits/output_T.dat")
wtab = 10**np.arange(0,4.0002,0.001)
mtab = 10**np.arange(-3,4.1,0.1)
def boundIAXO(Am2=2.3,ty=1,bk=1e-8,acm2=8*0.15,Lm=20,LowEthres=1e3):
	# for result
	bI = mtab*0
	counts = mtab*0
	# wbin widths
	dw = wtab[1:]-wtab[:-1]
	# wbin centers
	wm = (wtab[1:]+wtab[:-1])/2
	# area in m2 into cm2
	Acm2 = Am2*1e4
	# time in seconds
	ts   = ty*np.pi*1e7
	# background counts per eV, flat background
	nb = bk/1000*acm2*ts
	print('Nb = %f counts/keV'%(nb*1000))
	# integrate
	for im in range(len(mtab)):
		x  = mtab[im]**2*Lm/(1.97e-7)/(4*wm)
		x2 = x*x
		pp = 4*x2/(1+2*x2)
		f0 = atlas[im]
		f1 = (f0[1:]+f0[:-1])/2
		th = (wm>LowEthres)*(wm>mtab[im])*(1-mtab[im]**2/wm**2)
		
		# counts, missing chi^4
		ng = mtab[im]**4*(f1*dw*np.sqrt(th)*pp*Acm2*ts).sum()
		counts[im] = ng
		# chi limit from bg=nb
		bI[im] = 1/(ng/nb)**(1/4)
	return counts * ( (m2eV*m2eV / 2.3 ) * ( s2eV / (np.pi*1e7) ) ) 	# counts => eV3

dat = boundIAXO(LowEthres=1e0)
np.savetxt('data/atlas_mtab.dat', mtab, delimiter=',')
np.savetxt('data/atlas_1e0.dat', dat, delimiter=',')
dat = dat / top
ax2.plot(mtab,dat**(1/4),lw=2,ls='-',color='magenta',label='Atlas (1eV)')
print(len(mtab))
print(len(dat))

dat = boundIAXO(LowEthres=1e2)
dat = dat / top
ax2.plot(mtab,dat**(1/4),lw=2,ls='-',color='magenta',label='Atlas (100eV)')
np.savetxt('data/atlas_1e2.dat', dat, delimiter=',')


# prob
#fig = plt.figure(2)	# display is 1920 x 1080 (16:9)
#ax = fig.add_axes((.1,.1,.8,.8))
#ax.set(xlim=(1e-3,1e7), ylim=(0.,1.1))
#dat = loadtxt("data/highEprob.dat")
#ax.plot(dat[:,0],dat[:,1], label='conversion prob in vacuum')
#ax.set_xscale('log')

# axes
ax2.set_xlabel("Dark photon mass [eV]")
ax2.set_ylabel("Flux^(1/4) (normalised)")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/flux-comparison-Atlas3.jpg')
plt.show()
