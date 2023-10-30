# comparing different contributions for babyIAXO

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
#ax2.set(xlim=(1e-4,1e7), ylim=(1e-35,2e1))
#ax2.set(xlim=(1e-3,1e4), ylim=(1e-5,2e0))

s2eV = (6.582119569e-16)
m2eV = (1.973269804e-7)
sigma = 1e3		# eV
E0 = 559.1987036540328e3	# eV

def Lorentz(E1,theta):
	gamma = 1.046099691693158
	beta = 0.2935886968984385
	return gamma*E1*(1 + beta*np.cos(theta))	# eV

def Gauss(E):
	return np.exp(-0.5*((E-E0)/sigma)**2)

ECoM = np.linspace(E0-5*sigma,E0+5*sigma,1000)
thetaCoM = np.linspace(0,2*np.pi,1000)
Elab = []

for E in ECoM:
	for th in thetaCoM:
		Elab.append(Lorentz(E,th))

#print(Elab)
plt.hist(Elab,bins=1000)

# axes
ax2.set_xlabel("Lab frame photon energy [eV]")
ax2.set_ylabel("Flux")
#ax2.set_xscale('log')
#ax2.set_yscale('log')
#ax2.legend()

plt.savefig('plots/annihilation1.jpg')
plt.show()
