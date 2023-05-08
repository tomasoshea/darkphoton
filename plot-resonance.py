# plot from .dat files produced by C++ integration

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
#ax2.set( xlim=(1e0, 1e3), ylim=(-0.01, 1.01) )
#ax2.set( ylim = (1e0, 1e14) )


# get data and plot

# T plasmon
datT = loadtxt("data/flux_m2_X-11.dat", usecols=4)	# m = 100 eV
x = np.arange(0, 1, 0.01)
ax2.plot(x, datT, color='magenta', label = 'm = 100 eV, w = 500 eV')


# axes
ax2.set_xlabel("Solar radius fraction")
#ax2.set_ylabel("Flux (normalised)")
ax2.set_ylabel("Flux cm-2 s-1 keV-1")
##ax2.set_xscale('log')
##ax2.set_yscale('log')
#ax2.tick_params(axis='both')
ax2.legend()

fig2.savefig('plots/fluxplot1.jpg')

plt.show()

