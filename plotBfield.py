# plot solar B-field model from .dat files produced in c++ integration

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig = plt.figure(2, figsize = [16, 9], dpi = 120)
ax = fig.add_axes((.1,.1,.8,.8))
ax.set( xlim=(0,1), ylim=(-5., 3.01e3) )
ax2 = ax.twiny()

# import B-field
dat = loadtxt("data/Bfields.dat")
datR = loadtxt("data/Bfields-R.dat")
wp = dat[:,0]
r = datR[:,0]

# plot B-field
ax2.plot(dat[::,0], dat[:,1], alpha=0)
ax.plot(datR[1:-2,0], datR[1:-2,1], color="magenta")

ax2.set_xlabel("plasma frequency [eV]")
ax.set_xlabel("solar radius fraction")
ax.set_ylabel("B-field [T]")

# sort out secondary axis and ticks
xticks = [0, 43.1, 85.8, 121.7, 174.1, 220.7, 260.1, 289.5]
xticklabels = [290, 200, 100, 50, 20, 10, 5, 1]
ax2.set_xticks(xticks)
ax2.set_xticklabels(xticklabels)

fig.savefig('plots/Bfields.jpg')

plt.show()
