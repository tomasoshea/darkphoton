# Tom O'Shea 2023

import numpy as np

# import r & wp
r = np.loadtxt("data/rFrac.dat")
wp = np.loadtxt("data/wp.dat")

# combine
arr = np.array((r, wp)).swapaxes(0,1)

# save full array
np.savetxt("data/rVwp1", arr)

# save array at 0.01 values
red = np.zeros((100,2))
for i in np.arange(0,100):
    for j in range(0,len(r)):
        red[i][0] = round(i/100, ndigits=2)
        if round(r[j]*100, ndigits=2) == i+1:
            red[i][1] = wp[j]
np.savetxt("data/rVwp2", red)
