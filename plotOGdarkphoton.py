# plot from .dat files produced by C++ integration

from matplotlib import pyplot as plt
from PlotFuncs import *

plt.style.use("style.txt")	# import plot style

# setup plot
fig = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax = fig.add_axes((.1,.1,.8,.8))
ax.set(xlim=(3e-18,7e5), ylim=(1.0e-18,1e0))

# DPDM
DarkPhoton.DarkMatter(ax,text_on=0)

# Axion haloscopes
DarkPhoton.Haloscopes(ax)

# # LSW/Helioscopes
DarkPhoton.LSW(ax)
DarkPhoton.CAST(ax)
DarkPhoton.SHIPS(ax)

# Tests of coulomb law
DarkPhoton.Coulomb(ax)

# # Reactor neutrinos
DarkPhoton.TEXONO(ax)

# # Geomagnetic field
DarkPhoton.SuperMAG(ax)

# # DPDM searches
DarkPhoton.Xenon(ax)
DarkPhoton.DAMIC(ax)
DarkPhoton.SENSEI(ax)
DarkPhoton.SuperCDMS(ax)
DarkPhoton.FUNK(ax)
DarkPhoton.LAMPOST(ax)
DarkPhoton.Tokyo(ax)
DarkPhoton.SHUKET(ax)
DarkPhoton.DarkEfield(ax)
DarkPhoton.WISPDMX(ax)
DarkPhoton.SQuAD(ax)
DarkPhoton.DMPathfinder(ax)
DarkPhoton.ORPHEUS(ax)
DarkPhoton.MuDHI(ax)
DarkPhoton.DOSUE(ax)
DarkPhoton.FAST(ax)

# # Astrophysical boundse
DarkPhoton.StellarBounds(ax)
DarkPhoton.COBEFIRAS(ax)
DarkPhoton.Jupiter(ax)
DarkPhoton.Earth(ax)
DarkPhoton.Crab(ax)
DarkPhoton.IGM(ax)
DarkPhoton.LeoT(ax)
DarkPhoton.GasClouds(ax)
DarkPhoton.NeutronStarCooling(ax)


# ADDING IAXO BOUNDS!!
#DarkPhoton.IAXO(ax, pureL=True, text_on=False)

# BHSR
plt.fill_between([6.5e-15,2.9e-11],[1e-18,1e-18],y2=1,color='gray',edgecolor='none',zorder=-100,alpha=0.25)
plt.gcf().text(0.304,0.176,r'Black hole',fontsize=23,ha='center',rotation=0,color='k',alpha=0.8,path_effects=line_background(2,'w'))
plt.gcf().text(0.304,0.145,r'superradiance',fontsize=23,ha='center',rotation=0,color='k',alpha=0.8,path_effects=line_background(2,'w'))

# Final label
plt.arrow(0.435, 0.375, 0, -0.055, transform=fig.transFigure,figure=fig,
  length_includes_head=True,lw=2.5,
  head_width=0.012, head_length=0.028, overhang=0.13,
  edgecolor='k',facecolor='w',clip_on=False,zorder=-1)

plt.text(4e-9,0.8e-14,r'Dark',fontsize=27,ha='center')
plt.text(4e-9,0.15e-14,r'photon',fontsize=27,ha='center')
plt.text(4e-9,0.02e-14,r'DM',fontsize=27,ha='center')

# axes
ax.set_xlabel("Dark photon mass [eV]")
ax.set_ylabel("Kinetic Mixing Parameter")
ax.set_xscale('log')
ax.set_yscale('log')
ax.tick_params(axis='both')
#ax.legend()

fig.savefig('plots/darkphoton2.png')

plt.show()
