# plotting of dark photon limits using Ciaran O'Hares plotting script
# https://github.com/cajohare/AxionLimits
# reduced version zoomed in on IAXO region

from PlotFuncs import *

fig,ax = DarkPhoton.FigSetup(reduced=True)

# # DPDM searches
DarkPhoton.Xenon(ax)
DarkPhoton.DAMIC(ax)
DarkPhoton.SENSEI(ax)
DarkPhoton.SuperCDMS(ax)
DarkPhoton.SHUKET(ax)
DarkPhoton.DarkEfield(ax)
DarkPhoton.WISPDMX(ax)
DarkPhoton.SQuAD(ax)
DarkPhoton.DMPathfinder(ax)
DarkPhoton.ORPHEUS(ax)
DarkPhoton.DOSUE(ax)
DarkPhoton.FAST(ax)

# # Astrophysical boundse
DarkPhoton.StellarBounds(ax)

# ADDING IAXO BOUNDS!!
DarkPhoton.IAXO(ax)

MySaveFig(fig,'DarkPhoton-tPlasmon-reduced')
plt.show()
