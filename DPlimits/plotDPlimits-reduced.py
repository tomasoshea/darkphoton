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
DarkPhoton.IAXO(ax)#, text_on=False)
#plt.text(1e-1,5e-10,r'{\bf IAXO}',fontsize=30,color='white',rotation=-20,rotation_mode='anchor',ha='center',va='center', zorder=105.5)

MySaveFig(fig,'DarkPhoton-tPlasmon-clean2-reduced')
plt.show()
