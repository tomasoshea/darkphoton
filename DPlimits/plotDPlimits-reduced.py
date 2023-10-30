# plotting of dark photon limits using Ciaran O'Hares plotting script
# https://github.com/cajohare/AxionLimits
# reduced version zoomed in on IAXO region

from PlotFuncs import *

fig,ax = DarkPhoton.FigSetup(reduced=True)

# # DPDM searches
DarkPhoton.Xenon(ax,text_on=0)
DarkPhoton.DAMIC(ax,text_on=0)
DarkPhoton.SENSEI(ax,text_on=0)
DarkPhoton.SuperCDMS(ax,text_on=0)
DarkPhoton.SHUKET(ax,text_on=0)
DarkPhoton.DarkEfield(ax,text_on=0)
DarkPhoton.WISPDMX(ax,text_on=0)
DarkPhoton.SQuAD(ax,text_on=0)
DarkPhoton.DMPathfinder(ax,text_on=0)
DarkPhoton.ORPHEUS(ax,text_on=0)
DarkPhoton.DOSUE(ax,text_on=0)
DarkPhoton.FAST(ax,text_on=0)

# # Astrophysical bounds
DarkPhoton.StellarBounds(ax,text_on=0)

# ADDING IAXO BOUNDS!!
DarkPhoton.IAXO(ax,text_on=0)#, text_on=False)
#plt.text(1e-1,5e-10,r'{\bf IAXO}',fontsize=30,color='white',rotation=-20,rotation_mode='anchor',ha='center',va='center', zorder=105.5)

plt.legend()

MySaveFig(fig,'DarkPhoton-tPlasmonGas-stats-Atlas1e0-reduced3')
plt.show()
