# plotting of dark photon limits using Ciaran O'Hares plotting script
# https://github.com/cajohare/AxionLimits
# reduced version zoomed in on IAXO region

from PlotFuncs import *

fig,ax = DarkPhoton.FigSetup(reduced=True)

# # DPDM searches
DarkPhoton.Xenon(ax,text_on=False)
DarkPhoton.DAMIC(ax,text_on=False)
DarkPhoton.SENSEI(ax,text_on=False)
DarkPhoton.SuperCDMS(ax,text_on=False)
DarkPhoton.SHUKET(ax,text_on=False)
DarkPhoton.DarkEfield(ax,text_on=False)
DarkPhoton.WISPDMX(ax,text_on=False)
DarkPhoton.SQuAD(ax,text_on=False)
DarkPhoton.DMPathfinder(ax,text_on=False)
DarkPhoton.ORPHEUS(ax,text_on=False)
DarkPhoton.DOSUE(ax,text_on=False)
DarkPhoton.FAST(ax,text_on=False)

# # Astrophysical bounds
DarkPhoton.StellarBounds(ax,text_on=False)

# ADDING IAXO BOUNDS!!
DarkPhoton.IAXO(ax,text_on=False)
#plt.text(1e-1,5e-11,r'{\bf IAXO}',fontsize=20,color='white',rotation=-32,rotation_mode='anchor',ha='center',va='center', zorder=105.5)
#plt.text(1e-1,5e-10,r'{\bf IAXO}',fontsize=30,color='white',rotation=-20,rotation_mode='anchor',ha='center',va='center', zorder=105.5)

ax.legend()
MySaveFig(fig,'DarkPhoton-AtlasGas-30eV-5yr-reduced')
plt.show()
