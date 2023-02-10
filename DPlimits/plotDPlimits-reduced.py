# plotting of dark photon limits using Ciaran O'Hares plotting script
# https://github.com/cajohare/AxionLimits
# reduced version zoomed in on IAXO region

from PlotFuncs import *

fs = 20

fig,ax = DarkPhoton.FigSetup_reduced()

# # DPDM searches
DarkPhoton.Xenon(ax, text_on=False)
#DarkPhoton2.LAMPOST(ax)

# # Astrophysical bounds
DarkPhoton.StellarBounds(ax, text_on=False)

# ADDING IAXO BOUNDS!!
DarkPhoton.IAXO(ax, text_on=False)

# text
plt.text(2e-1,5e-12,r'{\bf babyIAXO}',fontsize=fs,color='black',rotation=-18,rotation_mode='anchor',ha='center',va='center', zorder=15.5)
plt.text(2e-1,2e-12,r'{\bf baselineIAXO}',fontsize=fs,color='black',rotation=-18,rotation_mode='anchor',ha='center',va='center', zorder=15.5)
plt.text(2e-1,6e-13,r'{\bf upgradedIAXO}',fontsize=fs,color='black',rotation=-18,rotation_mode='anchor',ha='center',va='center', zorder=15.5)
plt.text(2e-1,5e-11,r'{\bf Solar}',fontsize=fs,color='white',rotation=-0,rotation_mode='anchor',ha='center',va='center', zorder=15.5)
plt.text(2e-0,5e-13,r'{\bf Xenon}',fontsize=fs,color='crimson',rotation=-18,rotation_mode='anchor',ha='center',va='center', zorder=15.5)



MySaveFig(fig,'DarkPhoton-reduced')
plt.show()
