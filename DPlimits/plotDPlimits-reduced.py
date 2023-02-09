# plotting of dark photon limits using Ciaran O'Hares plotting script
# https://github.com/cajohare/AxionLimits
# reduced version zoomed in on IAXO region


from PlotFuncs import *

fig,ax = DarkPhoton.FigSetup_reduced()

# # DPDM searches
DarkPhoton.Xenon(ax, text_on=False)
#DarkPhoton2.LAMPOST(ax)

# # Astrophysical bounds
DarkPhoton.StellarBounds(ax, text_on=False)

# ADDING IAXO BOUNDS!!
DarkPhoton.IAXO(ax)

MySaveFig(fig,'DarkPhoton-reduced')
plt.show()

