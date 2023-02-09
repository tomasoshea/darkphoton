# plotting of dark photon limits using Ciaran O'Hares plotting script
# https://github.com/cajohare/AxionLimits
# reduced version zoomed in on IAXO region


from PlotFuncs import *
%matplotlib inline

fig,ax = DarkPhoton2.FigSetup()

# # DPDM searches
DarkPhoton2.Xenon(ax)
#DarkPhoton2.LAMPOST(ax)

# # Astrophysical bounds
DarkPhoton2.StellarBounds(ax)

# ADDING IAXO BOUNDS!!
DarkPhoton2.IAXO(ax)

MySaveFig(fig,'DarkPhoton-reduced-tPlasmon-5')

