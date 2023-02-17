# convert pressure from eV4 to Pa

J2eV = (1. / 1.602176634e-19)   # J => eV
m2eV = (1.973269804e-7) # m-1 => eV
Pa2atm = 9.86923e-6 # Pa => atm

peV = 6.4e9   # eV4
pPa = peV / (J2eV * m2eV**3)

print(pPa/1e9 , "e9 Pa")
print(pPa* Pa2atm/1e6 , "e6 atm")