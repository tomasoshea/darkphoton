# Tom O'Shea 2022

# script to find the upper limit on phi for a given number of observed events
# using a gamma distribution, for use in IAXO dark photon analysis

from numpy import *
from scipy.stats import gamma
from matplotlib import pyplot as plt


# confidence level
CL = 0.95

# increment
dN = 1e-5

# detection time in days
days = 4

# select detector
detector = 2	# 0-baby 1-baseline 2-upgraded


def survivalGamma(b, alpha = 1 - CL):
	#Right hand side of the Gamma distribution for signal+background. Alpha in (0, 1), b number of background signals expected.
	if alpha > 1 or alpha < 0: print("Value for alpha outside bounds (0, 1)")
	a = 1.
	nsig = 0.
	while a > alpha :
		a = gamma.sf(nsig+b, b)
		nsig+=dN
	return nsig, a

def cdfGamma(b):
	#Right hand side of the Gamma distribution for signal+background. Alpha in (0, 1), b number of background signals expected.
	if CL > 1 or CL < 0: print("Value for alpha outside bounds (0, 1)")
	a = 0.
	nsig = 0.
	while a < CL :
		a = gamma.cdf(nsig+b, b)
		nsig+=dN
	return nsig, a


# decide on detector
if detector==0:
	# babyIAXO parameters
	name="babyIAXO"
	A = 0.77	# detector area [m2]
	dE = 100	# energy range [keV]
	phiBg = 1e-7 * 1e4 * dE	# background flux [m-2 s-1]
	a = 0.6 * 1e-4	# XRay detection area [m2]
	t = days * 24 * 3600	# detection time [s]
	effD = 0.7	# detectior efficiency
	effO = 0.35	# optical efficiency
	effT = 0.5	# time efficiency (proportion pointed at sun)

elif detector==1:
	# IAXO baseline parameters
	name="IAXO baseline"
	A = 2.3	# detector area [m2]
	dE = 100	# energy range [keV]
	phiBg = 1e-8 * 1e4 * dE	# background flux [m-2 s-1]
	a = 1.2 * 1e-4	# XRay detection area [m2]
	t = days * 24 * 3600	# detection time [s]
	effD = 0.8	# detectior efficiency
	effO = 0.7	# optical efficiency
	effT = 0.5	# time efficiency (proportion pointed at sun)

elif detector==2:
	# IAXO upgraded parameters
	name="IAXO upgraded"
	A = 3.9	# detector area [m2]
	dE = 100	# energy range [keV]
	phiBg = 1e-9 * 1e4 * dE	# background flux [m-2 s-1]
	a = 1.2 * 1e-4	# XRay detection area [m2]
	t = days * 24 * 3600	# detection time [s]
	effD = 0.8	# detectior efficiency
	effO = 0.7	# optical efficiency
	effT = 0.5	# time efficiency (proportion pointed at sun)



# calculate background count
Nbg = phiBg * a * t * effT
print("Detector: {}\nExpected background count: {}".format(name, Nbg))

# use gamma function
Nd, prob = survivalGamma(Nbg)
#Nd, prob2 = cdfGamma(Nbg)
#prob = 1 - prob2
print("\n		N = {}\n		CL = {}".format(Nd, 1 - prob))

# calculate equivalent phi
phiLimit = Nd / ( A * effO * effD * effT * t )

print("\nBackground phi:	{} m-2 s-1".format(phiBg))
print("{}% CL phi:	{} m-2 s-1".format(round((1 - prob)*1e2), phiLimit))

