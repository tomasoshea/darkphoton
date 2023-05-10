# Tom O'Shea 2022

# script to find the upper limit on phi for a given number of observed events
# using a gamma distribution, for use in IAXO dark photon analysis

from numpy import *
from scipy.stats import gamma
#from scipy.stats import poisson
from scipy.special import gamma as G
from matplotlib import pyplot as plt


# confidence level
CL = 0.95

# increment
dN = 1

# limit on number of iterations
limit = int(1e6)

# detection time in days
days = 120

# E range
dE = 0.1	# in keV

# select detector
detector = 0	# 0-baby 1-baseline 2-upgraded,	-1 for CAST

# homemade gamma CDF
def pdf(n, mu):
	return mu**(n) * exp(-mu) / G(n+1)

def CDF(b):
	# integrate until CDF = CL
	N = 0
	F = 0
	count = 0
	while F < CL:
		if count > limit:
			print("\nExceeded count limit ({}) !!!\n".format(limit))
			break
		if dN > b:
			print("\nReduce increment size!!!\n")
			break
		F += (dN * ( pdf(N+dN+b,b) + pdf(N+b,b)) / 2)
		N += dN
		count += 1
		#print(pdf(N+dN+b,b))
	return N, F

def survivalGamma(b, alpha = 1 - CL):
	#Right hand side of the Gamma distribution for signal+background. Alpha in (0, 1), b number of background signals expected.
	if alpha > 1 or alpha < 0: print("Value for alpha outside bounds (0, 1)")
	a = 1.
	nsig = 0.
	while a > alpha :
		a = gamma.sf(nsig+b, b)
		nsig+=dN
	return nsig, a

"""def survivalGamma(b, alpha = 1 - CL):
	#Right hand side of the Gamma distribution for signal+background. Alpha in (0, 1), b number of background signals expected.
	if alpha > 1 or alpha < 0: print("Value for alpha outside bounds (0, 1)")
	a = 1.
	nsig = 0.
	while a > alpha :
		a = poisson.sf(nsig+b, b)
		nsig+=dN
	return nsig, a"""

def cdfGamma(b):
	#Right hand side of the Gamma distribution for signal+background. Alpha in (0, 1), b number of background signals expected.
	if CL > 1 or CL < 0: print("Value for alpha outside bounds (0, 1)")
	a = 0.
	nsig = 0.
	while a < CL :
		a = gamma.cdf(nsig+b, b)
		nsig+=dN
	return nsig, a

if detector != -1:

	# decide on detector
	if detector==0:
		# babyIAXO parameters
		name="babyIAXO"
		A = 0.77	# detector area [m2]
		#dE = 100	# energy range [keV]
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
		#dE = 100	# energy range [keV]
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
		#dE = 100	# energy range [keV]
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
	#Nd, prob = survivalGamma(Nbg)
	#Nd, prob2 = cdfGamma(Nbg)
	#Nd, prob2 = CDF(Nbg)
	#prob = 1 - prob2
	#print("\n		N = {}\n		CL = {}".format(Nd, 1 - prob))
	rng = random.default_rng()
	dist = rng.poisson(0.1)
	print(dist)



	# calculate equivalent phi
	#phiLimit = Nd / ( A * effO * effD * effT * t )

# CAST option
elif detector == -1:
	name="CAST"
	A = 2.9e-3	# detector area [m2]
	dE = 15	# energy range [keV]
	phiBg = 1e-1 * dE	# background flux [m-2 s-1]
	t = 108 * 3600	# detection time [s] (108h)
	eff = 0.5	# efficiency (assume 50% for now based on nothing)

	Nbg = phiBg * A * t
	print("Detector: {}\nExpected background count: {}".format(name, Nbg))
	Nd, prob = survivalGamma(Nbg)
	print("\n		N = {}\n		CL = {}".format(Nd, 1 - prob))
	phiLimit = Nd / ( A * t * eff )



#print("\nBackground phi:	{} m-2 s-1".format(phiBg))
#print("{}% CL phi:	{} m-2 s-1".format(round((1 - prob)*1e2), phiLimit))
