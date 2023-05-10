# Tom O'Shea 2022

# script to find the upper limit on phi for a given number of observed events
# using a poisson distribution at CL = 95%, for use in IAXO dark photon analysis

from numpy import *

# detection time in days
days = 5

# E range
dE = 0.07	# in keV

# size of sample
size = int(1e3)

# select detector
detector = 2	# 0-baby 1-baseline 2-upgraded,	-1 for CAST

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

	rng = random.default_rng()
	dist = rng.poisson(Nbg, size)

	# calculate weighted average from distro
	total = 0
	for i in dist:
		if i==0: total+=2.99573
		elif i==1: total+=4.49976
		elif i==2: total+=5.8279
		elif i==3: total+=7.0727
		else: print("get N=4 value!!")
	Nd = total / size
	
	# calculate equivalent phi
	phiLimit = Nd / ( A * effO * effD * effT * t )

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



print("\nBackground phi:	{} m-2 s-1".format(phiBg))
print("95% CL phi:	{} m-2 s-1".format(phiLimit))
