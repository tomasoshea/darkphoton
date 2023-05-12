# Tom O'Shea 2022

# script to find the upper limit on phi for a given number of observed events
# using a poisson distribution, for use in IAXO dark photon analysis

from numpy import *
from scipy.special import factorial
from scipy.optimize import minimize

# confidence level
CL = 0.95

# detection time in days
days = 5

# E range
dE = 0.07	# in keV

# size of sample
size = int(1e3)

# select detector
detector = 0	# 0-baby 1-baseline 2-upgraded,	-1 for CAST

# poisson distro
def poiss( N, m ):
	return m**(N) * exp(-m) / factorial(N)

# find mu for given N in poissonian
def get_mu( N ):
	p = 1
	m = N
	dm = 1e-3
	dp = 1e-3
	while (p - (1-CL)) > dp:
		p = poiss( N, m )
		m += dm
	return m

# function to minimise
def min( x, b, n, N ):
	mu = b + ( (x**4) * N )
	return mu - ( n * ( log(mu) + 1 - log(n) ) )

# function to minimise if n=0
def min0( x, b, N ):
	return b + ( (x**4) * N )


# import fluxes
dat = loadtxt("data/limits/babyIAXO-tPlasmonflux-gas.dat")
m = dat[:,0]
flux = dat[:,1]

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

	# setup chi array
	chi = []
	
	for i in range(len(m)):

		# calculate chi-less signal count
		Nsig = flux[i] * A * effO * effD * effT * t

		# minimise -ln(L) using dist and Nsig
		total = 0
		for n in dist:
			if n == 0: total += minimize( min0, 1e5, (Nbg, Nsig ) )['x'][0]
			else: total += minimize( min, 1e5, (Nbg, n, Nsig ) )['x'][0]
		chi.append( total / len(dist) )
		print(i)

print(chi)
