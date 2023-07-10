// Tom O'Shea 2023

// script to find the upper limit on phi for a given number of observed events
// using a poisson distribution, for use in IAXO dark photon analysis

#include <bits/stdc++.h>
#include "darkphoton.h"

using namespace std;

double CL = 0.95;	// confidence level
double days = 10;	// detection time
double dE = 0.07;	// E range [keV]
int samplesize = 2;	// size of random sample


// factorial
int factorial( int N ) {

	int total = 1;
	for ( int c = 1; c <= N; c++ ) { total *= c; }
	return total;
}


double backConversion( double m, double wpIAXO, vector<vector<double>> z2, vector<double> omegas, vector<double> fluxes, double L ) {

	double total = 0;	// initiate value of sum at 0
	for ( int c = 0; c < omegas.size(); c++ ) {

		if ( omegas[c] > m + 1e3 ) { continue; }	// set integral cutoff
		if ( omegas[c] <= m ) { continue; }	// only allow when energy greater than mass		

		double dA;
		if( m == wpIAXO ) {
			dA = 0.5 * pow(m,4) * ( ( PgasFull( omegas[c+1], m, L, z2 ) * fluxes[c+1] ) 
						+ ( PgasFull( omegas[c], m, L, z2 ) * fluxes[c] ) );
			dA *= (omegas[c+1] - omegas[c]);
		} else {
			dA = 0.5 * pow(m,4) * ( ( PgasFlux( omegas[c+1], m, wpIAXO, L, z2 ) * fluxes[c+1] ) 
						+ ( PgasFlux( omegas[c], m, wpIAXO, L, z2 ) * fluxes[c] ) );
			dA *= (omegas[c+1] - omegas[c]);
		}
		
		// only add if real
		if ( isnan(dA) ) { continue; }
		total += abs(dA);
	}

	return total;
}


double like( double x, double b, double n, double m, vector<vector<double>> z2,
			vector<double> omegas, vector<double> fluxes, double A,
				double phiBg, double area, double t, double effD, double effO,
				double effT, double len, double L ) {

	double item = 0;

	//double wpMultiple = 1.098541;	// 50
	double wpMultiple = 1.0476;		// 100
	
	for ( double wpIAXO = 1e-2; wpIAXO <= 1e0; wpIAXO *= wpMultiple ) {	// run over various densities

		// only do for neighbouring density steps
		if( abs(wpIAXO - m) > wpIAXO * wpMultiple ) { continue; }
	
		double s = backConversion( m, wpIAXO, z2, omegas, fluxes, L );
		
		s *= ( pow(m2eV, -2) / s2eV * A * effO * effD * effT * t );

		if ( n == 0 ) { item += ( b + pow(x,4)*s ); }
		else if ( n == 1 ) { item += ( ( b + pow(x,4)*s ) - log( b + pow(x,4)*s ) - 1 ); }
		else { item += ( ( b + pow(x,4)*s ) - n * ( log( b + pow(x,4)*s ) + 1 - log(n) ) ); }

	}
	return exp(-item);
}


// integral for getting CL
double integral( double b, double n, double m, vector<vector<double>> z2,
				vector<double> omegas, vector<double> fluxes, double A,
				double phiBg, double area, double t, double effD, double effO,
				double effT, double len, double L ) {

	double x = 1e-12;
	double x2 = x;
	//double dx = x;
	//cout << x2 << endl;
	double mx = 1.1;
	//double total = 0.5 * x * ( exp( - f( x, b, s, n ) ) + exp( - f( 0, b, s, n ) ) );
	double total= 0.5 * x * ( like( x, b, n, m, z2, omegas, fluxes, A, phiBg, area, t, effD, effO, effT, len, L )
				+ like( 0, b, n, m, z2, omegas, fluxes, A, phiBg, area, t, effD, effO, effT, len, L ) );
	double norm = total;

	// normalise with intg to inf
	while(true) {
		double dx = x2 * (mx - 1);
		double L1 = like( x2, b, n, m, z2, omegas, fluxes, A, phiBg, area, t, effD, effO, effT, len, L );
		//double L2 = L( x2*mx, b, s, n );
		double L2 = like( x2+dx, b, n, m, z2, omegas, fluxes, A, phiBg, area, t, effD, effO, effT, len, L );
		if ( isnan(L1 + L2) ) { continue; }
		else { norm += 0.5 * dx * ( L1 + L2 ); }
		if ( L2 + L1 == 0 ) { break; }
		//cout << "x2 = " << x2 << "	L1 + L2 = " << L1 + L2 << endl;
		//x2 += dx;
		x2 *= mx;
	}

	// integrate up to CL
	while ( ( total / norm ) < CL ) {
		double dx = x * (mx - 1);
		double L1 = like( x, b, n, m, z2, omegas, fluxes, A, phiBg, area, t, effD, effO, effT, len, L );
		double L2 = like( x+dx, b, n, m, z2, omegas, fluxes, A, phiBg, area, t, effD, effO, effT, len, L );
		//double L2 = L( x*mx, b, s, n );
		if ( isnan(L1 + L2) ) { continue; }
		else { total += 0.5 * dx * ( L1 + L2 ); }
		//cout << total * dx / norm << endl;
		x *= mx;
		//x2 += dx;
		//cout << x << endl;
	}
	//cout << x << endl;
	return x;
}


void chis( int detector ) {

	// read csv files to vectors
	vector<vector<double>> z2 = readGaunt("data/Z2.dat");	// gaunt factors for Z=2
	string loadname = "data/solarflux-70eV.dat";
	vector<double> omegas = loadtxt(loadname, 0);
	vector<double> fluxes = loadtxt(loadname, 1);
	
	// convert Gaunt factor Theta to T in eV
	for( int i = 1; i < 201; i++ ) { z2[0][i] = z2[0][i] * m_e; }

	// initialise parameters
	double A, phiBg, area, t, effD, effO, effT, len, b, L;
	vector<double> flux, mvec;
	string name;

	// choose detector

	if ( detector==0 ) {
		// babyIAXO parameters
		name="babyIAXO";
		A = 0.77;	// detector area [m2]
		phiBg = 1e-7 * 1e4 * dE;	// background flux [m-2 s-1]
		area = 0.6 * 1e-4;	// XRay detection area [m2]
		t = days * 24 * 3600;	// detection time [s]
		effD = 0.7;	// detectior efficiency
		effO = 0.35;	// optical efficiency
		effT = 0.5;	// time efficiency (proportion pointed at sun)
		L = 10 / m2eV;		// bore length [eV-1]
		}

	else if ( detector==1 ) {
		// baseline IAXO parameters
		name="baselineIAXO";
		A = 2.3;	// detector area [m2]
		phiBg = 1e-8 * 1e4 * dE;	// background flux [m-2 s-1]
		area = 1.2 * 1e-4;	// XRay detection area [m2]
		t = days * 24 * 3600;	// detection time [s]
		effD = 0.8;	// detectior efficiency
		effO = 0.7;	// optical efficiency
		effT = 0.5;	// time efficiency (proportion pointed at sun)
		L = 20 / m2eV;		// bore length [eV-1]
		}

	else if ( detector==2 ) {
		// upgraded IAXO parameters
		name="upgradedIAXO";
		A = 3.9;	// detector area [m2]
		phiBg = 1e-9 * 1e4 * dE;	// background flux [m-2 s-1]
		area = 1.2 * 1e-4;	// XRay detection area [m2]
		t = days * 24 * 3600;	// detection time [s]
		effD = 0.8;	// detectior efficiency
		effO = 0.7;	// optical efficiency
		effT = 0.5;	// time efficiency (proportion pointed at sun)
		L = 22 / m2eV;		// bore length [eV-1]
		}

	//cout << "detector = " << detector << endl;
	//cout << "L = " << L << endl;
	//cout << "t = " << t << endl;
	//cout << "area = " << A << endl;

	// calculate background count
	double Nbg = phiBg * area * t * effT;
	cout << "Detector: "<< name << endl;
	cout << "Expected background count: " << Nbg << endl << endl;

	// get poisson random number
	default_random_engine generator;
	poisson_distribution<int> distro( Nbg );

	vector<double> chi;	// initialise chi vector

	// get 95% chi for each m value my minimisation
	for ( double m = 1e-2; m <= 1e0; m *= 1.01 ) {
		
		double total = 0;	// keep total of all runs
		
		// get sample of random N from poisson
		for ( int i = 0; i < samplesize; i++ ) {
			double n = distro(generator);	// get random n from poisson
			//if ( n == 0 ) { continue; }
			//else { total += min( Nbg, Nsig, n ); }
			total += integral( b, n, m, z2, omegas, fluxes, A, phiBg, area, t, effD, effO, effT, len, L );
		}
			
		chi.push_back( total / samplesize );
		mvec.push_back(m);
		cout << "m = " << m << ":	" << total / samplesize << endl;
	}
	
	//cout << "chi length: " << chi.size() << "	m length: " << m.size() << endl;
	// write out
	string savename = "data/limits/newstats-70eV-100" + name + "-tPlasmonGas.dat";
	write2D( savename, mvec, chi );
}


int main(){
	
	// thread all 3 at same time
	thread t1(chis, 0); usleep(10000);	// baby
	thread t2(chis, 1); usleep(10000);	// baseline
	thread t3(chis, 2); usleep(10000);	// upgraded
	
	t1.join();
	t2.join();
	t3.join();
	
	cout << "\n¡¡complete!!" << endl;
	return 0;
}
