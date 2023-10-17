// Tom O'Shea 2023

// script to find the upper limit on phi for a given number of observed events
// using a poisson distribution, for use in IAXO dark photon analysis

#include <bits/stdc++.h>
#include "darkphoton.h"

using namespace std;

double CL = 0.95;	// confidence level
double days = 5;	// detection time
double dE = 0.3;	// E range [keV]
int samplesize = 10;		// size of random sample


// read 2nd column from 2 column datafile
vector<double> loadtxt( string name, int col ) {

	//cout << "Reading file " << name << "..." << endl;

	// open file defined in argument
	fstream file;
	file.open(name,ios::in);
	
	char delim('	');	// define delimiter for file parsing (tab)
	
	if ( file.is_open() ){   // checking whether the file is open
		string temp;	// define temporary storage string
		vector<double> row1, row2;	// define vector to store input values and return
		vector<string> v;
		
		while( getline(file, temp) ){ v.push_back( temp ); }

		for ( string i : v ) {

			stringstream X(i);
			string temp;
			vector<string> vec;
			while ( getline(X, temp, delim ) ) { vec.push_back(temp); }
			row1.push_back( stod(vec[0]) );
			row2.push_back( stod(vec[1]) );
		}
		
	file.close();   // close the file object.
	
	// choose column
	if ( col == 0 ) { return row1; }
	else if ( col == 1 ) { return row2; }
	else { cout << "put 0 or 1 for columns" << endl; return {69.}; }	
	}
	else{ cout << "file " << name << " doesn't exist" << endl; return {69.}; }
}


// factorial
int factorial( int N ) {

	int total = 1;
	for ( int c = 1; c <= N; c++ ) { total *= c; }
	return total;
}


double like( double x, double n, double m, vector<double> ne, vector<double> T, vector<double> wp,
	 vector<double> r, vector<double> nH, vector<double> nHe4, vector<double> nHe3, vector<vector<double>> z1, vector<vector<double>> z2,
	 double A, double phiBg, double area, double t, double effD, double effO, double effT, double len, double b, double L ) {

	double item = 0;
	
	for ( double wpIAXO = 1e-2; wpIAXO <= 1e-1; wpIAXO += 2e-2 ) {	// run over various densities
	
		double s = integrateGasFlux( m, wpIAXO, ne, T, wp, r, nH, nHe4, nHe3, L, z1, z2 );
		s *= ( pow(m2eV, -2) / s2eV * A * effO * effD * effT * t );

		if ( n == 0 ) { item += ( b + pow(x,4)*s ); }
		else if ( n == 1 ) { item += ( ( b + pow(x,4)*s ) - log( b + pow(x,4)*s ) - 1 ); }
		else { item += ( ( b + pow(x,4)*s ) - n * ( log( b + pow(x,4)*s ) + 1 - log(n) ) ); }

	}
	//cout << item << endl;
	//cout << item << endl;
	return exp(-item);
}


// integral for getting CL
double integral( double n, double m, vector<double> ne, vector<double> T, vector<double> wp,
	 vector<double> r, vector<double> nH, vector<double> nHe4, vector<double> nHe3, vector<vector<double>> z1, vector<vector<double>> z2,
	 double A, double phiBg, double area, double t, double effD, double effO, double effT, double len, double b, double L ) {

	double x = 1e-12;
	double x2 = x;
	//double dx = x;
	//cout << x2 << endl;
	double mx = 1.1;
	//double total = 0.5 * x * ( exp( - f( x, b, s, n ) ) + exp( - f( 0, b, s, n ) ) );
	double total= 0.5 * x * ( like( x, n, m, ne, T, wp, r, nH, nHe4, nHe3, z1, z2, A, phiBg, area, t, effD, effO, effT, len, b, L ) 
					+ like( 0, n, m, ne, T, wp, r, nH, nHe4, nHe3, z1, z2, A, phiBg, area, t, effD, effO, effT, len, b, L ) );
	double norm = total;

	// normalise with intg to inf
	while(true) {
		double dx = x2 * (mx - 1);
		double L1 = like( x2, n, m, ne, T, wp, r, nH, nHe4, nHe3, z1, z2, A, phiBg, area, t, effD, effO, effT, len, b, L );
		//double L2 = L( x2*mx, b, s, n );
		double L2 = like( x2+dx, n, m, ne, T, wp, r, nH, nHe4, nHe3, z1, z2, A, phiBg, area, t, effD, effO, effT, len, b, L );
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
		double L1 = like( x, n, m, ne, T, wp, r, nH, nHe4, nHe3, z1, z2, A, phiBg, area, t, effD, effO, effT, len, b, L );
		double L2 = like( x+dx, n, m, ne, T, wp, r, nH, nHe4, nHe3, z1, z2, A, phiBg, area, t, effD, effO, effT, len, b, L );
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
	vector<double> r = read("data/r.dat");	// sun radial distance [eV-1]
	vector<double> rFrac = read("data/rFrac.dat");	// sun radial distance as fraction
	vector<double> T = read("data/T.dat");	// solar temperature [eV]
	vector<double> ne = read("data/ne.dat");	// electron number density [eV3]
	vector<double> wp = read("data/wp.dat");	// plasma frequency [eV]
	vector<double> nH = read("data/nH.dat");	// H ion density [eV3]
	vector<double> nHe4 = read("data/nHe4.dat");	// He4 ion density [eV3]
	vector<double> nHe3 = read("data/nHe3.dat");	// He3 ion density [eV3]
	
	// get gaunt factors
	vector<vector<double>> z1 = readGaunt("data/Z1.dat");	// gaunt factors for Z=1
	vector<vector<double>> z2 = readGaunt("data/Z2.dat");	// gaunt factors for Z=2
	
	// convert Gaunt factor Theta to T in eV
	for( int i = 1; i < 201; i++ ) { z1[0][i] = z1[0][i] * m_e; }
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

	// calculate background count
	double Nbg = phiBg * area * t * effT;
	cout << "Detector: "<< name << endl;
	cout << "Expected background count: " << Nbg << endl << endl;

	// get poisson random number
	default_random_engine generator;
	poisson_distribution<int> distro( Nbg );

	vector<double> chi;	// initialise chi vector

	// get 95% chi for each m value my minimisation
	for ( double m = 1e-2; m <= 1e-1; m += 2e-2 ) {
		
		double total = 0;	// keep total of all runs
		
		// get sample of random N from poisson
		for ( int i = 0; i < samplesize; i++ ) {
			double n = distro(generator);	// get random n from poisson
			//if ( n == 0 ) { continue; }
			//else { total += min( Nbg, Nsig, n ); }
			total += integral(n, m, ne, T, wp, r, nH, nHe4, nHe3, z1, z2, A, phiBg, area, t, effD, effO, effT, len, b, L );
		}
			
		chi.push_back( total / samplesize );
		mvec.push_back(m);
		cout << (int)(m*1e2) << " out of " << 10 << ":	" << total / samplesize << endl;
	}
	
	//cout << "chi length: " << chi.size() << "	m length: " << m.size() << endl;
	// write out
	string savename = "data/limits/newstats-5yr-2-" + name + "-tPlasmonGas.dat";
	write2D( savename, mvec, chi );
}



int main(){
	
	// thread all 3 at same time
	thread t1(chis, 0); usleep(100);	// baby
	thread t2(chis, 1); usleep(100);	// baseline
	thread t3(chis, 2); usleep(100);	// upgraded
	
	t1.join();
	t2.join();
	t3.join();
	
	cout << "\n¡¡complete!!" << endl;
	return 0;
}
