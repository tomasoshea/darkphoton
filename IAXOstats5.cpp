// Tom O'Shea 2023

// script to find the upper limit on phi for a given number of observed events
// using a poisson distribution, for use in IAXO dark photon analysis

#include <bits/stdc++.h>

using namespace std;

double CL = 0.95;	// confidence level
double days = 500;	// detection time
double dE = 0.07;	// E range [keV]
int samplesize = 1e2;		// size of random sample


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


// write out 2 column datafile
void write2D( string name, vector<double> data1, vector<double> data2) {

	// delete old file if extant
	if ( remove(name.c_str()) == 0 ) {
		cout << "file " << name << " deleted." << endl;
		}
	
	// create new file
	fstream fout;
	fout.open(name, ios::out | ios::trunc);
	
	// get size of vectors
	int len1 = data1.size();
	int len2 = data2.size();
	
	// check vectors of equal length
	if ( len1 == len2 ) {
		
		// Read the input from vector
		for ( int c = 0; c < len1; c++ ) {
			// Insert the data to file
			fout << data1[c] << "	" << data2[c] << endl;
		}
		
		cout << "file " << name << " created succesfully." << endl;
	}
	
	// if vectors not equal return error
	else { cout << "ERROR - ensure vectors are of equal length" << endl; }
	
	fout.close();
}


// factorial
int factorial( int N ) {

	int total = 1;
	for ( int c = 1; c <= N; c++ ) { total *= c; }
	return total;
}


// 1/2 chi2 for 1 density step
double L( double x, double b, double s, int n, double m, vector<double> ne,
	vector<double> T, vector<double> wp, vector<double> r,
    vector<double> nH, vector<double> nHe4, vector<double> nHe3,
    double L, vector<vector<double>> z1, vector<vector<double>> z2  ) {

	double s = pow(x,4) * 
        integrateGas( m, ne, T, wp, r, nH, nHe4, nHe3, L, z1, z2 );
	
	double item;

	if ( n == 0 ) { item = ( b + pow(x,4)*s ); }
	else if ( n == 1 ) { item = ( ( b + pow(x,4)*s ) - log( b + pow(x,4)*s ) - 1 ); }
	else { item = ( ( b + pow(x,4)*s ) - n * ( log( b + pow(x,4)*s ) + 1 - log(n) ) ); }

	return exp(-item);
}


// integral for getting CL
double integral( double b, vector<double> s, vector<double> n ) {

	double x = 1e-10;
	double x2 = x;
	double mx = 1.01;
	//double dx = x;
	//cout << x2 << endl;
	//double total = 0.5 * x * ( exp( - f( x, b, s, n ) ) + exp( - f( 0, b, s, n ) ) );
	double total= 0.5 * x * ( L( x, b, s, n ) + L( 0, b, s, n ) );
	double norm = total;

	// normalise with intg to inf
	while(true) {
		double dx = x2 * (mx - 1);
		double L1 = L( x2, b, s, n );
		double L2 = L( x2*mx, b, s, n );
		//double L2 = L( x2+dx, b, s );
		if ( isnan(L1 + L2) ) { continue; }
		else { norm += 0.5 * dx * ( L1 + L2 ); }
		if ( L2 + L1 == 0 ) { break; }
		//x2 += dx;
		x2 *= mx;
		//cout << L1 << endl;
	}

	// integrate up to CL
	while ( ( total / norm ) < CL ) {
		double dx = x * (mx - 1);
		double L1 = L( x, b, s, n );
		//double L2 = L( x+dx, b, s );
		double L2 = L( x*mx, b, s, n );
		if ( isnan(L1 + L2) ) { continue; }
		else { total += 0.5 * dx * ( L1 + L2 ); }
		//cout << total * dx / norm << endl;
		x *= mx;
		//x2 += dx;
	}
	return x;
}


void chis( int detector ) {

	// import shite
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
	double A, phiBg, a, t, effD, effO, effT, len;
	string name;

	// choose detector

	if ( detector==0 ) {
		// babyIAXO parameters
		name="babyIAXO";
		A = 0.77;	// detector area [m2]
		phiBg = 1e-7 * 1e4 * dE;	// background flux [m-2 s-1]
		a = 0.6 * 1e-4;	// XRay detection area [m2]
		t = days * 24 * 3600;	// detection time [s]
		effD = 0.7;	// detectior efficiency
		effO = 0.35;	// optical efficiency
		effT = 0.5;	// time efficiency (proportion pointed at sun)
		m = loadtxt("data/limits/babyIAXO-tPlasmonflux-gas.dat",0);
		flux = loadtxt("data/limits/babyIAXO-tPlasmonflux-gas.dat",1);
		len = flux.size();
		}

	else if ( detector==1 ) {
		// baseline IAXO parameters
		name="baselineIAXO";
		A = 2.3;	// detector area [m2]
		phiBg = 1e-8 * 1e4 * dE;	// background flux [m-2 s-1]
		a = 1.2 * 1e-4;	// XRay detection area [m2]
		t = days * 24 * 3600;	// detection time [s]
		effD = 0.8;	// detectior efficiency
		effO = 0.7;	// optical efficiency
		effT = 0.5;	// time efficiency (proportion pointed at sun)
		m = loadtxt("data/limits/baselineIAXO-tPlasmonflux-gas.dat",0);
		flux = loadtxt("data/limits/baselineIAXO-tPlasmonflux-gas.dat",1);
		len = flux.size();
		}

	else if ( detector==2 ) {
		// upgraded IAXO parameters
		name="upgradedIAXO";
		A = 3.9;	// detector area [m2]
		phiBg = 1e-9 * 1e4 * dE;	// background flux [m-2 s-1]
		a = 1.2 * 1e-4;	// XRay detection area [m2]
		t = days * 24 * 3600;	// detection time [s]
		effD = 0.8;	// detectior efficiency
		effO = 0.7;	// optical efficiency
		effT = 0.5;	// time efficiency (proportion pointed at sun)
		m = loadtxt("data/limits/upgradedIAXO-tPlasmonflux-gas.dat",0);
		flux = loadtxt("data/limits/upgradedIAXO-tPlasmonflux-gas.dat",1);
		len = flux.size();
		}

	// calculate background count
	double Nbg = phiBg * a * t * effT;
	cout << "Detector: "<< name << endl;
	cout << "Expected background count: " << Nbg << endl << endl;

	// get poisson random number
	vector<double> n;
	default_random_engine generator;
	poisson_distribution<int> distro( Nbg );
	for ( int c = 0; c < samplesize; c++ ) { n.push_back( distro(generator) ); }


	//vector<double> chi;	// initialise chi vector
			
	// signal flux for chi = 1 for small dt
	vector<double> Nsig;
	for ( int c = 0; c < flux.size(); c++ ) { Nsig.push_back( flux[c] * A * effO * effD * effT * t ); }

	double total = 0;
	for ( int j = 0; j < samplesize; j++ ) {
		for ( int c = 0; c < flux.size(); c++ ) { n.push_back( distro(generator) ); }
		total += integral( Nbg, Nsig, n );
	}
			
		//chi.push_back( total );
		//cout << c+1 << " out of " << len << ":	" << total << endl;
	//cout << "chi length: " << chi.size() << "	m length: " << m.size() << endl;
	// write out
	string savename = "data/limits/stats-" + name + "5.dat";
	//write2D( savename, m, chi );
}



int main(){
	
	// thread all 3 at same time
	thread t1(chis, 0); usleep(100);	// baby
//	thread t2(chis, 1); usleep(100);	// baseline
//	thread t3(chis, 2); usleep(100);	// upgraded
	
	t1.join();
//	t2.join();
//	t3.join();
	
	cout << "\n¡¡complete!!" << endl;
	return 0;
}
