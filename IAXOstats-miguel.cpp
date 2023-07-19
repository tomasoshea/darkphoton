// Tom O'Shea 2023

// script to find the upper limit on phi for a given number of observed events
// using a poisson distribution, for use in IAXO dark photon analysis

#include <bits/stdc++.h>

using namespace std;

double CL = 0.95;	// confidence level
double dE = 14.5;	// E range [keV]
int samplesize = 1e3;		// size of random sample

// conversion factors
double s2eV = (6.582119569e-16);	// Hz to eV
double J2eV = (1. / 1.602176634e-19);	// Joules to eV (1 / e)
double m2eV = (1.973269804e-7);	// m-1 to eV
double K2eV = (8.617333262e-5);	// Kelvin to eV
double kg2eV = 5.609588604e35;	// from hbar/c2


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


double miguel_likelihood( double Nbg ){

	double Poisson = 0;
    double P[11] = {3.0, 4.74, 6.30, 7.75, 9.15, 10.51, 11.84, 13.15, 14.43, 15.71, 16.98};
    if ( (int)Nbg < 11 ) { Poisson = P[ (int)Nbg ] - Nbg; }
    else { Poisson = Nbg; }
    return Poisson;
}


void chis( int detector ) {

	// initialise parameters
	double A, phiBg, a, t, effD, effO, effT, len;
	string name, load;
	vector<double> flux, m;

	// choose detector

	if ( detector==0 ) {
		// babyIAXO parameters
		name="babyIAXO";
		A = 0.77;	// detector area [m2]
		phiBg = 1e-7 * 1e4 * dE;	// background flux [m-2 s-1]
		a = 0.6 * 1e-4;	// XRay detection area [m2]
		t = 1.5 * 365.25 * 24 * 3600;	// 1.5 years
		effD = 0.7;	// detectior efficiency
		effO = 0.35;	// optical efficiency
		effT = 0.5;	// time efficiency (proportion pointed at sun)
		load = "data/limits/babyIAXO-tPlasmon-miguel2.dat";
		m = loadtxt(load,0);
		flux = loadtxt(load,1);
		len = flux.size();
		}

	else if ( detector==1 ) {
		// baseline IAXO parameters
		name="baselineIAXO";
		A = 2.3;	// detector area [m2]
		phiBg = 1e-8 * 1e4 * dE;	// background flux [m-2 s-1]
		a = 1.2 * 1e-4;	// XRay detection area [m2]
		t = 3 * 365.25 * 24 * 3600;	// 3 years
		effD = 0.8;	// detectior efficiency
		effO = 0.7;	// optical efficiency
		effT = 0.5;	// time efficiency (proportion pointed at sun)
		load = "data/limits/baselineIAXO-tPlasmon-miguel2.dat";
		m = loadtxt(load,0);
		flux = loadtxt(load,1);
		len = flux.size();
		}

	else if ( detector==2 ) {
		// upgraded IAXO parameters
		name="upgradedIAXO";
		A = 3.9;	// detector area [m2]
		phiBg = 1e-9 * 1e4 * dE;	// background flux [m-2 s-1]
		a = 1.2 * 1e-4;	// XRay detection area [m2]
		t = 5 * 365.25 * 24 * 3600;	// 5 years
		effD = 0.8;	// detectior efficiency
		effO = 0.7;	// optical efficiency
		effT = 0.5;	// time efficiency (proportion pointed at sun)
		load = "data/limits/upgradedIAXO-tPlasmon-miguel2.dat";
		m = loadtxt(load,0);
		flux = loadtxt(load,1);
		len = flux.size();
		}

	// calculate background count
	double Nbg = phiBg * a * t * effT;
	cout << "Detector: "<< name << endl;
	cout << "Expected background count: " << Nbg << endl << endl;

	// get poisson random number
	default_random_engine generator;
	poisson_distribution<int> distro( Nbg );

	vector<double> chi;	// initialise chi vector

	// get 95% chi from miguel's poisson
	for ( int c = 0; c < len; c++ ) {
		
		// signal flux for chi = 1 for small dt
		double Nsig = ( flux[c] * pow(m2eV, -2) / s2eV ) * A * effO * effD * effT * t;
		double total = 0;	// keep total of all runs	

		double n = miguel_likelihood(Nbg);
		//cout << n << "		" << Nbg << "		" << Nsig << endl;
		double calcdChi = pow( n / Nsig , 0.25);

		chi.push_back( calcdChi );
		cout << c+1 << " out of " << len << ":	" << calcdChi << endl;
	}
	
	//cout << "chi length: " << chi.size() << "	m length: " << m.size() << endl;
	// write out
	string savename = "data/limits/stats-" + name + "-tPlasmon-miguel4.dat";
	write2D( savename, m, chi );
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
