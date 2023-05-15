// Tom O'Shea 2023

// script to find the upper limit on phi for a given number of observed events
// using a poisson distribution, for use in IAXO dark photon analysis

#include <bits/stdc++.h>

using namespace std;

double CL = 0.95;	// confidence level
double days = 5;	// detection time
double dE = 0.07;	// E range [keV]
int samplesize = 1e3;		// size of random sample


// read 2nd column from 2 column datafile
vector<double> loadtxt( string name, int col ) {

	//cout << "Reading file " << name << "..." << endl;

	// open file defined in argument
	fstream file;
	file.open(name,ios::in);
	
	char delim('	');	// define delimiter for file parsing (tab)
	
	if (file.is_open()){   // checking whether the file is open
		string temp;	// define temporary storage string
		vector<double> row1, row2;	// define vector to store input values and return
		
		int c = 0;	// use to separate columns
		
		while(getline(file, temp, delim)){  // put separated data into temp
			double rem = remainder(c,2);	// select which vector
			double item = stod(temp);	// convert string to double
			// divide columns
			if( rem == 0 ) { row1.push_back(item); }
			else { row2.push_back(item); }
			c++;
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
		

// poisson distro
double poisson( int N, double m ){
	return pow(m,N) * exp(-m) / factorial(N);
}


// find mu for given N in poissonian
double get_mu( int N ){
	double p = 0.9;
	double m = N;
	double dm = 1e-3;
	double dp = 1e-3;
	while( (p - (1-CL)) > dp) {
		p = poisson( N, m );
		m += dm;
	}
	return m;
}


// function to minimise
double f( double x, double b, double s, int n ){

	double mu = b + ( pow(x,4) * s );
	if( n == 0 ) { return mu; }
	else { return mu - ( n * ( log(mu) + 1 - log(n) ) ); }
}


// minimisation fn
double min( double b, double s, int n ){

	double x = 1e-1;	// starting chi
	int lim = 1e6;		// no. iterations cutoff
	bool done = false;	// completion marker
	double save = 0;	// save final result here
	int c = 0;	// for cutoff

	while( true ) {

		double val = f( x, b, s, n );
		double val2 = f( x/2, b, s, n );
		if( val2 < val ) { x /= 2; }
		else{ break; }
	}
	
	double dx = x/1000;
	while( true ) {
		
		c++;
		double minus = f( x-dx, b, s, n );
		double plus = f( x+dx, b, s, n );
		if( minus < plus ) { x -= dx; }
		else if( plus < minus ) { x += dx; }
		else { break; }
		
		if ( c > lim ) { break; }
	}
	//cout << x << endl;
	return x;
}

void chis( int detector ) {

	// initialise parameters
	double A, phiBg, a, t, effD, effO, effT, len;
	string name;
	vector<double> flux, m;

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
	default_random_engine generator;
	poisson_distribution<int> distro(Nbg);

	vector<double> chi;	// initialise chi vector

	// get 95% chi for each m value my minimisation
	for ( int c = 0; c < len; c++ ) {
		
		// signal flux for chi = 1
		double Nsig = flux[c] * A * effO * effD * effT * t;
		double total = 0;	// keep total of all runs
		
		// get sample of random N from poisson
		for ( int i = 0; i < samplesize; i++ ) {
			int n = distro(generator);	// get random n from poisson
			//cout << n << endl;
			//if ( n != 0 ) { cout << n << endl; }
			total += min( Nbg, Nsig, n );
		}
			
		chi.push_back( total / samplesize );
		//cout << c+1 << " out of " << len << ":	" << total / sample << endl;
	}
	
	//cout << "chi length: " << chi.size() << "	m length: " << m.size() << endl;
	// write out
	string savename = "data/stats-" + name + ".dat";
	write2D( savename, m, chi );
}



int main(){
	
	// thread all 3 at same time
	thread t1(chis, 0); usleep(1);	// baby
	thread t2(chis, 1); usleep(1);	// baseline
	thread t3(chis, 2); usleep(1);	// upgraded
	
	t1.join();
	t2.join();
	t3.join();
	
	cout << "\n¡¡Complete!!" << endl;
	return 0;
}
