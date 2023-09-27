// Tom O'Shea 2022

// script to generate a function of dark photon mixing parameter against dark photon mass
// to give a theoretical upper limit on these parameters that IAXO could achieve

// original t-plasmon integral for IAXO vacuum run

#include "darkphoton.h"	// home of the real code
#include <thread>

using namespace std;

// argc counts number of command line inputs (including the script itself)
// argv is a vector containing cmd args separated by space (including script itself)
int main( int argc, char** argv ) {

	// check suffix is given as argument
	if( argc == 1 ){
		cout << "enter filename descriptor as argument" << endl;
		cout << "usage: " << argv[0] << " -<1 word descriptor>" << endl;
		return 1;
	}
	
	else if(argc == 2) {
	
	cout << "running full limit script (multithread) for " << argv[1] << endl;
	
	// define suffix for version/energy/whatever
	string suffix = argv[1];
	
	string isotope = "57Fe";

	// read csv files to vectors
	vector<double> r = read("data/r.dat");	// solar temperature [eV]
	vector<double> T = read("data/T.dat");	// solar temperature [eV]
	vector<double> wp = read("data/wp.dat");	// plasma frequency [eV]
	
	// initialise nuclear variables
	double w, J1, J0, tau;
	vector<double> nI;
	
	if( isotope == "57Fe" ) {
		nI = read("data/n57Fe.dat");
		// 57Fe https://www.nndc.bnl.gov/ensdf/
		w = 14.4129e3;			// [eV]
		J1 = 1.5;
		J0 = 0.5;
		tau = 98.3e-9 / s2eV;	// [eV-1]
	}
	else if ( isotope == "55Mn" ) {
		nI = read("data/n55Mn.dat");
		// 55Mn https://www.nndc.bnl.gov/ensdf/
		w = 125.949e3;			// [eV]
		J1 = 3.5;
		J0 = 2.5;
		tau = 259e-12 / s2eV;	// [eV-1]
	}
	
	else {
		cout << "invalid isotope: " << isotope << endl;
		return 3;
	}

	// babyIAXO
	string nameBaby = "babyIAXO-" + isotope + suffix;
	double L1 = 10 / m2eV;	// in eV^-1
	thread t2( nuclearFlux, r, wp, T, nI, L1, w, J1, J0, tau, nameBaby );	

	//IAXO baseline
	double L2 = 20 / m2eV;
	string nameBaseline = "baselineIAXO-" + isotope + suffix;
	//thread t4( nuclearFlux, r, wp, T, nI, L2, w, J1, J0, tau, nameBaseline );	

	// IAXO upgraded
	string nameUpgraded = "upgradedIAXO" + isotope + suffix;
	double L3 = 22 / m2eV;	// in eV^-1
	//thread t6( nuclearFlux, r, wp, T, nI, L3, w, J1, J0, tau, nameUpgraded );	

	// wait until all threads are finished
	t2.join();
	//t4.join();
	//t6.join();
	
	cout << "\nProcess completed!\n" << endl;
	
	return(0);

	}
	
	else{
		cout << "only input filename suffix" << endl;
		return 2;
	}
}
