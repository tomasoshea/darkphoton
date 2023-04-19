// Tom O'Shea 2022

// script to generate a function of dark photon mixing parameter against dark photon mass
// to give a theoretical upper limit on these parameters that IAXO could achieve

// t-plasmon integral for IAXO gas run 

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

	// read csv files to vectors
	vector<double> r = read("data/r.dat");	// sun radial distance [eV-1]
	vector<double> rFrac = read("data/rFrac.dat");	// sun radial distance as fraction
	vector<double> T = read("data/T.dat");	// solar temperature [eV]
	vector<double> n = read("data/ne.dat");	// electron number density [eV3]
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
	
	// babyIAXO
	
	double Lraw = 10;	// m
	string nameBaby = "babyIAXO-tPlasmon" + suffix;
//	double phiRaw = 1.04e-4;	// m-2 s-1 (1 day)
	double phiRaw = 5.57e-5;	// m-2 s-1 (4 days)

	// convert raw values to eV
	double phi = phiRaw * m2eV * m2eV * s2eV;	// in eV^3 per keV
	double L = Lraw / m2eV;	// in eV^-1
	
	// multithread to run simultaneously
	thread t2( integrateTgas, n, T, wp, r, nH, nHe4, nHe3, L, z1, z2, phi*factor, nameBaby );	

	
	//IAXO baseline
	Lraw = 20;	// m
	string nameBaseline = "baselineIAXO-tPlasmon" + suffix;
	phiRaw = 3.24e-6;	// m-2 s-1 (4 days)

	// convert raw values to eV
	phi = phiRaw * m2eV * m2eV * s2eV;	// in eV^3 per keV
	L = Lraw / m2eV;	// in eV^-1

	// multithread to run simultaneously
	thread t4( integrateTgas, n, T, wp, r, nH, nHe4, nHe3, L, z1, z2, phi*factor, nameBaseline );

	/// IAXO upgraded
	Lraw = 22;	// m
	string nameUpgraded = "upgradedIAXO-tPlasmon" + suffix;
	phiRaw = 2.409e-08;	// m-2 s-1 (4 days)

	// convert raw values to eV
	phi = phiRaw * m2eV * m2eV * s2eV;	// in eV^3 per keV
	L = Lraw / m2eV;	// in eV^-1

	// multithread to run simultaneously
	thread t6( integrateTgas, n, T, wp, r, nH, nHe4, nHe3, L, z1, z2, phi*factor, nameUpgraded );


	// wait until all threads are finished
	t2.join();
	t4.join();
	t6.join();
	
	cout << "\nProcess completed!\n" << endl;
	
	return(0);

	}
	
	else{
		cout << "only input filename suffix" << endl;
		cout << "usage: " << argv[0] << " -<1 word descriptor>" << endl;
		return 2;
	}
}
