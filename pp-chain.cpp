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

	// read csv files to vectors
	vector<double> T = read("data/T.dat");	// solar temperature [eV]
	vector<double> wp = read("data/wp.dat");	// plasma frequency [eV]

	// babyIAXO
	string nameBaby = "babyIAXO-pp" + suffix;
	double L1 = 10 / m2eV;	// in eV^-1
	thread t2( ppchain, wp, T, L1, nameBaby );	

	//IAXO baseline
	double L2 = 20 / m2eV;
	string nameBaseline = "baselineIAXO-pp" + suffix;
	thread t4( ppchain, wp, T, L2, nameBaseline );	

	// IAXO upgraded
	string nameUpgraded = "upgradedIAXO-pp" + suffix;
	double L3 = 22 / m2eV;	// in eV^-1
	thread t6( ppchain, wp, T, L3, nameUpgraded );	

	// wait until all threads are finished
	t2.join();
	t4.join();
	t6.join();
	
	cout << "\nProcess completed!\n" << endl;
	
	return(0);

	}
	
	else{
		cout << "only input filename suffix" << endl;
		return 2;
	}
}
