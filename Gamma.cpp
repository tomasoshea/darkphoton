// Tom O'Shea 2023

// get detector Gamma to check

#include "darkphoton.h"
using namespace std;

int main() {
	double pressure = 13.48;						// [mbar]
	vector<double> Gamma,energy;
	for( double w = 1; w < 10; w+=0.1 ) {			// [keV]
		energy.push_back(w);
		Gamma.push_back(Gamma_IAXO(w, pressure));	// [m-1]
	}
	write2D("data/GammaIAXO13.dat", energy, Gamma);
return 0;
}
