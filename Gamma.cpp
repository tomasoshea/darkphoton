// Tom O'Shea 2023

// get detector Gamma to check

#include "darkphoton.h"
using namespace std;

/*
// plot G vs w for given p(m)
int main() {
	//double pressure = 13.48e3*(100*m2eV*s2eV*s2eV*kg2eV);	// [eV4]
	double T = 300*K2eV;									// [eV]
	double m = 1;											// [eV]
	double pressure = m*m*m_e*T/(8*pi*a);					// [eV4]
	vector<double> Gamma,energy;
	for( double w = 1e3; w < 10e3; w+=0.1e3 ) {				// [eV]
		energy.push_back(w/1e3);							// [keV]
		Gamma.push_back(Gamma_IAXO(w, pressure)/m2eV);		// [m-1]
	}
	write2D("data/GammaVsEnergy-1eV.dat", energy, Gamma);
return 0;
}*/

// plot G vs p(m) for a given w
int main() {
	//double pressure = 13.48e3*(100*m2eV*s2eV*s2eV*kg2eV);	// [eV4]
	double T = 300*K2eV;									// [eV]
	double w = 30;											// [eV]
	vector<double> Gamma,mass;
	for( double m = 1e-4; m < 1e4; m*=1.1 ) {				// [eV]
		double pressure = m*m*m_e*T/(8*pi*a);					// [eV4]
		mass.push_back(m);							// [keV]
		Gamma.push_back(Gamma_IAXO(w, pressure)/m2eV);		// [m-1]
	}
	write2D("data/GammaVsMass-30eV.dat", mass, Gamma);
return 0;
}
