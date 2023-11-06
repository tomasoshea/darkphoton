// Tom O'Shea 2023

// get detector Gamma to check

#include "darkphoton.h"
using namespace std;

// plot G vs p(m) for a given w
int main() {
	//double pressure = 13.48e3*(100*m2eV*s2eV*s2eV*kg2eV);	// [eV4]
	double T = 300*K2eV;									// [eV]
	double w = 5e3; 										// [eV]
	double L = 20/m2eV;										// [eV-1]
	vector<double> probVac,probGas,mass;
	for( double m = 1e-4; m < 1e4; m*=1.001 ) {				// [eV]
		double pressure = m*m*m_e*T/(8*pi*a);				// [eV4]
		mass.push_back(m);									// [keV]
		probVac.push_back(Pvacuum(w,m,L));
		probGas.push_back(Pgas(w,m,L,pressure));
	}
	write2D("data/probVac-5keV.dat", mass, probVac);
	write2D("data/probGas-5keV.dat", mass, probGas);
return 0;
}
