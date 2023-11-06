// Tom O'Shea 2023

// script to generate a function of dark photon mixing parameter against dark photon mass
// to give a theoretical upper limit on these parameters that IAXO could achieve

// original t-plasmon integral for IAXO vacuum & gas runs

#include "darkphoton.h"	// home of the real code
#include <thread>

using namespace std;


// T plasmon absorbtion length
double GammaOld( double w, double number, double T, double nH, double nHe4, double nHe3, double g1, double g2 ) {

	double p1 = 64 * pow(pi,2) * pow(a,3);
	double p2 = 3 * pow(m_e,2) * pow(w,3);
	double p3 = m_e * pow(number,2) / (2*pi*T);
	double p4 = 1 - exp(- w / T);
	double p5 = 8 * pi * pow(a,2) * number / (3 * pow(m_e,2) );		// Thompson
	
	// sum of ion densities
	double ions = (nH * g1) + g2 * ( (4 * nHe4) + (4 * nHe3) );

	double item = p1 * pow(p2, -1) * pow(p3, 0.5) * p4 * ions;		// bremsstrahlung
	
	item += p5;
	
	if (isnan(item)) { cout << "AHH!!" << endl; }
	
	return item;
}

// argc counts number of command line inputs (including the script itself)
// argv is a vector containing cmd args separated by space (including script itself)
int main( int argc, char** argv ) {

	vector<vector<double>> z2 = readGaunt("data/Z2.dat");	// gaunt factors for He	
	// convert Gaunt factor Theta to T in eV
	for( int i = 1; i < 201; i++ ) { z2[0][i] = z2[0][i] * m_e; }

	// detector params
	double T = 300*K2eV;	// room temp [eV]
	double pressure = 0.08 * (100*m2eV*s2eV*s2eV*kg2eV);	// mbar to eV4
	//double pressure = 13.43 * (100*m2eV*s2eV*s2eV*kg2eV);	// mbar to eV4
	double g1, nH, nHe3 = 0;
	double m = sqrt(8*pi*a*pressure/(m_e*T));
	double ne = m_e * pow(m,2) / (4 * pi * a);
	double nHe4 = ne / 2;	// 4He ion density [eV3]
	
	// select Gaunt factor from matrix
	int indexT2;	
	for( int i = 1; i < 200; i++ ) {
		if( z2[0][i] < T and z2[0][i+1] > T ) { indexT2 = i; }
	}

	vector<double> Gamma,energy;
	for( double w = 1e3; w < 10e3; w+=0.1e3 ) {			// [eV]
		int indexX2;
		for( int i = 1; i < 500; i++ ) {
		if( (z2[i][0] * T) < w and (z2[i+1][0] * T) > w ) { indexX2 = i; }
		}
		double g2 = z2[ indexT2 ][ indexX2 ];
		energy.push_back(w*1e-3);										// [keV]
		Gamma.push_back(GammaOld( w, ne, T, nH, nHe4, nHe3, g1, g2)/m2eV);	// [m-1]
	}
	//write2D("data/GammaIAXOold13.dat", energy, Gamma);
	write2D("data/GammaIAXOold08.dat", energy, Gamma);

	cout << "\nProcess completed!\n" << endl;
	
	return(0);
}
