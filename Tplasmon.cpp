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


// T-plasmon flux integrand
double Tintegrand( double w, double m, double ne, double T,
				double wp, double r, double G, double m2g ) {
	//if (isnan(G)) { cout << "AHH!!" << endl; }
	return pow(m*m*r/(pi*R),2)*G*w*sqrt(w*w - m*m)/((exp(w/T)-1)*(pow(m*m - m2g, 2) + pow(w*G,2)));
}


// integral over r
// inputs in eV
// outputs dPhi/dw [eV2]
double Ttrapeze( double w, double m, vector<double> ne, vector<double> T, vector<double> wp,
	 	vector<double> r, vector<double> nH, vector<double> nHe4, vector<double> nHe3,
	 	vector<vector<double>> z1, vector<vector<double>> z2 ) {
		
	int len = r.size();	// get length of vector
	double total = 0;	// initiate value of sum at 0
	// perform integration by looping over r values
	for ( int c = 0; c < len - 1; c++ ) {
	
		if ( w*w < m*m ) { continue; }
		
		double dr = r[c+1] - r[c];	// define trapezium spacing
		
		// select g(w, T) value from matrix
		int indexT1;
		int indexT2;
		int indexX1;
		int indexX2;
		
		for( int i = 1; i < 200; i++ ) {
			if( z1[0][i] < T[c] and z1[0][i+1] > T[c] ) { indexT1 = i; }
			if( z2[0][i] < T[c] and z2[0][i+1] > T[c] ) { indexT2 = i; }
		}
		
		for( int i = 1; i < 500; i++ ) {
			if( (z1[i][0] * T[c]) < w and (z1[i+1][0] * T[c]) > w ) { indexX1 = i; }
			if( (z2[i][0] * T[c]) < w and (z2[i+1][0] * T[c]) > w ) { indexX2 = i; }
		}
		
		double g1 = z1[ indexT1 ][ indexX1 ];
		double g2 = z2[ indexT2 ][ indexX2 ];
				
		double G = GammaOld(w, ne[c], T[c], nH[c], nHe4[c], nHe3[c], g1, g2 );
		
		double m2g = wp[c]*wp[c];
		if( (m>=wp[0]) and (pow(m*m - m2g , 2) < pow(w*G/10 , 2)) ) {	// resonance check
			if ( w < wp[c] ) { continue; }
			total = dr * pow(m2g*r[c]/(pi*R),2) * sqrt(w*w - m2g) / (w*G*(exp(w/T[c])-1));
			//cout << pow(m*m - m2g , 2) << "		" << pow(w*G , 2)/10 << endl;
			if(isnan(sqrt(w*w - m2g))) {cout << "AHH!!" << endl;}
			break;
		}
		else{
		double G2 =  GammaOld(w, ne[c+1], T[c+1], nH[c+1], nHe4[c+1], nHe3[c+1], g1, g2 );
		double m2g2 = wp[c+1]*wp[c+1];
		total += 0.5*dr*( Tintegrand(w, m, ne[c], T[c], wp[c], r[c],G,m2g)
			+ Tintegrand(w, m, ne[c+1], T[c+1], wp[c+1], r[c+1],G2,m2g2) );
		}
	}
	//cout << "alive" << endl;
	return total;
}


// full integration over omega
double Tintegrate( double m, vector<double> ne, vector<double> T, vector<double> wp,
	 	vector<double> r, vector<double> nH, vector<double> nHe4, vector<double> nHe3,
	 	vector<vector<double>> z1, vector<vector<double>> z2,
		bool gas, double L, double wmin ) {
		
	double T0 = 300*K2eV;	// [eV]
	
	double total = 0;	// initiate value of sum at 0
	for ( double w = 1e0; w < 1e4; w*=1.1 ) {
		if ( w <= m ) { continue; }	// only allow when energy greater than mass
		if ( w > m + 1e3 ) { continue; }	// set integral cutoff
		if ( w < wmin ) { continue; }
		double dw = 0.1*w;
		if ( gas == false ) {
			total += 0.5*dw*( ( Pvacuum(w+dw,m,L) * Ttrapeze(w+dw, m, ne, T, wp, r, nH, nHe4, nHe3, z1, z2 ) ) 
				+ ( Pvacuum(w,m,L) * Ttrapeze(w,m, ne, T, wp, r, nH, nHe4, nHe3, z1, z2) ) );
		}
		else if ( gas == true ) {
			double pressure = m*m*m_e*T0/(8*pi*a);		// [eV4]
			total += 0.5*dw*( ( Pgas(w+dw,m,L,pressure) * Ttrapeze(w+dw,m, ne, T, wp, r, nH, nHe4, nHe3, z1, z2) ) 
				+ ( Pgas(w,m,L,pressure) * Ttrapeze(w,m, ne, T, wp, r, nH, nHe4, nHe3, z1, z2) ) );
		}
	}
	return total;
}


// void functions for splitting jobs
void Tplasmon(  vector<double> ne, vector<double> T, vector<double> wp,
	 	vector<double> r, vector<double> nH, vector<double> nHe4, vector<double> nHe3,
	 	vector<vector<double>> z1, vector<vector<double>> z2,
		bool gas, double L, string name ) {
	
	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;
	
	double wmin = 100;	// [eV]

	// set path for writeout
	string path = "data/limits/";
	string ext = ".dat";
	if(gas) { ext = "-gas.dat"; }

	// get for many DP mass values
	for ( double m = 1e-6; m < 1e4; m*=1.1 ) {
			double phi = Tintegrate( m, ne, T, wp, r, nH, nHe4, nHe3, z1, z2, gas, L, wmin );
			chiIAXO.push_back(phi);
			massIAXO.push_back(m) ;
			cout << name << ":	m = " << m << "	phi X-4 = " << phi << endl;
	}

	// write out
	write2D( path + name + ext, massIAXO, chiIAXO );
}



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
	
	bool gas = true;
	

	// babyIAXO
	string nameBaby = "babyIAXO-tPlasmon" + suffix;
	double L1 = 10 / m2eV;	// in eV^-1
	// multithread to run simultaneously
	thread t2( Tplasmon, ne, T, wp, r, nH, nHe4, nHe3, z1, z2, gas, L1, nameBaby );	

	
	//IAXO baseline
	string nameBaseline = "baselineIAXO-tPlasmon" + suffix;
	double L2 = 20 / m2eV;	// in eV^-1
	// multithread to run simultaneously
	thread t4( Tplasmon, ne, T, wp, r, nH, nHe4, nHe3, z1, z2, gas, L2, nameBaseline );


	/// IAXO upgraded
	string nameUpgraded = "upgradedIAXO-tPlasmon" + suffix;
	double L3 = 22 / m2eV;	// in eV^-1
	// multithread to run simultaneously
	thread t6( Tplasmon,ne, T, wp, r, nH, nHe4, nHe3, z1, z2, gas, L3, nameUpgraded );


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
