// Tom O'Shea 2023

// script to generate a function of dark photon mixing parameter against dark photon mass
// to give a theoretical upper limit on these parameters that IAXO could achieve

// original t-plasmon integral for IAXO vacuum & gas runs

#include "darkphoton.h"	// home of the real code
#include <thread>

using namespace std;


// T plasmon absorbtion length
// from Redondo's HP Atlas https://arxiv.org/abs/1501.07292
double GammaT( double w, double number, double T, double nH, double nHe4, double nHe3, double g1, double g2 ) {

	double p1 = 64 * pow(pi,2) * pow(a,3);
	double p2 = 3 * pow(m_e,2) * pow(w,3);
	double p3 = m_e * pow(number,2) / (2*pi*T);
	double p4 = 1 - exp(- w / T);
	double p5 = 8 * pi * pow(a,2) * number / (3 * pow(m_e,2) );		// Thompson
	
	// sum of ion densities
	double ions = (nH * g1) + g2 * ( (4 * nHe4) + (4 * nHe3) );

	double item = p1 * pow(p2, -1) * pow(p3, 0.5) * p4 * ions;		// bremsstrahlung
	
	item += p5;
		
	return item;
}


// integrand for pure L gas conversion
double lMixingResIntegrand( double m, double ne, double T, double wp, double r, double nH, double nHe4, double nHe3, double g1, double g2 ) {

	// set values of solar B-field
	double eB;
	double rFrac = r / rSolar;

	// impoved model
	if ( rFrac <= r0 ) {
		eB = radiativeZone( rFrac ) * pow(m2eV,2) / s2eV;
		}
	else if ( rFrac > (r1 - d1) and rFrac < (r1 + d1)) {
		eB = tachocline( rFrac, r1, d1, B1 ) * pow(m2eV,2) / s2eV;
		}
	else if ( rFrac > (r2 - d2) and rFrac < (r2 + d2) ) {
		eB = tachocline( rFrac, r2, d2, B2 ) * pow(m2eV,2) / s2eV;
		}
	else { eB = 0; }
	
	if ( eB == 0 ) { return 0; }
	if ( m >= wp ) { return 0; }
	else{
		double wB = eB / m_e;
		double Gt = GammaT(wp, ne, T, nH, nHe4, nHe3, g1, g2);

		double p1 = pow( r/R , 2 ) / (8*pi);
		//double p2 = wp * pow(m,4) * pow(wB,2) * pow( pow(wp,2) - pow(m,2) , 0.5 );
		double p2 = pow(m,4) * pow(wB,4) * pow( pow(wp,2) - pow(m,2) , 0.5 ) / wp;
		double p3 = exp( wp / T ) - 1;
		//double p4 = pow( pow(m,2) - pow(wp,2) , 2 ) + pow(wp*Gt,2);
		double p4 = (pow(wp,4) + pow(wp*Gt,2))*pow(wp*Gt,2);
		
		double item = p1 * p2 / ( p3 * p4 );
		if(isnan(item)){cout<<"AHH!!"<<endl;}
		return item * pow(wp/m,4);
	}
}


double lMixingResIntegrate( double m, vector<double> ne, vector<double> T, vector<double> wp,
							vector<double> nH, vector<double> nHe4, vector<double> nHe3,
							vector<vector<double>> z1, vector<vector<double>> z2,
							vector<double> r, double L, double pressure, double wmin ) {

	double total = 0;	// initiate value of sum at 0
	int len = r.size();
	
	// integrate by trapezium rule over r array
	for ( int c = 0; c < len - 1; c++ ) {
	
		int j = len - c - 2;
	
		if ( wp[j] <= m ) { continue; }	// only allow when energy greater than mass
		if ( wp[j] < wmin ) { continue; }	// cutoff
				
		double w = wp[j];
	
		// select g(w, T) value from matrix
		int indexT1;
		int indexT2;
		int indexX1;
		int indexX2;
		
		for( int i = 1; i < 200; i++ ) {
			if( z1[0][i] < T[j] and z1[0][i+1] > T[j] ) { indexT1 = i; }
			if( z2[0][i] < T[j] and z2[0][i+1] > T[j] ) { indexT2 = i; }
		}
		
		for( int i = 1; i < 500; i++ ) {
			if( (z1[i][0] * T[j]) < w and (z1[i+1][0] * T[j]) > w ) { indexX1 = i; }
			if( (z2[i][0] * T[j]) < w and (z2[i+1][0] * T[j]) > w ) { indexX2 = i; }
		}
		
		double g1 = z1[ indexT1 ][ indexX1 ];
		double g2 = z2[ indexT2 ][ indexX2 ];
		
	
		double dr = abs( r[j+1] - r[j] );
				
		if( pressure == 0 ) {
		total += dr * 0.5 * ( ( Pvacuum(wp[j+1],m,L)* lMixingResIntegrand( m, ne[j+1], T[j+1], wp[j+1], r[j+1], nH[j+1], nHe4[j+1], nHe3[j+1], g1, g2 ))
			+ ( Pvacuum(wp[j],m,L) * lMixingResIntegrand( m, ne[j], T[j], wp[j], r[j], nH[j], nHe4[j], nHe3[j], g1, g2 ) ) );
		}
		else {
			total += dr * 0.5 * ( ( Pgas(wp[j+1],m,L,pressure) * lMixingResIntegrand(m,ne[j+1],T[j+1],wp[j+1],r[j+1],nH[j+1],nHe4[j+1],nHe3[j+1],g1,g2))
				+ ( Pgas(wp[j],m,L,pressure) * lMixingResIntegrand(m,ne[j],T[j],wp[j],r[j],nH[j],nHe4[j],nHe3[j],g1,g2) ) );
		}
	}
	return total;
}


// function to output E spectrum for L-mixing at resonance
void spectrumResL( double m ) {
	
	vector<double> r = read("data/r.dat");	// sun radial distance [eV-1]
	vector<double> rFrac = read("data/rFrac.dat");	// sun radial distance [eV-1]
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
	
	// initialise phi vector
	vector<double> phi;
	vector<double> E;
	vector<double> R;
	
	int lenWp = wp.size();
	
	// get phi values for each w = wp for resonance
	for ( int j = 0; j < lenWp - 1; j++ ) {

		if ( wp[j] <= m ) { continue; }	// only for m < wp
		
		else {
		
			double w = wp[j];
		
			// select g(w, T) value from matrix
			int indexT1;
			int indexT2;
			int indexX1;
			int indexX2;
			
			for( int i = 1; i < 200; i++ ) {
				if( z1[0][i] < T[j] and z1[0][i+1] > T[j] ) { indexT1 = i; }
				if( z2[0][i] < T[j] and z2[0][i+1] > T[j] ) { indexT2 = i; }
			}
			
			for( int i = 1; i < 500; i++ ) {
				if( (z1[i][0] * T[j]) < w and (z1[i+1][0] * T[j]) > w ) { indexX1 = i; }
				if( (z2[i][0] * T[j]) < w and (z2[i+1][0] * T[j]) > w ) { indexX2 = i; }
			}
			
			double g1 = z1[ indexT1 ][ indexX1 ];
			double g2 = z2[ indexT2 ][ indexX2 ];
		
			double dr = r[j+1] - r[j];
			double entry = lMixingResIntegrand( m, ne[j], T[j], wp[j], r[j], nH[j], nHe4[j], nHe3[j], g1, g2 )
										/(pow(100*m2eV,2)*s2eV*pow(m,4));	// cm-2 s-1 eV-1
			E.push_back(wp[j]);
			phi.push_back(entry);
			R.push_back(rFrac[j]);
		}
	}
	
	// write to file
	int mint = (int)log10(m);	// for labelling filename
	string name = "data/Espectrum-lMixing" + to_string( mint ) + ".dat";
	string nameR = "data/Espectrum-lMixing-R" + to_string( mint ) + ".dat";
	
	write2D( name , E, phi );
	write2D( nameR , R, phi );
		
	cout << "length of E vector: " << E.size() << endl;
	cout << "length of phi vector: " << phi.size() << "\n" << endl;

}

// now to run the integral
void lMixingRes( vector<double> ne, vector<double> T, vector<double> wp, vector<double> nH, vector<double> nHe4,
				vector<double> nHe3, vector<vector<double>> z1, vector<vector<double>> z2, vector<double> r,
				double L, string name ) {

	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;
	double pressure = 0;
	double wmin = 1;	// eV
	double T0 = 300*K2eV;

	// set path for writeout
	string path = "data/limits/";
	string ext = ".dat";

	// get for many DP mass values
	for ( double m = 1e-6; m < 1e4; m*=1.1 ) {
			//pressure = m*m*m_e*T0/(8*pi*a);		// [eV4]
			//pressure = 0.1*0.1*m_e*T/(8*pi*a);				// 0.1eV
			double phi = lMixingResIntegrate( m, ne, T, wp,nH, nHe4, nHe3,z1, z2,r, L, pressure, wmin);
			chiIAXO.push_back(phi);
			massIAXO.push_back(m) ;
			cout << name << ":	m = " << m << "	phi X-4 = " << phi << endl;
	}
	if(pressure!=0) { ext = "-gas.dat"; }
	// write out
	write2D( path + name + ext, massIAXO, chiIAXO );
}


// argc counts number of command line inputs (including the script itself)
// argv is a vector containing cmd args separated by space (including script itself)
int main( int argc, char** argv ) {

	// FOR SPECTRA
	//for( double m = 1e-6; m <= 1e4; m *= 10 ) {	spectrumResL(m) ; }

	// FOR FLUX
	
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
	vector<double> nHe4 = read("data/nHe4.dat");	// H ion density [eV3]
	vector<double> nHe3 = read("data/nHe3.dat");	// H ion density [eV3]
	
	// get gaunt factors
	vector<vector<double>> z1 = readGaunt("data/Z1.dat");	// gaunt factors for Z=1
	vector<vector<double>> z2 = readGaunt("data/Z2.dat");	// gaunt factors for Z=2
	
	// convert Gaunt factor Theta to T in eV
	for( int i = 1; i < 201; i++ ) { z1[0][i] = z1[0][i] * m_e; }
	for( int i = 1; i < 201; i++ ) { z2[0][i] = z2[0][i] * m_e; }
	

	// babyIAXO
	string nameBaby = "babyIAXO-lMixing" + suffix;
	double L1 = 10 / m2eV;	// in eV^-1
	// multithread to run simultaneously
	thread t2(lMixingRes,ne, T, wp,nH, nHe4,nHe3, z1, z2, r,L1,nameBaby);	

	
	//IAXO baseline
	string nameBaseline = "baselineIAXO-lMixing" + suffix;
	double L2 = 20 / m2eV;	// in eV^-1
	// multithread to run simultaneously
	thread t4(lMixingRes,ne, T, wp,nH, nHe4,nHe3, z1, z2, r,L2, nameBaseline );


	/// IAXO upgraded
	string nameUpgraded = "upgradedIAXO-lMixing" + suffix;
	double L3 = 22 / m2eV;	// in eV^-1
	/// multithread to run simultaneously
	thread t6(lMixingRes,ne, T, wp,nH, nHe4,nHe3, z1, z2, r, L3, nameUpgraded );


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
