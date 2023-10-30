// Tom O'Shea 2023

// script to generate a function of dark photon mixing parameter against dark photon mass
// to give a theoretical upper limit on these parameters that IAXO could achieve

// original t-plasmon integral for IAXO vacuum & gas runs

#include "darkphoton.h"	// home of the real code
#include <thread>

using namespace std;


// T-plasmon flux integrand
double Tintegrand( double w, double m, double ne, double T,
				double wp, double r, double G, double m2g ) {
	return pow(m*m*r/(pi*R),2)*G*w*sqrt(w*w - m*m)/((exp(w/T)-1)*(pow(m*m - m2g, 2) + pow(w*G,2)));
}


// integral over r
// inputs in eV
// outputs dPhi/dw [eV2]
double Ttrapeze( double w, double m, vector<double> ne, vector<double> T, vector<double> wp,
	 	vector<double> r, vector<double> nH0, vector<double> np, vector<double> nHminus ) {
		
	int len = r.size();	// get length of vector
	double total = 0;	// initiate value of sum at 0
	// perform integration by looping over r values
	for ( int c = 0; c < len - 1; c++ ) {
		double dr = r[c+1] - r[c];	// define trapezium spacing
		double G = Gamma( w, wp[c], T[c], 0, ne[c], np[c], 0 );
		double m2g = m2_gamma(w, T[c], wp[c], 0, np[c]);
		if( abs(m*m - m2g) < abs(w*G) ) {	// resonance check
			total = dr * pow(m2g*r[c]/(pi*R),2) * sqrt(w*w - m2g) / (w*G*(exp(w/T[c])-1));
			cout << "resonance" << endl;
			break;
		}
		else{
		double G2 = Gamma( w, wp[c+1], T[c+1], 0, ne[c+1], np[c+1], 0 );
		double m2g2 = m2_gamma(w, T[c+1], wp[c+1], 0, np[c+1]);
		total += 0.5*dr*( Tintegrand(w, m, ne[c], T[c], wp[c], r[c],G,m2g)
			+ Tintegrand(w, m, ne[c+1], T[c+1], wp[c+1], r[c+1],G2,m2g2) );
		}
	}
	cout << "alive" << endl;
	return total;
}


// full integration over omega
double Tintegrate( double m, vector<double> ne, vector<double> T, vector<double> wp,
	 	vector<double> r, vector<double> nH, vector<double> np, vector<double> nHminus,
		bool gas, double L ) {
	
	double total = 0;	// initiate value of sum at 0
	double dw = 1e2;
	for ( double w = 1e2; w < 1e5 - dw; w+=dw ) {
		if ( w <= m ) { continue; }	// only allow when energy greater than mass
		if ( w > m + 1e3 ) { continue; }	// set integral cutoff
		if ( gas == false ) {
			total += 0.5*dw*( ( Pvacuum(w+dw,m,L) * Ttrapeze(w+dw,m,ne,T,wp,r,nH,np,nHminus) ) 
				+ ( Pvacuum(w,m,L) * Ttrapeze(w,m,ne,T,wp,r,nH,np,nHminus) ) );
		}
		else if ( gas == true ) {
			total += 0.5*dw*( ( Pgas(w+dw,m,L) * Ttrapeze(w+dw,m,ne,T,wp,r,nH,np,nHminus) ) 
				+ ( Pgas(w,m,L) * Ttrapeze(w,m,ne,T,wp,r,nH,np,nHminus) ) );
		}
	}
	return total;
}


// void functions for splitting jobs
void Tplasmon( vector<double> ne, vector<double> T, vector<double> wp,
	 	vector<double> r, vector<double> nH, vector<double> np, vector<double> nHminus,
		bool gas, double L, string name ) {
	
	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;

	// set path for writeout
	string path = "data/limits/";
	string ext = ".dat";
	if(gas) { ext = "-gas.dat"; }

	// get for many DP mass values
	for ( double m = 1e-6; m < 1e5; m*=1.1 ) {
			double phi = Tintegrate(m,ne,T,wp,r,nH,np,nHminus,gas,L);
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
	vector<double> np = read("data/nH.dat");	// H ion density [eV3]

	vector<double> nHminus = {0};	// FOR NOW - change
	vector<double> nH = {0};
	
	// get gaunt factors
	vector<vector<double>> z1 = readGaunt("data/Z1.dat");	// gaunt factors for Z=1
	vector<vector<double>> z2 = readGaunt("data/Z2.dat");	// gaunt factors for Z=2
	
	// convert Gaunt factor Theta to T in eV
	for( int i = 1; i < 201; i++ ) { z1[0][i] = z1[0][i] * m_e; }
	for( int i = 1; i < 201; i++ ) { z2[0][i] = z2[0][i] * m_e; }
	

	// babyIAXO
	string nameBaby = "babyIAXO-tPlasmon" + suffix;
	double L1 = 10 / m2eV;	// in eV^-1
	// multithread to run simultaneously
	thread t2( Tplasmon, ne, T, wp, r, nH, np, nHminus, false, L1, nameBaby );	

	
	//IAXO baseline
	string nameBaseline = "baselineIAXO-tPlasmon" + suffix;
	double L2 = 20 / m2eV;	// in eV^-1
	// multithread to run simultaneously
	thread t4( Tplasmon, ne, T, wp, r, nH, np, nHminus, false, L2, nameBaseline );


	/// IAXO upgraded
	string nameUpgraded = "upgradedIAXO-tPlasmon" + suffix;
	double L3 = 22 / m2eV;	// in eV^-1
	// multithread to run simultaneously
	thread t6( Tplasmon, ne, T, wp, r, nH, np, nHminus, false, L3, nameUpgraded );


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
