// Tom O'Shea 2023

// calculation of IAXO sensitivity to solar hidden photons
// header file containing all the bits

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h> 
#include <filesystem>
#include <bits/stdc++.h>

using namespace std;

// constants
double pi = 3.14159265358979323846;	// pi
double m_e = 510998.950;		// m_e [eV]
double a = (1. / 137.);		// alpha
double Ry = 13.605693;			// Rydberg energy [eV]
double aBohr = 5.291772109e-11;	//Bohr radius [m]
double c0 = 299792458;			// speed of light [ms-1]

// solar params
double R_raw = 149.5978707e9;	// mean earth-sun distance [m]
double rSolar_raw = 6.957e8;	// solar radius [m]
double B0 = 3e3;	// radiative zone max B [T] 200;//
double B1 = 50;	// tachocline max B [T] 4;//
double B2 = 4;	// outer region max B [T] 3;//
double r0 = 0.712;	// [R0]
double r1 = 0.732;	// [R0]
double d1 = 0.02;	// [R0]
double r2 = 0.96;	// [R0]
double d2 = 0.035;	// [R0]

// conversion factors
double s2eV = (6.582119569e-16);	// Hz to eV
double J2eV = (1. / 1.602176634e-19);	// Joules to eV (1 / e)
double m2eV = (1.973269804e-7);	// m-1 to eV
double K2eV = (8.617333262e-5);	// Kelvin to eV
double kg2eV = 5.609588604e35;	// from hbar/c2

// converted values
double R = R_raw / m2eV;	// eV-1
double rSolar = rSolar_raw / m2eV;	// eV-1

// utility constants
double wRange = 1e3;	// range of w integral
bool savenquit = false;	// for error catching
//int line = 0;	// for REST writeout


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// READ & WRITE ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// read in datafiles in csv form
vector<double> read( string name ) {

	//cout << "Reading file " << name << "..." << endl;

	// open file defined in argument
	fstream file;
	file.open(name,ios::in);
	
	char delim('\n');	// define delimiter for file parsing
	
	if (file.is_open()){	// checking whether the file is open
		string temp;		// define temporary storage string
		vector<double> row;	// define vector to store input values and return
		
		while(getline(file, temp, delim)){  // read data from file object and put it into string
			double item = stod(temp);	// convert string to double
			row.push_back(item);	// append double to vector
		}
		
	file.close();   // close the file object.
	
	return row;	// returns vector of values
	
	}
	else{ cout << "couldn't find file: " << name << endl ; vector<double> err = {69.} ; return err ; }
}


// read in gaunt factors from matlab matrix files
vector<vector<double>> readGaunt( string name ) {

	//cout << "Reading file " << name << "..." << endl;

	// open file defined in argument
	fstream file;
	file.open(name,ios::in);
	
	char delim(' ');	// define delimiter for file parsing
	
	if (file.is_open()){   // checking whether the file is open
		string temp;	// define temporary storage string
		vector<vector<double>> g2D; // 2D matrix to store rows
		vector<double> row;	// define vector to store input values and return
		vector<double> x;	// define vector to store input values and return
		
		int c = 0;	// define counter for elements
		int line = 0;	// counter for lines
		
		while( !file.eof() ) {  // loop until end of file
		
			getline(file, temp, delim);	// get one row
			
			// check data is not empty
			if( temp == "\n" ) { continue; }
			
			double item = stod(temp);	// convert string to double

			// add a zero on first line
			if( c == 0 ) { row.push_back(0.); }
			
			row.push_back(item);	// append double to vector
			c++;
			
			// when row is full append to 2D vector and wipe row				
			if( row.size() == 201 ) {
			
				g2D.push_back(row);
				row.clear();
				
				line++;
			}

			if( line == 501 ) { break; }

		}
		
	file.close();   // close the file object.
	
	return g2D;	// returns vector of values
	}
	
	else{ cout << "couldn't find file: " << name << endl ; vector<vector<double>> err = {{69.}} ; return err ; }
}


// read column from 2 column datafile
vector<double> loadtxt( string name, int col ) {

	//cout << "Reading file " << name << "..." << endl;

	// open file defined in argument
	fstream file;
	file.open(name,ios::in);
	
	char delim(' ');	// define delimiter for file parsing (tab)
	
	if ( file.is_open() ){   // checking whether the file is open
		string temp;	// define temporary storage string
		vector<double> row1, row2;	// define vector to store input values and return
		vector<string> v;
		
		while( getline(file, temp) ){ v.push_back( temp ); }

		for ( string i : v ) {

			stringstream X(i);
			string temp;
			vector<string> vec;
			while ( getline(X, temp, delim ) ) { vec.push_back(temp); }
			row1.push_back( stod(vec[0]) );
			row2.push_back( stod(vec[1]) );
		}
		
	file.close();   // close the file object.
	
	// choose column
	if ( col == 0 ) { return row1; }
	else if ( col == 1 ) { return row2; }
	else { cout << "put 0 or 1 for columns" << endl; return {69.}; }
	}
	else{ cout << "file " << name << " doesn't exist" << endl; return {69.}; }
}


// write out simple datafile
void write( string name, vector<double> data ) {

	// delete old file if extant
	if ( remove(name.c_str()) == 0 ) {
		cout << "file " << name << " deleted." << endl;
		}
	
	// file pointer
	fstream fout;

	// creates new csv file
	fout.open(name, ios::out | ios::trunc);

	// Read the input from vector
	for ( double item : data ) {
		// Insert the data to file
		fout << item << endl;
	}
	cout << "file " << name << " created succesfully." << endl;
	
	fout.close();
}


// write out 2 column datafile
void write2D( string name, vector<double> data1, vector<double> data2) {

	// delete old file if extant
	if ( remove(name.c_str()) == 0 ) {
		cout << "file " << name << " deleted." << endl;
		}
	
	// create new file
	fstream fout;
	fout.open(name, ios::out | ios::trunc);
	
	// get size of vectors
	int len1 = data1.size();
	int len2 = data2.size();
	
	// check vectors of equal length
	if ( len1 == len2 ) {
		
		// Read the input from vector
		for ( int c = 0; c < len1; c++ ) {
			// Insert the data to file
			fout << data1[c] << "	" << data2[c] << endl;
		}
		
		cout << "file " << name << " created succesfully." << endl;
	}
	
	// if vectors not equal return error
	else { cout << "ERROR - ensure vectors are of equal length" << endl; }
	
	fout.close();
}

// 2D writeout for REST flux
// https://github.com/rest-for-physics/axionlib-data/tree/cae7d4df4ca2f0bef5a2268602cf090fdf4a208d/solarFlux
void writeREST( string name, vector<double> data, int line ){	// data vector needs values of w for constant r

	//int piece = 0;	// for newlining

	// delete old file if extant
	if ( line == 0 ) {
	if ( remove(name.c_str()) == 0 ) {
		cout << "file " << name << " deleted." << endl;
		}
	}

	// file pointer
	fstream fout;

	// creates new csv file
	fout.open(name, ios::out | ios::app);

	// Read the input from vector
	for ( double item : data ) {
		// Insert the data to file
		fout << item << "	";
	}

	// newline at end of w range
	fout << "\b" << endl;	// first remove last tab

	if( line == 0 ){ cout << "file " << name << " created...\n" << endl; }
	
	fout.close();
}


// merge and write out integration output
void merge( string name ) {
	
	// read in data
	vector<double> m1 = read("data/" + name + "-mass-suppressed.dat");
	vector<double> m2 = read("data/" + name + "-mass-resonance.dat");
	vector<double> m3 = read("data/" + name + "-mass-unsuppressed.dat");
	
	vector<double> chi1 = read("data/" + name + "-chi-suppressed.dat");
	vector<double> chi2 = read("data/" + name + "-chi-resonance.dat");
	vector<double> chi3 = read("data/" + name + "-chi-unsuppressed.dat");
	
	// merge
	vector<double> m;
	m.insert( m.end(), m1.begin(), m1.end() );
	m.insert( m.end(), m2.begin(), m2.end() );
	m.insert( m.end(), m3.begin(), m3.end() );
	
	vector<double> chi;
	chi.insert( chi.end(), chi1.begin(), chi1.end() );
	chi.insert( chi.end(), chi2.begin(), chi2.end() );
	chi.insert( chi.end(), chi3.begin(), chi3.end() );
	
	// write out
	string path = "data/limits/";
	string ext = ".dat";
	write2D( path + name + ext, m, chi );

}

// complex sqrt of (a + ib), starting points at sqrts and results stored at pointa, pointb
void csqrt( double a, double b, double* pointa, double* pointb ) {

	double x = sqrt(a);
	double y = sqrt(b);
	int n = 0;

	while( n < 10 ) {
		x = 0.5 * ( x + (a*x + b*y) / ( pow(x,2) + pow(y,2) ) );
		y = 0.5 * ( y + (b*x - a*y) / ( pow(x,2) + pow(y,2) ) );
		n++;
	}

	*pointa = x;
	*pointb = y;
}

// allow writeout of data upon interruption
void interrupt( int sig ) {

	// give warning message
	cout << "Keyboard interrupt" << endl;
	cout << "writing out current data" << endl;

 	savenquit = true;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// RE[PI] , IM[PI] ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// Debye screening scale k_D
double debye( double T, double ne, double np ) {
	double nHepp = 0;
	double nHep = 0;
	return sqrt(4*pi*a*(ne+np+4*(nHep + nHepp))/T);
}


// CDF of Holtsmark distro
// taken from D. G. Hummer (1986)
double QHoltsmark( double B ) {
	double sum = 0;
	for( int n = 0; n < 100; n++ ) {
		double bn = (3/(2*n + 3))*tgamma(4*n/3 + 2)/tgamma(2*n + 2)*pow(-1,n);
		sum += bn*pow(B,2*n);
	}
	return (4/(9*pi))*B*B*B*sum;
}

// cross section for H- photoionisation
// taken from A. Ohmura AND H. Ohmura (1960)
double sigmaHminus( double w ) {	// w in eV
	double k = (c0/(w/s2eV))/aBohr;	// wavenumber in terms of a0
	double gamma = 0.2355883;
	double rho = 2.646;
	return (gamma*k*k*k/((1 - gamma*rho)*pow(gamma*gamma + k*k, 3))) * 6.8475e-18*pow(100/m2eV,2);	// [eV-2]
}

// prob of finding atom in state n, Zn
double probZn( double n, double T, double np ) {
	double Kn = 0;
	if(n>3) { Kn = (16*n*n*(n+(7/6)))/(3*pow(n+1,2)*(n*n + n + 0.5)); }
	else{ Kn = 1; }
	double En = Ry/(n*n);	// [eV]
	double wn = QHoltsmark( (Kn*En*En/(4*a*a))*pow(4*pi*np/3,-2/3) );
	return 2*n*n*wn*exp(En/T);
}

double partFunc( double T, double np ) {
	double sum = 0;
	for( int n = 1; n < 10; n++ ) { sum += probZn(n,T,np); }
	return sum;
}


// free-free Gaunt factor
double GauntFF( double w, double T, double ne, double np ) {
	double kD = debye(T,ne,np);
	double y = kD*sqrt(2*m_e*T);
	double total = 0;
	double dx = 1;
	for( double x1 = 0; x1!=-1; x1+=dx ) {
		double x2 = x1+dx;
		double tot1 = 0;
		double tot2 = 0;
		for( double t = sqrt(x1*x1 + (w/T)) - x1; t <= sqrt(x1*x1 + (w/T)) + x1; t += x1/50 ) {
			if(x1==0) { tot1 = 0; break; }
			tot1 += t*t*t/pow(t*t + y*y, 2);
		}
		for( double t = sqrt(x2*x2 + (w/T)) - x2; t <= sqrt(x2*x2 + (w/T)) + x2; t += x2/50 ) {
			if(x2==0) { tot1 = 0; break; }
			tot2 += t*t*t/pow(t*t + y*y, 2);
		}
		double I1 = 0.5*x1*exp(-x1*x1)*(sqrt(x1*x1+(w/T))/x1)*((1-exp(-2*pi*a*sqrt(m_e/(2*T*(x1*x1+(w/T))))))/(1-exp(-2*pi*a*sqrt(m_e/(2*T*x1*x1)))))*tot1;
		double I2 = 0.5*x2*exp(-x2*x2)*(sqrt(x2*x2+(w/T))/x2)*((1-exp(-2*pi*a*sqrt(m_e/(2*T*(x2*x2+(w/T))))))/(1-exp(-2*pi*a*sqrt(m_e/(2*T*x2*x2)))))*tot2;
		
		if( abs(I1-I2) < 1e-100 ) { break; }
		total += (dx/2)*(I1+I2);
	}
	if(isnan(total)) { return 0; }
	else{ return total; }
}

// free electron Gamma
double GammaFree( double ne ) {
	return 8*pi*a*a*ne/(3*m_e*m_e);
}

// free-free Gamma
double GammaFF( double w, double T, double ne, double np ) {
	double Fff = GauntFF(w,T,ne,np);
	return (64*pi*pi*a*a*a/(3*m_e*m_e*w*w*w))*sqrt(m_e/(2*pi*T))*(1-exp(-w/T))*ne*np*Fff;
}


// bound-free Gamma for H
double GammaBF( double w, double T, double nH0, double np, double nHminus ) {
	double sum = 0;
	double Ztilde = partFunc(T,np);
	for( int n = 0; n < 10; n++ ) {
		double En = Ry/(n*n);	// [eV]
		if( w < En ) { continue; }
		double Zn = probZn(n,T,np)/Ztilde;
		sum += Zn / pow(n,5);
	}
	return (8*pi*m_e*pow(a,5)/(3*sqrt(3)*w*w*w))*(1-exp(-w/T))*nH0*sum + (1-exp(-w/T))*nHminus*sigmaHminus(w);
}


// bound-bound Gamma
double GammaBB( double w, double T, double np, double nH0 ) {
	double sum = 0;
	double Ztilde = partFunc(T,np);
	for( int n = 0; n < 10; n++ ) {
		double Zn = probZn(n,T,np)/Ztilde;
		for( int n1 = 0; n < 10; n++ ) {
			if(n <= n1) { continue; }
			double wr = Ry*(pow(n,-2) - pow(n1,-2));
			double gr = 2*a*wr*w/(3*m_e);
			double fnn = 1;	// ASK!!
			sum += fnn*w*w*gr/( pow(w*w - wr*wr, 2) + (w*w*gr*gr) );
		}
	}
	return (4*pi*a/m_e)*nH0*sum;
}


// bound-bound m^2_gamma
double m2BB( double w, double T, double nH0, double np ) {
	double sum = 0;
	double Ztilde = partFunc(T,np);
	for( int n = 0; n < 10; n++ ) {
		double Zn = probZn(n,T,np)/Ztilde;
		for( int n1 = 0; n < 10; n++ ) {
			if(n <= n1) { continue; }
			double wr = Ry*(pow(n,-2) - pow(n1,-2));
			double gr = 2*a*wr*w/(3*m_e);
			double fnn;
			sum += fnn*w*w*(w*w - wr*wr)/( pow(w*w - wr*wr, 2) + (w*w*gr*gr) );
		}
	}
	return (4*pi*a/m_e)*nH0*sum;
}


// bound-free m^2_gamma
double m2BF( double w, double T, double nH0, double np ) {
	double sum = 0;
	double Ztilde = partFunc(T,np);
	for( int n = 0; n < 10; n++ ) {
	double En = Ry/(n*n);	// [eV]
		double Zn = probZn(n,T,np)/Ztilde;
		sum += Zn*pow(n,-5)*( pow(w/En,2) - log(En*En/abs(En*En - w*w)) );
	}
	return (8*pi*m_e*pow(a,5)/(3*sqrt(3)*w*w))*nH0*sum;
}


// total m^2_gamma
double m2_gamma( double w, double T, double wp, double nH0, double np ) {
	double total = 0;
	total += (wp*wp);
	total += m2BB(w,T,nH0,np);
	total += m2BF(w,T,nH0,np);
	return total;
}

// total Gamma
double Gamma( double w, double wp, double T, double nH0, double ne, double np, double nHminus ) {
	double total = 0;
	total += GammaFree(wp);
	total += GammaFF(w, T, ne, np);
	total += GammaBF(w, T, nH0, np, nHminus);
	total += GammaBB(w,T,np,nH0);
	return total;
}

// Gamma[eV] for 4He in IAXO buffer gas
// taken from Julia's thesis equation (8.5)
// CERN-THESIS-2009-042
double Gamma_IAXO( double w, double pressure ) {	// p[eV4] and w[eV]
	w /= 1000;				// [keV]
	pressure /= (100*m2eV*s2eV*s2eV*kg2eV);		// [mbar]	
	//if( pressure > 1e4 ) { return 1e200; }
	double logw = log10(w);
	double logG = 0.014*pow(logw,6) + 0.166*pow(logw,5)
				+ 0.464*pow(logw,4) + 0.473*pow(logw,3)
				- 0.266*pow(logw,2) - 3.241*logw
				- 0.760 + log10(pressure) - log10(300/1.8);
	return pow(10,logG) * m2eV;		// Gamma[eV]
}


// bound-free m^2_gamma fpr He in IAXO buffer gas
double m2IAXO( double w, double T, double nHe0 ) {
	double En = 4*Ry;	// [eV]
	return (8*pi*m_e*pow(a,5)/(3*sqrt(3)*w*w))*nHe0*(pow(w/En,2) - log(En*En/abs(En*En - w*w)));
}


///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// IAXO CONVERSION //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// IAXO conversion probability for vacuum run
double Pvacuum( double w, double m, double L ) {

	// calculate Delta P
	double dP = w - sqrt(pow(w,2) - pow(m,2));
	if( isnan(dP) ) { return 0; }
	
	// plug into sin2
	double arg = dP * L / 2;
	if( arg > 20 * pi ) { return 2; }	// average
	else{ 
		if( arg < 0.1 ) { return pow(m*m*L/(2*w),2); }
		else{ return 4 * pow( sin(arg) , 2 ); }
	}
}


// conversion probability for gas resonance
double Pgas( double w, double m, double L, double pressure ) {
	
	if( w < m ) { return 0; }
	// define detector parameters for Gamma_t
	double T = 300*K2eV;							// eV
	double G = Gamma_IAXO(w, pressure);				// eV
	double mp2 = 8*pi*a*pressure/(m_e*T);			// eV2
	double dP = sqrt(w*w - mp2) - sqrt(w*w - m*m);

	if( abs(mp2 - m*m) > G*w/10 ) {	// non-resonant conversion
		//cout <<( mp2-(m*m))/mp2 << endl;
		return (m*m*m*m/(pow(mp2 - m*m, 2) + pow(w*G,2))) * (1 + exp(- G*L) - (2 * exp(- G*L / 2) * cos(dP*L)));
		//return  pow(m*m/mp2,2)*(1 + exp(- G*L) - (2 * exp(- G*L / 2) * cos(dP*L)));
	}
	else{
		if( G*L < 1e-3 ) { return pow( m * m * L / (2 * w), 2 ); }
		else{ return pow( m * m / (G * w), 2 ) * (1 + exp(- G*L) - (2 * exp(- G*L / 2))); }
	}
}

// contribution from l-DPs converted in IAXO B-field
/// in limit wp_gas -> m
double PpureL_B( double m, double w, double B, double L, double pressure ) {

	// photon absorption
	double G = Gamma_IAXO(w, pressure);
	
	// set B-field value
	double eB = B * pow(m2eV,2) / s2eV;
	double wB = eB / m_e;

	// conversion prob (book 6, Tue 10/10/2023)
	if( G*L < 1e-3 ) { 
		return pow( m * m / w , 4 ) * pow(wB*L/2, 2) 
				/ ( pow( w*w - m*m , 2 ) + w*w*G*G );
	}
	else{
		return pow( m * m / (w*w) , 4 ) * pow(wB*wB/G, 2)
					* ( 1 + exp(-G*L) - 2*exp(-G*L/2) )
				/ ( pow( w*w - m*m , 2 ) + w*w*G*G );
	}
}

double PpureL_ideal( double m, double w, double L, double pressure ) {
	double G = Gamma_IAXO(w, pressure);
	if( G * L < 1e-6 ) { return pow( m * L / 2 , 2 ); }
	else { return pow( m / G , 2 ) * ( 1 + exp(-2*G*L) - (2 * exp( -G*L)) ); }

}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// SOLAR B FIELDS //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// B field in solar radiative zone [T]
double radiativeZone( double r ) {
	double l = (10 * r0) + 1;
	double K = (1+l) * pow( 1 + pow(l,-1) , l ) * B0;
	double part1 = pow( r/r0 , 2 );
	double part21 = 1 - pow( r/r0 , 2 );
	double part2 = pow( part21, l );
	
	return K * part1 * part2;
}


// B field in tachocline and outer regions [T]
double tachocline( double r, double rmax, double d, double B ) {
	return B * ( 1 - pow( (r - rmax)/d , 2 ) );
}


///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// T-PLASMON /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// POST-DISCOVERY //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// full conversion prob for logL of gas run
double PgasFlux( double w, double m, double mg2IAXO, double L ) {
	
	// define detector parameters for Gamma_t
	double nH, nHe3 = 0;	// only 4He is used
	double g1 = 0;
	double T = 300 * K2eV;	// detector at room temp [eV]
	double pressure = m;		// WRONG - change
	double ne = m_e * pow(m,2) / (4 * pi * a);
	double nHe0 = ne / 2;	// 4He ion density [eV3]
	double G = Gamma_IAXO(w, pressure );

	return (pow(m,4) / ( pow( pow(mg2IAXO,2) - pow(m,2), 2 ) + pow(w*G,2) ))
				*(1 + exp(- G*L) - (2 * exp(- G*L / 2)));
}

/*
// gas integration over omega for logL
double integrateGasFlux( double m, double wpIAXO, vector<double> n, vector<double> T, vector<double> wp,
	 vector<double> r, vector<double> nH, vector<double> nHe4, vector<double> nHe3, double L, vector<vector<double>> z1, vector<vector<double>> z2 ) {
	
	double total = 0;	// initiate value of sum at 0
	double dw = 10.;
	for ( double w = 1e2; w < 2e4 - dw; w+=dw ) {
		if ( w > m + 1e3 ) { continue; }	// set integral cutoff
		if ( w <= m ) { continue; }	// only allow when energy greater than mass
		total += 0.5 * ( ( PgasFlux(w+dw,m,mg2IAXO,L) * Ttrapeze(w+dw,m,ne,T,wp,r,nH,np,nHminus) ) 
			+ ( PgasFlux(w,m,mg2IAXO,L) * Ttrapeze(w,m,ne,T,wp,r,nH,np,nHminus) ) );
	}
	return total;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////// NUCLEAR /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// 5.49 MeV DPs from nuclear fusion
// also used for e+ e- annihilation 511 keV DPs
// use flux from https://arxiv.org/pdf/2305.14420.pdf 
// for now, using values from centre of solar core
void ppchain( vector<double> plasmaFreq, vector<double> temperature, double L, string name ) {
	
	// define constants
	double rate = 1.7e38 * s2eV;			// rate of deuterium fusion [eV]
//	double w = 5.49e6;						// fusion energy [eV]
	double w = 5.11e5;						// electron rest mass [eV]
	double T = temperature[0];
	double wp = plasmaFreq[0];
	double G = compton(w, T, wp);			// compton absorbtion rate [eV]
	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;
	//vector<double> sine;

	// set path for writeout
	string path = "data/limits/";
	string ext = ".dat";

	for ( double m = 1e1; m < w; m*=1.01 ) {
		//double prob = P( w, m, L );		// back-conversion prob
		double prob = 0.5;	// averages out to 1/2 after m ~ 1eV
		double phi = prob * rate * pow(w*w - m*m, 1.5) * pow(w,-3) * pow(m,4) / ( (4*pi*R*R) * ( pow( m*m - wp*wp, 2 ) + pow(w*G,2) ) );
		chiIAXO.push_back( 2* phi );	// times 2 for diphoton emmision in e+ e- mode
		massIAXO.push_back( m );
		//sine.push_back(prob);
	}
	write2D( path + name + ext, massIAXO, chiIAXO );
	//write2D( "data/" + name + "-highEprob.dat", massIAXO, sine );
}


// get energy loss bounds for pp-chain
void ppchainLuminosity( vector<double> plasmaFreq, vector<double> temperature, string name ) {

	// define constants
	double rate = 1.7e38 * s2eV;			// rate of deuterium fusion [eV]
	double w = 5.49e6;						// fusion energy [eV]
//	double w = 5.11e5;						// electron rest mass [eV]
	double T = temperature[0];
	double wp = plasmaFreq[0];
	double G = compton(w, T, wp);			// compton absorbtion rate [eV]
	double L0 = 3.828e26 * J2eV * s2eV;			// solar luminosity [eV2]
	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;

	// set path for writeout
	string path = "data/limits/";
	string ext = "-Eloss.dat";

	for ( double m = 1e1; m < w; m*=1.01 ) {
		double Lt = rate * pow(w*w - m*m, 1.5) * pow(w,-2) * pow(m,4) / ( ( pow( m*m - wp*wp, 2 ) + pow(w*G,2) ) )/L0;
		chiIAXO.push_back( pow(Lt,-0.5) );	// times 2 for diphoton emmision in e+ e- mode
		massIAXO.push_back( m );
	}
	write2D( path + name + ext, massIAXO, chiIAXO );
}

/*
// integrand for nuclear deexcitation DPs
double nuclearIntegrand( double m, double w, double T, double wp, double nI, double r, double J1, double J0 ) {
	double G = compton(w, T, wp);
	return nI*r*r*( (2*J1 + 1)/(2*J0 + 1) )*exp(-w/T)/( pow(wp*wp - m*m, 2) + pow(w*G,2));
}


// integrate over r for deexcitation DPs
double nuclearIntegral( double m, double w, vector<double> r, vector<double> wp, vector<double> T, vector<double> nI, double J1, double J0, double tau ) {

	double intg = 0.;
	for(int c = 0; c < r.size(); c++ ) {
		intg += 0.5 * (r[c+1] - r[c]) * (nuclearIntegrand(m,w,T[c+1],wp[c+1],nI[c+1],r[c+1],J1,J0) + nuclearIntegrand(m,w,T[c],wp[c],nI[c],r[c],J1,J0));
	}
	return pow(m,4) * pow(sqrt(w*w - m*m)/w, 3) * intg / ( tau * R*R );
}


// get flux for deexcitation DPs
void nuclearFlux( vector<double> r, vector<double> wp, vector<double> T, vector<double> nI, double L, 
					double w, double J1, double J0, double tau, string name ) {

	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;

	// set path for writeout
	string path = "data/limits/";
	string ext = ".dat";
	
	double min = *min_element( wp.begin(), wp.end() );
	double max = *max_element( wp.begin(), wp.end() );
	int len = wp.size();

	for( double m = 1e-3; m < min; m *= 1.01 ) {
	//	double prob = P( w, m, L );		// back-conversion prob
		double prob = 0.5;	// averages out to 1/2 after m ~ 1eV
		double phi = prob * nuclearIntegral( m, w, r, wp, T, nI, J1, J0, tau );
		massIAXO.push_back(m);
		chiIAXO.push_back(phi);
	}
	
	for( int c = 0; c < len; c++ ) {
	//	double prob = P( w, m, L );		// back-conversion prob
		double m = wp[len - c - 1];
		double prob = 0.5;	// averages out to 1/2 after m ~ 1eV
		double phi = prob * nuclearIntegral( m, w, r, wp, T, nI, J1, J0, tau );
		massIAXO.push_back(m);
		chiIAXO.push_back(phi);
	}
	
	for( double m = max; m < w; m *= 1.001 ) {
	//	double prob = P( w, m, L );		// back-conversion prob
		double prob = 0.5;	// averages out to 1/2 after m ~ 1eV
		double phi = prob * nuclearIntegral( m, w, r, wp, T, nI, J1, J0, tau );
		massIAXO.push_back(m);
		chiIAXO.push_back(phi);
	}
	
	write2D( path + name + ext, massIAXO, chiIAXO );
}


///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// L-PLASMON MIXING /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


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
	else{
		double wB = eB / m_e;
		double Gt = Gamma(wp, ne, T, nH, nHe4, nHe3, g1, g2);

		double p1 = pow( r/R , 2 ) / (8*pi);
		double p2 = wp * pow(m,4) * pow(wB,2) * pow( pow(wp,2) - pow(m,2) , 0.5 );
		double p3 = exp( wp / T ) - 1;
		double p4 = pow( pow(m,2) - pow(wp,2) , 2 ) + pow(wp*Gt,2);
		
		double item = p1 * p2 / ( p3 * p4 );
		return item;
	}
}


double lMixingResIntegrate( double m, vector<double> ne, vector<double> T, vector<double> wp, vector<double> nH, vector<double> nHe4, vector<double> nHe3, vector<vector<double>> z1, vector<vector<double>> z2, vector<double> r, double L ) {

	double total = 0;	// initiate value of sum at 0
	int len = r.size();
	
	// integrate by trapezium rule over r array
	for ( int c = 0; c < len - 1; c++ ) {
	
		int j = len - c - 2;
	
		if ( wp[j] <= m ) { continue; }	// only allow when energy greater than mass
		//if ( wp[j+1] < 100 ) { continue; }	// only detectable above 0.1 keV
		if ( wp[j+1] > 300 ) { continue; }	// only detectable below 100 eV
		
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
			
		
			double dr = abs( r[j+1] - r[j] );
			double height = 0.5 * ( ( P( wp[j+1], m, L )* lMixingResIntegrand( m, ne[j+1], T[j+1], wp[j+1], r[j+1], nH[j+1], nHe4[j+1], nHe3[j+1], g1, g2 ))
				+ ( P( wp[j], m, L ) * lMixingResIntegrand( m, ne[j], T[j], wp[j], r[j], nH[j], nHe4[j], nHe3[j], g1, g2 ) ) );
			double dA = dr * height;
			
			// only add if real
			if ( isnan(dA) ) { continue; }
			else { total += dA; }
			}
	}
	return total;

}

double lMixingResGasIntegrate( double m, vector<double> ne, vector<double> T, vector<double> wp, vector<double> nH,
								vector<double> nHe4, vector<double> nHe3, vector<vector<double>> z1,
								vector<vector<double>> z2, vector<double> r, double L ) {

	double total = 0;	// initiate value of sum at 0
	int len = r.size();
	
	// integrate by trapezium rule over r array
	for ( int c = 0; c < len - 1; c++ ) {
	
		int j = len - c - 2;
	
		if ( wp[j] <= m ) { continue; }	// only allow when energy greater than mass
		if ( wp[j+1] < 30 ) { continue; }	// only detectable above 0.1 keV
		//if ( wp[j+1] > 300 ) { continue; }	// **EXPT** only detectable below 0.1 keV

		
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
			
		
			double dr = abs( r[j+1] - r[j] );
			double height = 0.5 * ( ( PgasFull( wp[j+1], m, L, z2 ) * lMixingResIntegrand( m, ne[j+1], T[j+1], wp[j+1], r[j+1], nH[j+1], nHe4[j+1], nHe3[j+1], g1, g2 ))
				+ ( PgasFull( wp[j], m, L, z2 ) * lMixingResIntegrand( m, ne[j], T[j], wp[j], r[j], nH[j], nHe4[j], nHe3[j], g1, g2 ) ) );
			double dA = dr * height;
			
			// only add if real
			if ( isnan(dA) ) { continue; }
			else { total += dA; }
			}
	}
	return total;

}


// now to run the integral
void lMixingRes( vector<double> n, vector<double> T, vector<double> wp,	vector<double> nH, vector<double> nHe4,
				vector<double> nHe3, vector<vector<double>> z1, vector<vector<double>> z2, vector<double> r, double L, string name ) {
	
	// implement new interrupt with save
	//signal( SIGINT, interrupt );

	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;

	string path = "data/limits/";
	string ext = "-lMixingRes-flux.dat";
	
	for ( double m = 1e-4; m < 1e3; m*=1.001 ) {

		// check if interrupt
		if( savenquit ){
		// write out
			write2D( path + name + "-INT" + ext, massIAXO, chiIAXO );
			// wait for other processes to save too
			cout << "waiting 10s to let other threads finish" << endl;
			sleep(10);
			exit(SIGINT);
		}

		double entryIAXO = lMixingResIntegrate( m, n, T, wp, nH, nHe4, nHe3, z1, z2, r, L );
		chiIAXO.push_back( entryIAXO );
		massIAXO.push_back( m );
		//cout << name << ":	m = " << m << "	flux X-4 = " << entryIAXO << endl;
	}

	// write out
	write2D( path + name + ext, massIAXO, chiIAXO );
	
}


void lMixingResGas( vector<double> n, vector<double> T, vector<double> wp,vector<double> nH, vector<double> nHe4,
				vector<double> nHe3, vector<vector<double>> z1, vector<vector<double>> z2, vector<double> r, double L, string name ) {
	
	// implement new interrupt with save
	//signal( SIGINT, interrupt );

	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;
	
	string path = "data/limits/";
	string ext = "-lMixingResGas-flux.dat";

	for ( double m = 1e-4; m < 5; m*=1.1 ) {

		// check if interrupt
		if( savenquit ){
		// write out
			write2D( path + name + "-INT" + ext, massIAXO, chiIAXO );
			// wait for other processes to save too
			cout << "waiting 2s to let other threads finish" << endl;
			sleep(2);
			exit(SIGINT);
		}

		double entryIAXO = lMixingResGasIntegrate( m, n, T, wp, nH, nHe4, nHe3, z1, z2, r, L );
		chiIAXO.push_back( entryIAXO );
		massIAXO.push_back( m );
		cout << name << ":	m = " << m << "	flux X-4 = " << entryIAXO << endl;
	}

	// write out
	write2D( path + name + ext, massIAXO, chiIAXO );
	
}


///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// PURE L-PLASMON ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// resonant gas conversion in IAXOdarkphoton-pureL.cpp
double PpureL( double m, double wp, double L ) {

	// define detector parameters for Gamma_t
	double nH, nHe3 = 0;	// only 4He is used
	double T = 300 * K2eV;	// detector at room temp [eV]
	double ne = m_e * pow(m,2) / (4 * pi * a);
	double nHe4 = ne / 2;	// 4He ion density [eV3]
	
	double Gl = GammaLfull( wp, T, ne, nH, nHe4, nHe3, wp, m );
	
	double p1 = pow( m / Gl , 2 );
	double p2 = 1 + exp(-Gl*L) - (2 * exp( -0.5*Gl*L) );
		
	if( Gl * L < 1e-6 ) { return pow( m * L / 2 , 2 ); }
	else { double item = p2 * p1; return item; }

}


double PpureLConversion( double w, double m, double L, vector<vector<double>> z2, double B ) {

	// define detector parameters for Gamma_t
	double nH, nHe3 = 0;	// only 4He is used
	double g1 = 1;
	double T = 300 * K2eV;	// detector at room temp [eV]
	double ne = m_e * pow(m,2) / (4 * pi * a);
	double nHe4 = ne / 2;	// 4He ion density [eV3]
	
	// select Gaunt factor from matrix
	int indexT2 = 60;
	int indexX2;
	for( int i = 1; i < 500; i++ ) {
		if( (z2[i][0] * T) < w and (z2[i+1][0] * T) > w ) { indexX2 = i; break; }
	}
	double g2 = z2[ indexT2 ][ indexX2 ];
	
	// get gammas
	double Gt = Gamma(w, ne, T, nH, nHe4, nHe3, g1, g2);
	double Gl = GammaLfull( w, T, ne, nH, nHe4, nHe3, m, m );

	// set B-field value
	double eB = B * pow(m2eV,2) / s2eV;
	double wB = eB / m_e;

	// full version
	double conv = pow(m,8) * pow(wB,2) * ( 1 + exp(-Gt*L) - 2*exp(-0.5*Gt*L) ) / ( 2 * pow(w,4) * pow(Gt,2) * ( pow( pow(w,2) - pow(m,2) , 2 ) + pow(w*Gl,2) ) );

	// low m limit
	double convLow = pow(m,8) * pow(wB,2) * pow(L,2) / ( 4 * pow(w,4) * ( pow(w,2) + pow(Gl,2) ) );

	if ( Gt*L < 1e-6 ) { return convLow; }
	else { return conv; }
}

// resonant potential crystal conversion approx
double PpureL_crystal( double w, double m, double L ) {

	// refractive indices
	double nx = 1.380;	// 1.544;	// quartz n_o
	double nl = 1.385;	// 1.553;	// quartz n_e
	double nplus = pow(nl,2) + pow(nx,2) - 2;
	double nminus = pow(nl,2) - pow(nx,2);

	double dP = ( sqrt( pow(w,2) - pow(m,2) ) - w * sqrt( 1 - nplus/2 ) );
	double dm = pow(w,2) * nplus;

	double item = ( pow(w,2) * pow(m,6) * pow(nminus,2) / ( 4 * pow(dm,4) ) );// * ( 1 + 8 * pow( cos( 0.5 * dP * L ) , 2 ) );

	//cout << pow(nminus,2) << endl;
	return item;
}




// integrand for pure L gas conversion
double pureLintegrand( double m, double T, double wp, double r ) {

	double p1 = pow( r / (2*R) , 2 ) / pi;
	double p2 = wp * pow(m,2) * sqrt( pow(wp,2) - pow(m,2) );
	double p3 = exp( wp / T ) - 1;
	
	double item = p1 * p2 / p3;
	return item;

}

// calculation for induced current (returns squared value)
double current( double m, double w, double L ) {
	
	double k =  sqrt(pow(w,2) - pow(m,2));
	return pow( 2 * sin(k*L/2) / k , 2);
}

double pureLintegrate( double m, vector<vector<double>> z2, vector<double> T,
						vector<double> wp, vector<double> r, double L, double B ) {

	double total = 0;	// initiate value of sum at 0
	int len = r.size();
	
	// integrate by trapezium rule over r array
	for ( int j = 0; j < len - 1; j++ ) {
	
		if ( wp[j] <= m ) { continue; }	// only allow when energy greater than mass
		//if ( wp[j+1] > 10 ) { continue; }	// only res conversion up to 5eV allowed
		
		else {
			
			if ( wp[j+1] > 300 ) { cout << wp[j+1] << endl; }

			double dr = abs(r[j+1] - r[j]);
			double height = 0.5 * ( (  PpureL( m, wp[j+1], L ) * pureLintegrand( m, T[j+1], wp[j+1], r[j+1] ))
				+ ( PpureL( m, wp[j], L) * pureLintegrand( m, T[j], wp[j], r[j] ) ) );
			//double height = 0.5 * ( ( PpureL_B( m, wp[j+1], B, L, z2 ) * pureLintegrand( m, T[j+1], wp[j+1], r[j+1] ))
			//	+ ( PpureL_B( m, wp[j], B, L, z2) * pureLintegrand( m, T[j], wp[j], r[j] ) ) );
			//double height = 0.5 * ( wp[j+1] * pureLintegrand( m, T[j+1], wp[j+1], r[j+1] ) * current( m, wp[j+1], L ) )
			//	+ ( ( wp[j] * pureLintegrand( m, T[j], wp[j], r[j] ) * current( m, wp[j], L ) ) );
			double dA = dr * height;
			
			// only add if real
			if ( isnan(dA) ) { continue; }
			else { total += dA; }
			}
	}
	return total;
}


// now to run the integral
void pureL( vector<vector<double>> z2, vector<double> T,
			vector<double> wp, vector<double> r, double L, double B, string name ) {

	// implement new interrupt with save
	//signal( SIGINT, interrupt );
	
	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;

	string path = "data/limits/";
	string ext = "-pureL.dat";

	//double P0 = 1.e-12 * J2eV * pow(m2eV,2) * s2eV;	// 1 W m-2 in eV
	//double J0 = 1e-9;		// 1 uA threshold

	for ( double m = 1e-6; m < 1e3; m*=1.1 ) {

		// check if interrupt
		if( savenquit ){
		// write out
			write2D( path + name + "-INT" + ext, massIAXO, chiIAXO );
			// wait for other processes to save too
			cout << "waiting 10s to let other threads finish" << endl;
			sleep(10);
			exit(SIGINT);
		}

		double fluxIAXO = pureLintegrate( m, z2, T, wp, r, L, B );		// E0 / R
		//cout << name << ":	m = " << m << "	flux = " << fluxIAXO << endl;
		chiIAXO.push_back(fluxIAXO);
		massIAXO.push_back(m);
		if( fluxIAXO == 0. ) { break; }
	}

	// write out
	write2D( path + name + ext, massIAXO, chiIAXO );
}


///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// REST DATA /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// T-plasmon flux integrand in REST approx (w >> m)
double integrandREST( double w, double T, double wp, double r ) {

	//double G = Gamma(w, number, T, nH, nHe4, nHe3, g1, g2);
	
	double p1 = w * pow(r,2) * pow(pi * R, -2);
	//double p2 = pow(w,2) * G;
	double p3 = exp(w/T) - 1;
	//double p5 = pow( pow(wp,2) - pow(m,2), 2 ) + pow(w * G, 2);

	double p5, p2 = 1;
	//if ( sel == 0 ) { p5 = pow( pow(wp,2), 2) + pow(w * G, 2); }	// sup
	//else if ( sel == 1 ) { p5 = pow(w * G, 2); }	// res
	//else if ( sel == 2 ) { p5 = 1; }	// unsup (divide by m4)

	double item = p1 / p3;
	
	return item;
}

// suppressed section (approx m << wp, m << w) (missing chi2 m4)
double supREST( double w, double n, double T, double wp, double r, double nH, double nHe4, double nHe3, double g1, double g2 ) {

	double G = Gamma(w, n, T, nH, nHe4, nHe3, g1, g2);
	cout << pow( wp*wp / ( w * G ), 2 ) << endl;
	return pow(r*w/(pi*R*wp*wp), 2) * G / ( exp(w/T) - 1 );
}

// resonant section (approx m ~ wp, m << w) (missing chi2 m4)
double resREST( double w, double n, double T, double wp, double r, double nH, double nHe4, double nHe3, double g1, double g2 ) {

	double G = Gamma(w, n, T, nH, nHe4, nHe3, g1, g2);
	//cout << "r = " << r/rSolar << "	w = " << w << "		" << pow( wp*wp / ( w * G ), 2 ) << endl;

	return pow(r/(pi*R), 2) / G / ( exp(w/T) - 1 );
}

// unsuppresed section (approx m >> wp, m << w) (missing chi2)
double unsupREST( double w, double n, double T, double wp, double r, double nH, double nHe4, double nHe3, double g1, double g2 ) {

	double G = Gamma(w, n, T, nH, nHe4, nHe3, g1, g2);
	return pow(r*w/(pi*R), 2) * G / ( exp(w/T) - 1 );
}

// integral over r
double trapezeREST( double w, vector<double> T, vector<double> wp,
	  vector<double> r, double low ) {

	int len = r.size();	// get length of vector

	double total = 0;	// initiate value of sum at 0
	double high = low + 0.001*rSolar;
	
	// perform integration by looping over r values
	for ( int c = 0; c < len - 1; c++ ) {

		if ( low == 0 ) {
			r2 = r[c+1] - 0.001*rSolar;
			r1 = r[c] - 0.001*rSolar;
		}
		else { r2 = r[c+1]; r1 = r[c]; }

		if ( r1 < low ) { continue; }
		else if ( r1 > high ) { continue; }

		double dr = r2 - r1;	// define trapezium spacing
		double height = 0.5 * ( integrandREST(w, T[c], wp[c], r1 ) 
			+ integrandREST(w, T[c+1], wp[c+1], r2 ) );
		double dA = abs(dr * height);
		
		// only add if real
		if ( isnan(dA) ) { continue; }
		else { total += dA; }
	}
		
	return total;
		
}


// full integration over omega
double integrateREST( vector<double> T, vector<double> wp,
	  vector<double> r, double wlow, double rlow ) {
	
	double total = 0;	// initiate value of sum at 0

	// integrate by trapezium rule over w array
	//double dw = 1e2;
	//for ( double w = 1e2; w < 1e5 - dw; w+=dw ) {
	double whigh = wlow + 100;	// eV
	double dw = 1e0;
	for ( double w = wlow; w <= whigh - dw; w+=dw ) {
	
		//if ( w <= m ) { continue; }	// only allow when energy greater than mass
		
		//else {
		
			double height = 0.5 * ( trapezeREST( w+dw, T, wp, r, rlow ) 
				+  trapezeREST( w, T, wp, r, rlow ) );
			double dA = abs(dw * height);
			
			// only add if real
			if ( isnan(dA) ) { continue; }
			else { total += dA; }
		//}
	}
	return total;
}


// m in eV
//void fluxREST( int select ) {
void fluxREST() {

	int line = 0;	// for REST writeout

	// read csv files to vectors
	vector<double> r = read("data/r.dat");	// sun radial distance [eV-1]
	vector<double> rFrac = read("data/rFrac.dat");	// sun radial distance as fraction
	vector<double> T = read("data/T.dat");	// solar temperature [eV]
	vector<double> wp = read("data/wp.dat");	// plasma frequency [eV]

	int len = r.size();	// get length of vector
	double i = 1;	// to get solar radius fraction

	// initialise vector etc
	vector<double> flux;
	string name = "RESTstuff/HiddenPhotonFlux_OShea_202309(2).dat";
	//string name = "data/flux_m" + to_string(intm) + "_X" + to_string(intchi) + "-again.dat";

	for ( double rlow = 0; rlow < 0.999*rSolar; rlow += 0.001*rSolar ) {
		for ( double wlow = 1; wlow < 20e3; wlow += 100 ) {
	
//			double avg = integrateREST2(n, T, wp, r, nH,
//						 nHe4, nHe3, z1, z2, wlow, rlow, select ) / 0.1;	// eV3 keV-1
			double avg = integrateREST( T, wp, r, wlow, rlow ) / 0.1;	// eV3 keV-1			
			// convert eV2 to cm-2 s-1 keV-1
			avg *= pow( m2eV, -2 ) * 1e-4 / s2eV;
			flux.push_back(avg);
			//if ( avg < 1e-35 ) { flux.push_back(0); }
			//else { flux.push_back( avg ); }
		}
		writeREST( name, flux, line );
		line++;
		flux.clear();
		// note here:	\33[A moves up a line
		// 				\33[2K deletes line
		// 				\r goes to start of line
		//cout << "\33[A\33[2K\r" << line/1000 << "\% complete..." << endl;
		cout << "Line " << line << " out of 1000" << endl;

	}
	cout << "	creation of file " << name << " completed!" << endl;
}


// integral over r
double trapezemx( double w, double m, vector<double> T, vector<double> wp,
	  vector<double> r, vector<double> ne, vector<double> nH, vector<double> nHe4,
	  vector<double> nHe3, vector<vector<double>> z1, vector<vector<double>> z2, double low ) {

	int len = r.size();	// get length of vector

	double total = 0;	// initiate value of sum at 0
	double high = low + 0.001*rSolar;
	
	// perform integration by looping over r values
	for ( int c = 0; c < len - 1; c++ ) {

		if ( r[c] < low ) { continue; }
		else if ( r[c] > high ) { continue; }
		
		// select g(w,T) value from matrix
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

		double dr = r[c+1] - r[c];	// define trapezium spacing
		double height = 0.5 * ( integrand(w, m, ne[c], T[c], wp[c], r[c], nH[c], nHe4[c], nHe3[c], g1, g2) 
			+ integrand(w, m, ne[c+1], T[c+1], wp[c+1], r[c+1], nH[c+1], nHe4[c+1], nHe3[c+1], g1, g2) );
		double dA = abs(dr * height);
		
		// only add if real
		if ( isnan(dA) ) { continue; }
		else { total += dA; }
	}
		
	return total;
		
}


// full integration over omega
double integratemx( double m, vector<double> T, vector<double> wp,
	  vector<double> r, vector<double> ne, vector<double> nH, vector<double> nHe4,
	  vector<double> nHe3, vector<vector<double>> z1, vector<vector<double>> z2,
	  double wlow, double rlow ) {
	
	double total = 0;	// initiate value of sum at 0

	// integrate by trapezium rule over w array
	//double dw = 1e2;
	//for ( double w = 1e2; w < 1e5 - dw; w+=dw ) {
	double whigh = wlow + 100;	// eV
	double dw = 1e0;
	for ( double w = wlow; w < whigh - dw; w+=dw ) {
	
		//if ( w <= m ) { continue; }	// only allow when energy greater than mass
		
		//else {
		
			double height = 0.5 * ( trapezemx( w+dw, m, T, wp, r, ne, nH, nHe4, nHe3, z1, z2, rlow ) 
				+  trapezemx( w, m, T, wp, r, ne, nH, nHe4, nHe3, z1, z2, rlow ) );
			double dA = abs(dw * height);
			
			// only add if real
			if ( isnan(dA) ) { continue; }
			else { total += dA; }
		//}
	}
	return total;
}


void fluxmx( double m ) {

	int line = 0;	// for REST writeout

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

	int len = r.size();	// get length of vector
	double i = 1;	// to get solar radius fraction

	// initialise vector etc
	vector<double> flux;
	int intm = (int)log10(m);
	string name = "data/flux_m" + to_string(intm) + "-newbins.dat";

	for ( double rlow = 0; rlow < 0.999*rSolar; rlow += 0.001*rSolar ) {
		for ( double wlow = 1; wlow < 20e3; wlow += 100 ) {
	
//			double avg = integrateREST2(n, T, wp, r, nH,
//						 nHe4, nHe3, z1, z2, wlow, rlow, select ) / 0.1;	// eV3 keV-1
			double avg = integratemx( m, T, wp, r, ne, nH, nHe4, nHe3, z1, z2, wlow, rlow ) / 0.1;	// eV3 keV-1		
			// convert eV2 to cm-2 s-1 keV-1
			avg *= pow( m2eV, -2 ) * 1e-4 / s2eV;
			//flux.push_back(avg);
			flux.push_back( avg );
		}
		writeREST( name, flux, line );
		line++;
		flux.clear();
		// note here:	\33[A moves up a line
		// 				\33[2K deletes line
		// 				\r goes to start of line
		cout << "\33[A\33[2K\r" << (line)/10 << "\% complete..." << endl;
		//cout << "\33[2K\r" << line << "\% complete...";

	}
	cout << "	creation of file " << name << " completed!" << endl;
}


// average Gamma over r bin
double trapezeG( double w, vector<double> n, vector<double> T,
	  vector<double> r, vector<double> nH, vector<double> nHe4, vector<double> nHe3,
	  vector<vector<double>> z1, vector<vector<double>> z2, double low ) {

	int len = r.size();	// get length of vector

	double total = 0;	// initiate value of sum at 0
	double high = low + 0.001*rSolar;
	
	// perform integration by looping over r values
	for ( int c = 0; c < len - 1; c++ ) {

		if ( r[c] < low ) { continue; }
		else if ( r[c] > high ) { continue; }

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

		double dr = r[c+1] - r[c];	// define trapezium spacing
		double height = 0.5 * ( Gamma( w, n[c], T[c], nH[c], nHe4[c], nHe3[c], g1, g2 )
						+ Gamma( w, n[c+1], T[c+1], nH[c+1], nHe4[c+1], nHe3[c+1], g1, g2 ) );
		double dA = abs(dr * height);
		
		// only add if real
		if ( isnan(dA) ) { cout << "nan" << endl; continue; }
		else { total += dA; }
	}
	//cout << low << "	" << w << "	" << total << endl;
	return total;	// return average over bin
}

// average Gamma over energy bin
double integrateG( vector<double> n, vector<double> T,
	  vector<double> r, vector<double> nH, vector<double> nHe4, vector<double> nHe3,
	  vector<vector<double>> z1, vector<vector<double>> z2, double wlow, double rlow ) {
	
	double total = 0;	// initiate value of sum at 0

	// integrate by trapezium rule over w
	double whigh = wlow + 100;	// eV
	double dw = 10;
	for ( double w = wlow; w < whigh - dw; w+=dw ) {

		double height = 0.5 * ( w * trapezeG( w+dw, n, T, r, nH, nHe4, nHe3, z1, z2, rlow ) 
			+ ( (w+dw) * trapezeG( w, n, T, r, nH, nHe4, nHe3, z1, z2, rlow ) ) );
		double dA = abs(dw * height);
		
		// only add if real
		if ( isnan(dA) ) { cout << "nan" << endl; continue; }
		else { total += dA; }
	}
	return total;	// return average over energy bin
}


// get average omega*Gamma
void wGout() {

	int line = 0;	// for REST writeout

	vector<double> r = read("data/rFrac.dat");	// sun radial distance [eV-1]
	vector<double> T = read("data/T.dat");	// solar temperature [eV]
	vector<double> n = read("data/ne.dat");	// electron number density [eV3]
	vector<double> nH = read("data/nH.dat");	// H ion density [eV3]
	vector<double> nHe4 = read("data/nHe4.dat");	// He4 ion density [eV3]
	vector<double> nHe3 = read("data/nHe3.dat");	// He3 ion density [eV3]
	
	// get gaunt factors
	vector<vector<double>> z1 = readGaunt("data/Z1.dat");	// gaunt factors for Z=1
	vector<vector<double>> z2 = readGaunt("data/Z2.dat");	// gaunt factors for Z=2

	// convert Gaunt factor Theta to T in eV
	for( int i = 1; i < 201; i++ ) { z1[0][i] = z1[0][i] * m_e; }
	for( int i = 1; i < 201; i++ ) { z2[0][i] = z2[0][i] * m_e; }

	// initialise vector etc
	vector<double> wG;
	string name = "data/wG-avg3.dat";

	for ( double rlow = 0; rlow < 0.999*rSolar; rlow += 0.001*rSolar ) {
		for ( double wlow = 1; wlow < 20e3; wlow += 100 ) {

			double avg = integrateG( n, T, r, nH,
						 nHe4, nHe3, z1, z2, wlow, rlow ) / 100 / 0.001*rSolar;	// eV			
			wG.push_back( avg );
			//cout << avg << endl;
		}
		writeREST( name, wG, line );
		line++;
		wG.clear();
		// note here:	\33[A moves up a line
		// 				\33[2K deletes line
		// 				\r goes to start of line
		cout << "\33[A\33[2K\r" << line/10 << "\% complete..." << endl;
	}
	cout << "	creation of file " << name << " completed!" << endl;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// SPECTRA ETC ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// function to output E spectrum for T-plasmon
void spectrum( double m ) {
	
	//set chi
	double chi = 1;

	vector<double> r = read("data/r.dat");	// sun radial distance [eV-1]
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
	
	// initialise phi vector
	vector<double> phi;
	vector<double> E;
	
	for ( double w = 2e2; w < 2e6; w *=1.001 ) {

		if ( w <= m ) { continue; }	// only for m < wp
		
		// compute dphi/dE
		double entry = pow(chi,2) * trapeze( w, m, n, T, wp, r, nH, nHe4, nHe3, z1, z2 );
		E.push_back(w);
		phi.push_back(entry);
	}
	

	int mint = (int)log10(m);
	//string name = "data/Espectrum" + to_string( mint ) + ".dat";
	string name = "data/Espectrum-max.dat";
	
	write2D( name , E, phi );
	cout << "length of E vector: " << E.size() << endl;
	cout << "length of phi vector: " << phi.size() << "\n" << endl;
	
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
			double entry = lMixingResIntegrand( m, ne[j], T[j], wp[j], r[j], nH[j], nHe4[j], nHe3[j], g1, g2 )/pow(m,4);
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


void spectrumPureL( double m ) {
	
	vector<double> r = read("data/r.dat");	// sun radial distance [eV-1]
	vector<double> rFrac = read("data/rFrac.dat");	// sun radial distance [eV-1]
	vector<double> T = read("data/T.dat");	// solar temperature [eV]
	vector<double> n = read("data/ne.dat");	// electron number density [eV3]
	vector<double> wp = read("data/wp.dat");	// plasma frequency [eV]
	vector<double> nH = read("data/nH.dat");	// H ion number density [eV3]
	vector<double> nHe4 = read("data/nHe4.dat");	// He4 ion number density [eV3]
	vector<double> nHe3 = read("data/nHe3.dat");	// He3 ion number density [eV3]
	
	// initialise phi vector
	vector<double> phi;
	vector<double> E;
	vector<double> R;
	
	// get phi values for each w = wp for resonance
	int len = wp.size();
	for ( int c = 0; c < len; c++ ) {
	
		int j = len - c - 2;
	
		if ( wp[j] <= m ) { continue; }	// only allow when energy greater than mass
		if ( wp[j] > 1000 ) { continue; }
		
		else {
			double entry = pureLintegrand( m, T[j], wp[j], r[j] ) ;
			E.push_back(wp[j]);
			phi.push_back(entry);
		}
	}
	
	// write to file
	int mint = (int)log10(m);	// for labelling filename
	string name = "data/Espectrum-pureL" + to_string( mint ) + ".dat";
	
	write2D( name , E, phi );
		
	cout << "length of E vector: " << E.size() << endl;
	cout << "length of phi vector: " << phi.size() << "\n" << endl;

}


// function to output B-fields
void Bfields() {
	
	vector<double> r = read("data/rFrac.dat");	// sun radial distance [eV-1]
	vector<double> T = read("data/T.dat");	// solar temperature [eV]
	vector<double> n = read("data/ne.dat");	// electron number density [eV3]
	vector<double> wp = read("data/wp.dat");	// plasma frequency [eV]
	
	// initialise phi vector
	vector<double> B;
	vector<double> E;
	vector<double> R;
	
	int lenWp = wp.size();
	
	// get phi values for each w = wp for resonance
	for ( int j = 0; j < lenWp - 1; j++ ) {
	
		double rFrac = r[j];

				// impoved model
				if ( rFrac <= r0 ) {
					double entry = radiativeZone( rFrac );
					E.push_back(wp[j]);
					R.push_back(rFrac);
					B.push_back(entry);
					}
				else if ( rFrac > (r1 - d1) and rFrac < (r1 + d1) ) {
					double entry = tachocline( rFrac, r1, d1, B1 );
					E.push_back(wp[j]);
					R.push_back(rFrac);
					B.push_back(entry);
					}
				else if ( rFrac > (r2 - d2) and rFrac < (r2 + d2) ) {
					double entry = tachocline( rFrac, r2, d2, B2 );
					E.push_back(wp[j]);
					R.push_back(rFrac);
					B.push_back(entry);
					}
		}
	
	// write to file
	string name = "data/Bfields.dat";
	string nameR = "data/Bfields-R.dat";
	
	write2D( name , E, B );
	write2D( nameR , R, B );

}

// calculate photon mfp from gamma
void mfp(){

	vector<double> r = read("data/rFrac.dat");	// sun radial distance [eV-1]
	vector<double> T = read("data/T.dat");	// solar temperature [eV]
	vector<double> ne = read("data/ne.dat");	// electron number density [eV3]
	vector<double> nH = read("data/nH.dat");	// H ion density [eV3]
	vector<double> nHe4 = read("data/nHe4.dat");	// He4 ion density [eV3]
	vector<double> nHe3 = read("data/nHe3.dat");	// He3 ion density [eV3]
	
	// get gaunt factors
	vector<vector<double>> z1 = readGaunt("data/Z1.dat");	// gaunt factors for Z=1
	vector<vector<double>> z2 = readGaunt("data/Z2.dat");	// gaunt factors for Z=2

	// convert Gaunt factor Theta to T in eV
	for( int i = 1; i < 201; i++ ) { z1[0][i] = z1[0][i] * m_e; }
	for( int i = 1; i < 201; i++ ) { z2[0][i] = z2[0][i] * m_e; }
	
	// initialise vectors
	vector<double> mfp, GammaT;
	double w = 1e1;	// set at 1keV for now
	
	int len = r.size();
	
	// get phi values for each w = wp for resonance
	for ( int j = 0; j < len; j++ ) {
		
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

		double G = Gamma( w, ne[j], T[j], nH[j], nHe4[j], nHe3[j], g1, g2 );	// [eV]
		mfp.push_back( m2eV / G );	// [m]
		GammaT.push_back(G);
	}

	// write out
	string name1 = "data/mfp.dat";
	string name2 = "data/GammaT.dat";
	string name = "data/mfp-m.dat";

	write(name1, mfp);

	//write2D( name1, r, mfp );
	//write2D( name2, r, GammaT );
}


// get gamma array for resonance width
void gammaOut(){

	vector<double> r = read("data/rFrac.dat");	// sun radial distance [eV-1]
	vector<double> T = read("data/T.dat");	// solar temperature [eV]
	vector<double> ne = read("data/ne.dat");	// electron number density [eV3]
	vector<double> nH = read("data/nH.dat");	// H ion density [eV3]
	vector<double> nHe4 = read("data/nHe4.dat");	// He4 ion density [eV3]
	vector<double> nHe3 = read("data/nHe3.dat");	// He3 ion density [eV3]
	
	// get gaunt factors
	vector<vector<double>> z1 = readGaunt("data/Z1.dat");	// gaunt factors for Z=1
	vector<vector<double>> z2 = readGaunt("data/Z2.dat");	// gaunt factors for Z=2

	// convert Gaunt factor Theta to T in eV
	for( int i = 1; i < 201; i++ ) { z1[0][i] = z1[0][i] * m_e; }
	for( int i = 1; i < 201; i++ ) { z2[0][i] = z2[0][i] * m_e; }
	
	// initialise vectors
	vector<double> mfp, GammaT, wGammaT;

	int len = r.size();
	int line = 0;
	int rval = 0;
	
	// names for writeout
	string name1 = "data/mfp3.dat";
	string name2 = "data/GammaT3.dat";
	string name3 = "data/wGammaT3.dat";

	// get phi values for each w = wp for resonance
	for ( int j = 0; j < len; j++ ) { 
	
	if ( 0.001 * rval > r[j] ) { continue; }
	else{		
	//mfp.push_back(r[j]);
	//GammaT.push_back(r[j]);
	//wGammaT.push_back(r[j]);

	for ( double w = 50; w <= 20e3; w += 100 ) {

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

		double G = Gamma( w, ne[j], T[j], nH[j], nHe4[j], nHe3[j], g1, g2 );	// eV
		mfp.push_back( m2eV * 1e2 / G );	// cm
		GammaT.push_back(G);	// eV
		wGammaT.push_back(w*G);	// eV2

	}

	writeREST( name1, mfp, line );
	writeREST( name2, GammaT, line );
	writeREST( name3, wGammaT, line );
	line++;
	rval++;
	// note here:	\33[A moves up a line
	// 				\33[2K deletes line
	// 				\r goes to start of line
	cout << "\33[A\33[2K\r" << ( rval +1 ) / 10 << "\% complete..." << endl;
//	cout << "\33[A\33[2K\r" << ( 100 * line / len ) << "\% complete..." << endl;
	mfp.clear();
	GammaT.clear();
	wGammaT.clear();
	}
	}
	
	cout << "adding zeros..." << endl << endl;

	while ( rval < 1000 ) {
		vector<double> temp;
		for ( int c = 0; c < 200; c++ ) { temp.push_back(0); }
		writeREST( name1, temp, line );
		writeREST( name2, temp, line );
		writeREST( name3, temp, line );
		line++;
		rval++;
		cout << "\33[A\33[2K\r" << ( rval +1 ) / 10 << "\% complete..." << endl;
	}

	cout << "creation of files:\n	" << name1 << "\n	" << name2 << "\n	" << name3 << "\ncompleted!" << endl;
}
*/
// the end :)
