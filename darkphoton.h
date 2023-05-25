// Tom O'Shea 2022

// script to generate a function of dark photon mixing parameter against dark photon mass
// to give a theoretical upper limit on these parameters that IAXO could achieve

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

// solar params
double R_raw = 149.5978707e9;	// mean earth-sun distance [m]
double rSolar_raw = 6.957e8;	// solar radius [m]
double B0 = 3e3;	// radiative zone max B [T]
double B1 = 50;	// tachocline max B [T]
double B2 = 3;	// outer region max B [T]
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
	
	if (file.is_open()){   // checking whether the file is open
		string temp;	// define temporary storage string
		vector<double> row;	// define vector to store input values and return
		
		while(getline(file, temp, delim)){  // read data from file object and put it into string
			double item = stod(temp);	// convert string to double
			row.push_back(item);	// append double to vector
		}
		
	file.close();   // close the file object.
	
	return row;	// returns vector of values
	
	}
	else{ cout << "some error i guess..." << endl ; vector<double> err = {69.} ; return err ; }
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
	
	else{ cout << "some error i guess..." << endl ; vector<vector<double>> err = {{69.}} ; return err ; }
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
	fout << endl;

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
///////////////////////////////////////// GENERAL FUNCS ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// T plasmon production rate
double Gamma( double w, double number, double T, double nH, double nHe4, double nHe3, double g1, double g2 ) {

	double p1 = 16 * pow(pi,2) * pow(a,3);
	double p2 = 3 * pow(m_e,2) * pow(w,3);
	double p3 = 2 * pi * m_e * pow(number,2) / (3 * T);
	double p4 = 1 - exp(- w / T);
	double p5 = 8 * pi * pow(a,2) * number / (3 * pow(m_e,2) );
	
	// sum of ion densities
	double ions = (nH * g1) + g2 * ( (4 * nHe4) + (4 * nHe3) );

	double item = p1 * pow(p2, -1) * pow(p3, 0.5) * p4 * ions;
	
	item += p5;

	return item;
}


// L plasmon production rate - valid around resonance (w = wp)
double GammaLfull( double w, double T, double ne, double nH, double nHe4, double nHe3, double wp, double m ) {
	//if( w > 300 ) { cout << w << endl; }
	double p1 = 64 * pow(pi,2) * pow(a,3) * ne;
	double p2 = 3 * pow( 2 * pi * T , 0.5 ) * pow(m_e,1.5) * pow(w,3);
	double F = cyl_bessel_k(0, w/2) * sinh(w/2);
	double sum = nH + ( 4 * (nHe4 + nHe3) );
	
	double item = ( p1 * sum * F / p2 );	// + ( p3 / p4 );
	return item;
}


// A_l => S_t mixing rate gamma_lt^4
double mixing( double w, double wp, double eB ) {

	double item = pow( (wp * eB) / (m_e * w) , 4 ) / 8;
	return item;
}


// IAXO conversion probability for vacuum run
double P( double w, double m, double L ) {

	// calculate Delta P
	double dP = w - pow( pow(w,2) - pow(m,2) , 0.5 );
	
	// plug into sin2
	double arg = dP * L / 2;
	
	double item = 4 * pow( sin(arg) , 2 );
	
	return item;
}


// conversion probability for gas resonance in low m limit
double Pgas( double w, double m, double L ) {

	double item = pow( m * m * L / (2 * w), 2 );
	
	return item;
}


// conversion probability for gas resonance
double PgasFull( double w, double m, double L, vector<vector<double>> z2 ) {
	
	// define detector parameters for Gamma_t
	double nH, nHe3 = 0;	// only 4He is used
	double g1 = 0;
	double T = 25 * K2eV;	// detector at room temp [eV]
	double ne = m_e * pow(m,2) / (4 * pi * a);
	double nHe4 = ne / 2;	// 4He ion density [eV3]
	
	// select Gaunt factor from matrix
	int indexT2;
	int indexX2;
	for( int i = 1; i < 200; i++ ) {
		if( (z2[0][i] < T) and (z2[0][i+1] > T) ) { indexT2 = i; break; }
	}
	if( indexT2 == 1 ) { cout << "bad gaunt!!" << endl;}
	else if( indexT2 == 199 ) { cout << "bad gaunt!!" << endl;}


	for( int i = 1; i < 500; i++ ) {
		if( (z2[i][0] * T) < w and (z2[i+1][0] * T) > w ) { indexX2 = i; break; }
	}
	if( indexX2 == 1 ) { cout << "bad gaunt!!" << endl;}
	else if( indexX2 == 499 ) { cout << "bad gaunt!!" << endl;}


	double g2 = z2[ indexT2 ][ indexX2 ];
	
	// calculate detector Gamma
	double G = Gamma(w, ne, T, nH, nHe4, nHe3, g1, g2);

	// calculate conversion prob
	double p1 = pow( m * m / (G * w), 2 );
	//double p2 = 1 + exp(- G*L) - (2 * exp(- G*L / 2));
	double p2 = 1 + exp(- G*L) - (2 * exp(- G*L / 2));
	
	double item = p1 * p2;
	
	return item;
}


// B field in solar radiative zone [T]
double radiativeZone( double r ) {

	double l = (10 * r0) + 1;
	double K = (1+l) * pow( 1 + pow(l,-1) , l ) * B0;
	
	double part1 = pow( r/r0 , 2 );
	double part21 = 1 - pow( r/r0 , 2 );
	double part2 = pow( part21, l );
	
	double a = K * part1 * part2;
	
	return a;	
}


// B field in tachocline and outer regions [T]
double tachocline( double r, double rmax, double d, double B ) {

	double a = B * ( 1 - pow( (r - rmax)/d , 2 ) );
	return a;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// T-PLASMON /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// T-plasmon flux integrand
double integrand( double w, double m, double number, double T, double wp, double r, double nH, double nHe4, double nHe3, double g1, double g2 ) {

	double G = Gamma(w, number, T, nH, nHe4, nHe3, g1, g2);
	
	double p1 = pow(r,2) * pow(pi * R, -2);
	double p2 = w * pow( pow(w,2) - pow(m,2) , 0.5 );
	double p3 = exp(w/T) - 1;
	double p4 = pow(m,4) * G;
	double p5 = pow( pow(wp,2) - pow(m,2), 2) + pow(w * G, 2);

	double item = p1 * p2 * p4 / (p3 * p5);
	
	return item;
}


// integral over r
double trapeze( double w, double m, vector<double> n, vector<double> T, vector<double> wp,
	 vector<double> r, vector<double> nH, vector<double> nHe4, vector<double> nHe3, vector<vector<double>> z1, vector<vector<double>> z2 ) {

	int len = r.size();	// get length of vector

	double total = 0;	// initiate value of sum at 0
	
	// perform integration by looping over r values
	for ( int c = 0; c < len - 1; c++ ) {

		//if ( r[c] > 0.1 * rSolar ) { continue; }
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
		double height = 0.5 * ( integrand(w, m, n[c], T[c], wp[c], r[c], nH[c], nHe4[c], nHe3[c], g1, g2) 
			+ integrand(w, m, n[c+1], T[c+1], wp[c+1], r[c+1], nH[c+1], nHe4[c+1], nHe3[c+1], g1, g2) );
		double dA = abs(dr * height);
		
		// only add if real
		if ( isnan(dA) ) { continue; }
		else { total += dA; }
	}
		
	return total;
		
}


// full integration over omega
double integrate( double m, vector<double> n, vector<double> T, vector<double> wp,
	 vector<double> r, vector<double> nH, vector<double> nHe4, vector<double> nHe3, double L, vector<vector<double>> z1, vector<vector<double>> z2 ) {
	
	double total = 0;	// initiate value of sum at 0

	// integrate by trapezium rule over w array
	//double dw = 1e2;
	//for ( double w = 1e2; w < 1e5 - dw; w+=dw ) {
	double dw = 1e1;
	for ( double w = 1e2; w < 1e4 - dw; w+=dw ) {
	
		if ( w <= m ) { continue; }	// only allow when energy greater than mass
		else if ( w > m + wRange ) { continue; }
		
		else {

			double height = 0.5 * ( ( P( w+dw, m, L ) * trapeze( w+dw, m, n, T, wp, r, nH, nHe4, nHe3, z1, z2 ) ) 
				+ ( P( w, m, L ) * trapeze( w, m, n, T, wp, r, nH, nHe4, nHe3, z1, z2 ) ) );
			//double height = 0.5 * ( ( trapeze( w+dw, m, n, T, wp, r, nH, nHe4, nHe3, z1, z2 ) ) 
			//	+ ( trapeze( w, m, n, T, wp, r, nH, nHe4, nHe3, z1, z2 ) ) );
			double dA = abs(dw * height);
			
			// only add if real
			if ( isnan(dA) ) { continue; }
			else { total += dA; }
		}
	}
	return total;
}


// gas integration over omega
double integrateGas( double m, vector<double> n, vector<double> T, vector<double> wp,
	 vector<double> r, vector<double> nH, vector<double> nHe4, vector<double> nHe3, double L, vector<vector<double>> z1, vector<vector<double>> z2 ) {
	
	double total = 0;	// initiate value of sum at 0

	// integrate by trapezium rule over w array
	//double dw = 1e2;
	//for ( double w = 1e2; w < 1e5 - dw; w+=dw ) {
	//double dw = 100;
	//for ( double w = 100; w < 1e5 - dw; w+=dw ) {
	double dw = 1e1;
	for ( double w = 1e2; w < 1e4 - dw; w+=dw ) {
	
		if ( w > m + 1e3 ) { continue; }	// set integral cutoff
		if ( w <= m ) { continue; }	// only allow when energy greater than mass
		else if ( w > m + wRange ) { continue; }
		
		else {
		
			double height = 0.5 * ( ( PgasFull( w+dw, m, L, z2 ) * trapeze( w+dw, m, n, T, wp, r, nH, nHe4, nHe3, z1, z2 ) ) 
				+ ( PgasFull( w, m, L, z2 ) * trapeze( w, m, n, T, wp, r, nH, nHe4, nHe3, z1, z2 ) ) );
			double dA = abs(dw * height);
			
			// only add if real
			if ( isnan(dA) ) { continue; }
			else { total += dA; }
		}
	}
	return total;
}


// void functions for splitting jobs
void integrateT( vector<double> n, vector<double> T, vector<double> wp,
		vector<double> r, vector<double> nH, vector<double> nHe4, vector<double> nHe3, double L,
		vector<vector<double>> z1, vector<vector<double>> z2, string name ) {
	
	// implement new interrupt with save
	//signal( SIGINT, interrupt );

	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;

	// set path for writeout
	string path = "data/limits/";
	string ext = ".dat";
	
	double min = *min_element( wp.begin(), wp.end() );
	double max = *max_element( wp.begin(), wp.end() );
	
	// suppressed section
	for ( double m = 1e-4; m < min; m*=1.1 ) {

		// check if interrupt
		if( savenquit ){
		// write out
			write2D( path + name + "-INT" + ext, massIAXO, chiIAXO );
			// wait for other processes to save too
			cout << "waiting 10s to let other threads finish" << endl;
			sleep(10);
			exit(SIGINT);
		}

		else{
			double phi = integrate( m, n, T, wp, r, nH, nHe4, nHe3, L, z1, z2 );
			chiIAXO.push_back( phi );
			massIAXO.push_back( m );
			cout << name << ":	m = " << m << "	phi X-4 = " << phi << endl;
			//cout << " flux: " << entryIAXO * sqrt(chi4IAXO) * pow(m2eV,-2) / s2eV << endl;
		}
	}
	

	// resonant sector
	int len = wp.size();
	for ( int i = 0; i < len; i = i+2 ) {

		// check if interrupt
		if( savenquit ){
		// write out
			write2D( path + name + "-INT" + ext, massIAXO, chiIAXO );
			// wait for other processes to save too
			cout << "waiting 10s to let other threads finish" << endl;
			sleep(10);
			exit(SIGINT);
		}
	
		int j = len - i - 1;
		double m = wp[j];
		//if( m > 1 ){ break; }
		
		double phi = integrate( m, n, T, wp, r, nH, nHe4, nHe3, L, z1, z2 );
		chiIAXO.push_back( phi );
		massIAXO.push_back( m );
		cout << name << ":	m = " << m << "	phi X-4 = " << phi <<endl;
		//cout << " flux: " << entryIAXO * sqrt(chi4IAXO) * pow(m2eV,-2) / s2eV << endl;
	}

		// unsuppressed section
	for ( double m = max; m < 1e4; m*=1.1 ) {

		// check if interrupt
		if( savenquit ){
		// write out
			write2D( path + name + "-INT" + ext, massIAXO, chiIAXO );
			// wait for other processes to save too
			cout << "waiting 10s to let other threads finish" << endl;
			sleep(10);
			exit(SIGINT);
		}

		else{
		double phi = integrate( m, n, T, wp, r, nH, nHe4, nHe3, L, z1, z2 );
		chiIAXO.push_back( phi );
		massIAXO.push_back( m );
		cout << name << ":	m = " << m << "	phi X-2 = " << phi << endl;
		}
	}

	// write out
	write2D( path + name + ext, massIAXO, chiIAXO );
}


void integrateTgas( vector<double> n, vector<double> T, vector<double> wp,
		vector<double> r, vector<double> nH, vector<double> nHe4, vector<double> nHe3, double L,
		vector<vector<double>> z1, vector<vector<double>> z2, string name ) {
	
	// implement new interrupt with save
	//signal( SIGINT, interrupt );

	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;

	// set path for writeout
	string path = "data/limits/";
	string ext = "-gas.dat";
	
	double min = *min_element( wp.begin(), wp.end() );
	
	// suppressed section
	for ( double m = 1e-3; m < min; m*=1.1 ) {

		// check if interrupt
		if( savenquit ){
		// write out
			write2D( path + name + "-INT" + ext, massIAXO, chiIAXO );
			// wait for other processes to save too
			cout << "waiting 2s to let other threads finish" << endl;
			sleep(2);
			exit(SIGINT);
		}

		double phi = integrateGas( m, n, T, wp, r, nH, nHe4, nHe3, L, z1, z2 );
		chiIAXO.push_back( phi );
		massIAXO.push_back( m );
		cout << name << ":	m = " << m << "	phi X-2 = " << phi <<endl;
	}
	
	// resonant sector
	int len = wp.size();
	for ( int i = 0; i < len; i++ ) {

		// check if interrupt
		if( savenquit ){
		// write out
			write2D( path + name + "-INT" + ext, massIAXO, chiIAXO );
			// wait for other processes to save too
			cout << "waiting 2s to let other threads finish" << endl;
			sleep(2);
			exit(SIGINT);
		}
	
		int j = len - i - 1;
		double m = wp[j];
		if( m > 5 ){ break; }

		double phi = integrateGas( m, n, T, wp, r, nH, nHe4, nHe3, L, z1, z2 );
		chiIAXO.push_back( phi );
		massIAXO.push_back( m );
		cout<< name << ":	m = " << m << "	phi X-2 = " << phi <<endl;
	}

		// write out
	write2D( path + name + ext, massIAXO, chiIAXO );
	
}


///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// L-PLASMON MIXING /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// integrand for pure L gas conversion
double lMixingResIntegrand( double m, double ne, double T, double wp, double r, double rFrac, double nH, double nHe4, double nHe3, double g1, double g2 ) {

	// set values of solar B-field
	double eB;

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
	
	double wB = eB / m_e;
	double Gt = Gamma(wp, ne, T, nH, nHe4, nHe3, g1, g2);

	double p1 = pow( r/R , 2 ) / (8*pi);
	double p2 = wp * pow(m,4) * pow(wB,2) * pow( pow(wp,2) - pow(m,2) , 0.5 );
	double p3 = exp( wp / T ) - 1;
	double p4 = pow( pow(m,2) - pow(wp,2) , 2 ) + pow(wp*Gt,2);
	
	double item = p1 * p2 / ( p3 * p4 );
	return item;

}


double lMixingResIntegrate( double m, vector<double> ne, vector<double> T, vector<double> wp, vector<double> rFrac, vector<double> nH, vector<double> nHe4, vector<double> nHe3, vector<vector<double>> z1, vector<vector<double>> z2, vector<double> r, double L ) {

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
			double height = 0.5 * ( ( P( wp[j+1], m, L )* lMixingResIntegrand( m, ne[j+1], T[j+1], wp[j+1], r[j+1], rFrac[j+1], nH[j+1], nHe4[j+1], nHe3[j+1], g1, g2 ))
				+ ( P( wp[j], m, L ) * lMixingResIntegrand( m, ne[j], T[j], wp[j], r[j], rFrac[j], nH[j], nHe4[j], nHe3[j], g1, g2 ) ) );
			double dA = dr * height;
			
			// only add if real
			if ( isnan(dA) ) { continue; }
			else { total += dA; }
			}
	}
	return total;

}

double lMixingResGasIntegrate( double m, vector<double> ne, vector<double> T, vector<double> wp, vector<double> rFrac, vector<double> nH, vector<double> nHe4, vector<double> nHe3, vector<vector<double>> z1, vector<vector<double>> z2, vector<double> r, double L ) {

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
			double height = 0.5 * ( ( PgasFull( wp[j+1], m, L, z2 ) * lMixingResIntegrand( m, ne[j+1], T[j+1], wp[j+1], r[j+1], rFrac[j+1], nH[j+1], nHe4[j+1], nHe3[j+1], g1, g2 ))
				+ ( PgasFull( wp[j], m, L, z2 ) * lMixingResIntegrand( m, ne[j], T[j], wp[j], r[j], rFrac[j], nH[j], nHe4[j], nHe3[j], g1, g2 ) ) );
			double dA = dr * height;
			
			// only add if real
			if ( isnan(dA) ) { continue; }
			else { total += dA; }
			}
	}
	return total;

}


// now to run the integral
void lMixingRes( vector<double> n, vector<double> T, vector<double> wp, vector<double> rFrac, vector<double> nH, vector<double> nHe4, vector<double> nHe3, vector<vector<double>> z1, vector<vector<double>> z2, vector<double> r, double L, double phi, string name ) {
	
	// implement new interrupt with save
	signal( SIGINT, interrupt );

	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;

	string path = "data/limits/";
	string ext = "-lMixingRes.dat";
	
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

		double entryIAXO = lMixingResIntegrate( m, n, T, wp, rFrac, nH, nHe4, nHe3, z1, z2, r, L );
		double chi4IAXO = phi / entryIAXO;
		chiIAXO.push_back( pow( chi4IAXO, 0.25 ) );
		massIAXO.push_back( m );
		cout << name << ":	m = " << m << "	chi = " << pow( chi4IAXO, 0.25 ) << endl;
	}

	// write out
	write2D( path + name + ext, massIAXO, chiIAXO );
	
}


void lMixingResGas( vector<double> n, vector<double> T, vector<double> wp, vector<double> rFrac, vector<double> nH, vector<double> nHe4, vector<double> nHe3, vector<vector<double>> z1, vector<vector<double>> z2, vector<double> r, double L, double phi, string name ) {
	
	// implement new interrupt with save
	signal( SIGINT, interrupt );

	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;
	
	string path = "data/limits/";
	string ext = "-lMixingResGas.dat";

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

		double entryIAXO = lMixingResGasIntegrate( m, n, T, wp, rFrac, nH, nHe4, nHe3, z1, z2, r, L );
		double chi4IAXO = phi / entryIAXO;
		chiIAXO.push_back( pow( chi4IAXO, 0.25 ) );
		massIAXO.push_back( m );
		cout << name << ":	m = " << m << "	chi = " << pow( chi4IAXO, 0.25 ) << endl;
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
		
	if( Gl * L < 1e-6 ) { double item = pow( m * L / 2 , 2 ); return item; }
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
	int indexT2;
	int indexX2;
	for( int i = 1; i < 200; i++ ) {
		if( (z2[0][i] < T) and (z2[0][i+1] > T) ) { indexT2 = i; break; }
	}
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
	double P = pow(m,8) * pow(wB,2) * ( 1 + exp(-Gt*L) - 2*exp(-0.5*Gt*L) ) / ( 2 * pow(w,4) * pow(Gt,2) * ( pow( pow(w,2) - pow(m,2) , 2 ) + pow(w*Gl,2) ) );

	// low m limit
	double Plow = pow(m,8) * pow(wB,2) * pow(L,2) / ( 4 * pow(w,4) * ( pow(w,2) + pow(Gl,2) ) );

	if ( Gt*L < 1e-6 ) { return Plow; }
	else { return P; }
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

double pureLintegrate( double m, vector<vector<double>> z2, vector<double> T,
						vector<double> wp, vector<double> r, double L, double B ) {

	double total = 0;	// initiate value of sum at 0
	int len = r.size();
	
	// integrate by trapezium rule over r array
	for ( int c = 0; c < len - 1; c++ ) {
	
		int j = len - c - 2;
	
		if ( wp[j] <= m ) { continue; }	// only allow when energy greater than mass
		//if ( wp[j+1] > 10 ) { continue; }	// only res conversion up to 5eV allowed
		
		else {
			
			if ( wp[j+1] > 300 ) { cout << wp[j+1] << endl; }

			double dr = abs(r[j+1] - r[j]);
			//double height = 0.5 * ( ( PpureL( m, wp[j+1], L) * pureLintegrand( m, T[j+1], wp[j+1], r[j+1] ))
			//	+ ( PpureL( m, wp[j], L) * pureLintegrand( m, T[j], wp[j], r[j] ) ) );
			double height = 0.5 * ( ( PpureLConversion( wp[j+1], m, L, z2, B ) * pureLintegrand( m, T[j+1], wp[j+1], r[j+1] ))
				+ ( PpureLConversion( wp[j], m, L, z2, B ) * pureLintegrand( m, T[j], wp[j], r[j] ) ) );
			double dA = dr * height;
			
			// only add if real
			if ( isnan(dA) ) { continue; }
			else { total += dA; }
			}
	}
	return total;

}


// now to run the integral
void pureL( vector<vector<double>> z2, vector<double> T, vector<double> wp,
	 vector<double> r, double L, double phi, string name, double B ) {

	// implement new interrupt with save
	//signal( SIGINT, interrupt );
	
	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;

	string path = "data/limits/";
	string ext = "-pureL.dat";
	
	for ( double m = 1e-4; m < 1e1; m*=2 ) {

		// check if interrupt
		if( savenquit ){
		// write out
			write2D( path + name + "-INT" + ext, massIAXO, chiIAXO );
			// wait for other processes to save too
			cout << "waiting 10s to let other threads finish" << endl;
			sleep(10);
			exit(SIGINT);
		}

		double entryIAXO = pureLintegrate( m, z2, T, wp, r, L, B );
		double chi4IAXO = phi / entryIAXO;
		cout << name << ":	m = " << m << "	chi = " << pow( chi4IAXO, 0.25 ) << endl;
		chiIAXO.push_back( pow( chi4IAXO, 0.25 ) );
		massIAXO.push_back( m );
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
	double high = low + 0.01*rSolar;
	
	// perform integration by looping over r values
	for ( int c = 0; c < len - 1; c++ ) {

		if ( r[c] < low ) { continue; }
		else if ( r[c] > high ) { continue; }

		double dr = r[c+1] - r[c];	// define trapezium spacing
		double height = 0.5 * ( integrandREST(w, T[c], wp[c], r[c] ) 
			+ integrandREST(w, T[c+1], wp[c+1], r[c+1] ) );
		double dA = abs(dr * height);
		
		// only add if real
		if ( isnan(dA) ) { continue; }
		else { total += dA; }
	}
		
	return total;
		
}

// integral over r
double trapezeREST2( double w, vector<double> n, vector<double> T, vector<double> wp,
	  vector<double> r, vector<double> nH, vector<double> nHe4, vector<double> nHe3,
	  vector<vector<double>> z1, vector<vector<double>> z2, double low, int select ) {

	int len = r.size();	// get length of vector

	double total = 0;	// initiate value of sum at 0
	double high = low + 0.01*rSolar;
	
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

		// choose version
		if( select == 0 ) {
		double height = 0.5 * ( supREST(w, n[c], T[c], wp[c], r[c], nH[c], nHe4[c], nHe3[c], g1, g2) 
			+ supREST(w, n[c+1], T[c+1], wp[c+1], r[c+1], nH[c+1], nHe4[c+1], nHe3[c+1], g1, g2) );
		double dA = abs(dr * height);		
		// only add if real
		if ( isnan(dA) ) { continue; }
		else { total += dA; }
		}

		else if( select == 1 ) {
		total = 0.01 * rSolar * resREST(w, n[c], T[c], wp[c], r[c], nH[c], nHe4[c], nHe3[c], g1, g2);
		break;
		}

		else if( select == 2 ) {
		double height = 0.5 * ( unsupREST(w, n[c], T[c], wp[c], r[c], nH[c], nHe4[c], nHe3[c], g1, g2) 
			+ unsupREST(w, n[c+1], T[c+1], wp[c+1], r[c+1], nH[c+1], nHe4[c+1], nHe3[c+1], g1, g2) );
		double dA = abs(dr * height);		
		// only add if real
		if ( isnan(dA) ) { continue; }
		else { total += dA; }
		}

		else { cout << "error: select 0, 1 or 2..." << endl; }

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
	for ( double w = wlow; w < whigh - dw; w+=dw ) {
	
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

// full integration over omega
double integrateREST2( vector<double> n, vector<double> T, vector<double> wp,
	  vector<double> r, vector<double> nH, vector<double> nHe4, vector<double> nHe3,
	  vector<vector<double>> z1, vector<vector<double>> z2, double wlow, double rlow, int select ) {
	
	double total = 0;	// initiate value of sum at 0

	// integrate by trapezium rule over w array
	//double dw = 1e2;
	//for ( double w = 1e2; w < 1e5 - dw; w+=dw ) {
	double whigh = wlow + 100;	// eV
	double dw = 1e0;
	for ( double w = wlow; w < whigh - dw; w+=dw ) {
			double height = 0.5 * ( trapezeREST2( w+dw, n, T, wp, r, nH, nHe4, nHe3, z1, z2, rlow, select ) 
				+  trapezeREST2( w, n, T, wp, r, nH, nHe4, nHe3, z1, z2, rlow, select ) );
			double dA = abs(dw * height);

			// only add if real
			if ( isnan(dA) ) { continue; }
			else { total += dA; }
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
	string name = "data/flux_full.dat";
	//string name = "data/flux_m" + to_string(intm) + "_X" + to_string(intchi) + "-again.dat";

	for ( double rlow = 0; rlow < 0.99*rSolar; rlow += 0.01*rSolar ) {
		for ( double wlow = 1; wlow < 20e3; wlow += 100 ) {
	
//			double avg = integrateREST2(n, T, wp, r, nH,
//						 nHe4, nHe3, z1, z2, wlow, rlow, select ) / 0.1;	// eV3 keV-1
			double avg = integrateREST( T, wp, r, wlow, rlow ) / 0.1;	// eV3 keV-1			
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
		cout << "\33[A\33[2K\r" << line << "\% complete..." << endl;
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
	double high = low + 0.01*rSolar;
	
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
	return total;	// return average over bin
}

// average Gamma over energy bin
double integrateG( vector<double> n, vector<double> T,
	  vector<double> r, vector<double> nH, vector<double> nHe4, vector<double> nHe3,
	  vector<vector<double>> z1, vector<vector<double>> z2, double wlow, double rlow ) {
	
	double total = 0;	// initiate value of sum at 0

	// integrate by trapezium rule over w
	double whigh = wlow + 100;	// eV
	double dw = 1;
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
	string name = "data/wG-avg2.dat";

	for ( double rlow = 0; rlow < 0.99*rSolar; rlow += 0.01*rSolar ) {
		for ( double wlow = 1; wlow < 20e3; wlow += 100 ) {

			double avg = integrateG( n, T, r, nH,
						 nHe4, nHe3, z1, z2, wlow, rlow ) / 100 / 0.01*rSolar;	// eV			
			wG.push_back( avg );
			//cout << avg << endl;
		}
		writeREST( name, wG, line );
		line++;
		wG.clear();
		// note here:	\33[A moves up a line
		// 				\33[2K deletes line
		// 				\r goes to start of line
		cout << "\33[A\33[2K\r" << line << "\% complete..." << endl;
	}
	cout << "	creation of file " << name << " completed!" << endl;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// SPECTRA ETC ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// function to output E spectrum for T-plasmon
void spectrum( double m ) {
	
	//set chi
	double chi = 1e-11;

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
	
	for ( double w = 0; w < 2e4; w += 100 ) {

		if ( w <= m ) { continue; }	// only for m < wp
		
		// compute dphi/dE
		double entry = pow(chi,2) * trapeze( w, m, n, T, wp, r, nH, nHe4, nHe3, z1, z2 );
		E.push_back(w);
		phi.push_back(entry);
		
	}
	

	int mint = (int)log10(m);
	string name = "data/Espectrum" + to_string( mint ) + ".dat";
	
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
			double entry = lMixingResIntegrand( m, ne[j], T[j], wp[j], r[j], rFrac[j], nH[j], nHe4[j], nHe3[j], g1, g2 );
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

		double G = Gamma( w, ne[j], T[j], nH[j], nHe4[j], nHe3[j], g1, g2 );	// eV
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
	string name1 = "data/mfp2.dat";
	string name2 = "data/GammaT2.dat";
	string name3 = "data/wGammaT2.dat";

	// get phi values for each w = wp for resonance
	for ( int j = 0; j < len; j++ ) { 
	
	if ( 0.01 * rval > r[j] ) { continue; }
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
	cout << "\33[A\33[2K\r" << ( rval +1 ) << "\% complete..." << endl;
//	cout << "\33[A\33[2K\r" << ( 100 * line / len ) << "\% complete..." << endl;
	mfp.clear();
	GammaT.clear();
	wGammaT.clear();
	}
	}

	while ( rval < 100 ) {
		vector<double> temp;
		for ( int c = 0; c < 200; c++ ) { temp.push_back(0); }
		writeREST( name1, temp, line );
		writeREST( name2, temp, line );
		writeREST( name3, temp, line );
		line++;
		rval++;
	}

	cout << "creation of files:\n	" << name1 << "\n	" << name2 << "\n	" << name3 << "\ncompleted!" << endl;
}

// the end :)
