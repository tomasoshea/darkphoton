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


// IAXO conversion probability for vacuum run
double P( double w, double m, double L ) {

	// calculate Delta P
	double dP = w - pow( pow(w,2) - pow(m,2) , 0.5 );
	
	// plug into sin2
	double arg = dP * L / 2;
	
	double item = 4 * pow( sin(arg) , 2 );
	
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


// conversion probability for gas resonance
double PgasStats( double w, double m, double L, vector<vector<double>> z2, double den ) {
	
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
	double p1 = pow(m,4) / ( pow( pow(den,2) - pow(m,2), 2 ) + pow(w*G,2) );
	double p2 = 1 + exp(- G*L) - (2 * exp(- G*L / 2));
	
	double item = p1 * p2;
	
	return item;
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
	double dw = 1e0;
	for ( double w = 1e3; w < 1e5 - dw; w+=dw ) {
	
		if ( w > m + 1000 ) { continue; }	// set integral cutoff
		if ( w <= m ) { continue; }	// only allow when energy greater than mass
		else if ( w > m + wRange ) { continue; }
		
		else {
		
			double height = 0.5 * ( ( P( w+dw, m, L ) * trapeze( w+dw, m, n, T, wp, r, nH, nHe4, nHe3, z1, z2 ) ) 
				+ ( P( w, m, L ) * trapeze( w, m, n, T, wp, r, nH, nHe4, nHe3, z1, z2 ) ) );
			double dA = abs(dw * height);
			
			// only add if real
			if ( isnan(dA) ) { continue; }
			else { total += dA; }
			//cout << "hola" << endl;
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
	double dw = 1e0;
	for ( double w = 30; w < 1e2 - dw; w+=dw ) {
	
		//if ( w > m + 1e3 ) { continue; }	// set integral cutoff
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
		vector<vector<double>> z1, vector<vector<double>> z2, double phi, string name ) {

	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;

	// set path for writeout
	string path = "data/limits/";
	string ext = ".dat";
	
	double min = *min_element( wp.begin(), wp.end() );
	double max = *max_element( wp.begin(), wp.end() );

	
	// suppressed section
	for ( double m = 1e-4; m < min; m*=2 ) {
        
		double entryIAXO = integrate( m, n, T, wp, r, nH, nHe4, nHe3, L, z1, z2 );
		double chi4IAXO = phi / entryIAXO;
		chiIAXO.push_back( pow( chi4IAXO, 0.25 ) );
		massIAXO.push_back( m );
		cout << name << ":	m = " << m << "	chi = " << pow( chi4IAXO, 0.25 ) << endl;
	}
	
	// resonant sector
	int len = wp.size();
	for ( int i = 0; i < len; i = i+5 ) {

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
		
		double entryIAXO = integrate( m, n, T, wp, r, nH, nHe4, nHe3, L, z1, z2 );
		double chi4IAXO = phi / entryIAXO;
		chiIAXO.push_back( pow( chi4IAXO, 0.25 ) );
		massIAXO.push_back( m );
		cout << name << ":	m = " << m << "	chi = " << pow( chi4IAXO, 0.25 ) <<endl;
	}

		// unsuppressed section
	for ( double m = max; m < 1e5; m*=1.01 ) {

		double entryIAXO = integrate( m, n, T, wp, r, nH, nHe4, nHe3, L, z1, z2 );
		double chi4IAXO = phi / entryIAXO;
		chiIAXO.push_back( pow( chi4IAXO, 0.25 ) );
		massIAXO.push_back( m );
		cout << name << ":	m = " << m << "	chi = " << pow( chi4IAXO, 0.25 ) << endl;
		}


	// write out
	write2D( path + name + ext, massIAXO, chiIAXO );
	
}


void integrateTgas( vector<double> n, vector<double> T, vector<double> wp,
		vector<double> r, vector<double> nH, vector<double> nHe4, vector<double> nHe3, double L,
		vector<vector<double>> z1, vector<vector<double>> z2, double phi, string name ) {

	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;

	// set path for writeout
	string path = "data/limits/";
	string ext = "-gas.dat";
	
	double min = *min_element( wp.begin(), wp.end() );
	
	// suppressed section
	for ( double m = 1e-3; m < min; m*=2 ) {

		double entryIAXO = integrateGas( m, n, T, wp, r, nH, nHe4, nHe3, L, z1, z2 );
		double chi4IAXO = phi / entryIAXO;
		chiIAXO.push_back( pow( chi4IAXO, 0.25 ) );
		massIAXO.push_back( m );
		cout << name << ":	m = " << m << "	chi = " << pow( chi4IAXO, 0.25 ) <<endl;
	}
	
	// resonant sector
	int len = wp.size();
	for ( int i = 0; i < len; i=i+5 ) {
	
		int j = len - i - 1;
		double m = wp[j];
		if( m > 5 ){ break; }
		
		double entryIAXO = integrateGas( m, n, T, wp, r, nH, nHe4, nHe3, L, z1, z2 );
		double chi4IAXO = phi / entryIAXO;
		chiIAXO.push_back( pow( chi4IAXO, 0.25 ) );
		massIAXO.push_back( m );
		cout<< name << ":	m = " << m << "	chi = " << pow( chi4IAXO, 0.25 ) <<endl;
	}

		// write out
	write2D( path + name + ext, massIAXO, chiIAXO );
	
}
