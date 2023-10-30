// Tom O'Shea 2023

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h> 
#include <filesystem>
#include <bits/stdc++.h>


using namespace std;

double K2eV = (8.617333262e-5);	// Kelvin to eV
double T = 300 * K2eV;	// detector at room temp [eV]
double Ry = 13.605693;	// Rydberg energy [eV]

int main( int argc, char** argv ) {

	// check suffix is given as argument
	if( argc =! 2 ){
		cout << "enter n_max as argument" << endl;
		return 1;
	}
	int nmax = stoi(argv[1]);

	long double sum = 0;
	for( int n = 1; n <= nmax; n++ ) {
		double Kn = 0;
		if(n>3) { Kn = (16*n*n*(n+(7/6)))/(3*pow(n+1,2)*(n*n + n + 0.5)); }
		else{ Kn = 1; }
		double En = Ry/(n*n);	// [eV]
		double wn = 1;
		sum += 2*n*n*wn*exp(En/T);
	}
	cout << sum << endl;
}
