// Tom O'Shea 2023

// script to generate a function of dark photon mixing parameter against dark photon mass
// to give a theoretical upper limit on these parameters that IAXO could achieve

// original t-plasmon integral for IAXO vacuum & gas runs

#include "darkphoton.h"	// home of the real code
#include <thread>

using namespace std;


// read in gaunt factors from matlab matrix files
vector<vector<double>> readAtlas( string name ) {

	//cout << "Reading file " << name << "..." << endl;

	// open file defined in argument
	fstream file;
	file.open(name,ios::in);
	
	char delim(' ');	// define delimiter for file parsing
	
	if (file.is_open()){   // checking whether the file is open
		vector<vector<double>> g2D; // 2D matrix to store rows
		vector<double> row;	// define vector to store input values and return
		int line = 0;	// counter for lines
		
		while( !file.eof() ) {  // loop until end of file'
			int c = 0;	// define counter for elements
			string temp;	// define temporary storage string
			vector<double> row;
			while(getline(file, temp, delim)){	// get one row
				if( temp == "\n" ) { continue; }
				row.push_back(stod(temp));	// append double to vector
				temp.clear();
				c++;
				if(c==4000) {
					g2D.push_back(row);
					row.clear();
					c=0;
					line++;
				}
			}
		}
		file.close();   // close the file object.
		return g2D;	// returns vector of values
	}
	else{ cout << "ERROR: couldn't open file " << name << endl ; vector<vector<double>> err = {{69.}} ; return err ; }
}


// integration over omega for gas run
double integrateAtlas( vector<double> flux, double m, vector<double> wtab, double L, double pressure, double wmin ) {
	
	double total = 0;	// initiate value of sum at 0
	for ( int c = 0; c < wtab.size()-1; c++ ) {
		if ( wtab[c] <= m ) { continue; }	// only allow when energy greater than mass
		if ( wtab[c] < wmin ) { continue; }
		double dw = wtab[c+1] - wtab[c];
		if( pressure == 0 ) {
			total += 0.5*dw*( ( Pvacuum(wtab[c+1],m,L) * flux[c+1] )
				+ ( Pvacuum(wtab[c],m,L) * flux[c] ) );
		}
		else {
			total += 0.5*dw*( ( Pgas(wtab[c+1],m,L,pressure) * flux[c+1] )
				+ ( Pgas(wtab[c],m,L,pressure) * flux[c] ) );
		}
	}
	return total*m*m*m*m*pow(m2eV*100,2)*s2eV;
}


// void functions for splitting jobs
void Atlas( vector<vector<double>> flux, vector<double> mtab, vector<double> wtab, double L, string name ) {
	
	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;
	double pressure = 0;
	double T = 300*K2eV;
	double wmin = 100;	// eV

	// set path for writeout
	string path = "data/limits/";
	string ext = ".dat";

	// get for many DP mass values
	for ( int c = 0; c < mtab.size(); c++ ) {
			pressure = mtab[c]*mtab[c]*m_e*T/(8*pi*a);		// [eV4]
			//pressure = 0.1*0.1*m_e*T/(8*pi*a);				// 0.1eV
			double phi = integrateAtlas(flux[c],mtab[c],wtab,L,pressure,wmin);
			chiIAXO.push_back(phi);
			massIAXO.push_back(mtab[c]) ;
			cout << name << ":	m = " << mtab[c] << "	phi X-4 = " << phi << endl;
	}
	if(pressure!=0) { ext = "-gas.dat"; }
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
	vector<double> wtab = read("data/wtab.dat");
	vector<double> mtab = read("data/mtab.dat");
	vector<vector<double>> flux = readAtlas("data/limits/output_T.dat");	// cm-2 s-1

	// babyIAXO
	string nameBaby = "babyIAXO-Atlas" + suffix;
	double L1 = 10 / m2eV;	// in eV^-1
	// multithread to run simultaneously
	thread t2( Atlas, flux, mtab, wtab, L1, nameBaby );	

	
	//IAXO baseline
	string nameBaseline = "baselineIAXO-Atlas" + suffix;
	double L2 = 20 / m2eV;	// in eV^-1
	// multithread to run simultaneously
	thread t4(  Atlas, flux, mtab, wtab, L2, nameBaseline );


	/// IAXO upgraded
	string nameUpgraded = "upgradedIAXO-Atlas" + suffix;
	double L3 = 22 / m2eV;	// in eV^-1
	// multithread to run simultaneously
	thread t6(  Atlas, flux, mtab, wtab, L3, nameUpgraded );


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
