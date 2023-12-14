// Tom O'Shea 2023

// script to generate limits from IAXO photon flux
// using maximum likelihood method
// from Julia's thesis


#include "darkphoton.h"	// home of the real code
#include <thread>

using namespace std;

// set params
bool gas = false;
int samplesize = 1e3;
double days = 5;
double dE = 1e4;
double wmin = 1;	// eV


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
				cout << c << endl;
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
double integrateAtlas( vector<double> flux, double m, vector<double> wtab, double L, double pressure ) {
	
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
	//if( isnan(total) ) { return 0; }
	return pow(m,4)*total*pow(m2eV*100,2)*s2eV;
}


double likelihood( double g, double m, double n, double bg, vector<double> flux, vector<double> wtab,
					double L, double A, double t, double effO, double effD ) {
	double logL, mu;
	// only sum over densities for gas run
	if( gas == false ) {
		double N = pow(g,4)*integrateAtlas( flux, m, wtab, L, 0 )*A*t*effO*effD;
		//cout << N << endl;
		if( N == 0 ) { return 0; }
		mu = bg + N;
		if(n==0) { logL = -mu; }
		else if (n==1) { logL = log(mu) - mu + 1; }
		else{ logL = n * (log(mu) - log(n)) - mu + n; }
	}
	else{
		logL = 0;
		double TIAXO = 300*K2eV;
		for( double pm = 0.9*m; pm <= 1.1*m; pm+=0.02*m ) {
			double pressure = pm*pm*m_e*TIAXO/(8*pi*a);		// eV4
			double N = pow(g,4)*integrateAtlas( flux, m, wtab, L, pressure )*A*t*effO*effD;
			if( N == 0 ) { return 0; }
			mu = bg + N;
			if(n==0) { logL = -mu; }
			else if (n==1) { logL = log(mu) - mu + 1; }
			else{ logL = n * (log(mu) - log(n)) - mu + n; }
		}
	}
	//cout << "alive" << endl;
	return exp(logL);
}


double CL95( double m, double n, double bg, vector<double> flux, vector<double> wtab, double L,
			 double A, double t, double effO, double effD ) {
	// first, get total area of intg for normalisation
	double norm = 0;
	double mg = 1.1;
	double thresh = 1e-100;
	double g = 1e-20;
	while(true) {
		double dg = g*(mg-1);
		double L1 = likelihood(g, m, n, bg, flux, wtab, L, A, t, effO, effD);
		double L2 = likelihood(g*mg, m, n, bg, flux, wtab, L, A, t, effO, effD);
		norm += 0.5*dg*(L2+L1);
		if( 0.5*dg*(L2+L1) < thresh ) { break; }
		g*=mg;
		//cout << L1 + L2 << endl;
	}
	if( norm == 0 ) { return 1e200; }
	// now, get 95% CL
	double total = 0;
	g = 1e-20;
	/*while(true) {
		double dg = g*9;
		double L1 = likelihood(g, m, n, bg, flux, wtab, L, A, t, effO, effD);
		double L2 = likelihood(g*10, m, n, bg, flux, wtab, L, A, t, effO, effD);
		total += 0.5*dg*(L2+L1);
		if( total > 0.95*norm ) {
			total -= 0.5*dg*(L2+L1);
			break;
		}
		else { g*=10; }
	}*/
	while(true) {
		double dg = g*(mg-1);
		double L1 = likelihood(g, m, n, bg, flux, wtab, L, A, t, effO, effD);
		double L2 = likelihood(g*mg, m, n, bg, flux, wtab, L, A, t, effO, effD);
		total += 0.5*dg*(L2+L1);
		g*=mg;
		if( total > 0.95*norm ) { break; }
	}
	
	return g;
}



// void functions for splitting jobs
void statsAtlas( vector<vector<double>> flux, vector<double> mtab, vector<double> wtab,
				string name, int detector ) {
	
	// define vectors
	vector<double> massIAXO;
	vector<double> chiIAXO;

	// set path for writeout
	string path = "data/limits/statsAtlas-";
	string ext = ".dat";
	
	// select IAXO version
	double L,A,area,bkg,t,effD,effO;
	if( detector == 0 ) {			// babyIAXO
		L = 10/m2eV;								// eV-1
		A = 0.77/(m2eV*m2eV);						// eV-2
		area = 0.6/(1e4*m2eV*m2eV);					// eV-2
		bkg = 1e-7*(1e4*1e-3*(m2eV*m2eV*s2eV));		// eV2
		if(gas) { t = days*24*3600/s2eV; }			// eV-1
		else{ t = 1.5*365*24*3600/s2eV; }			// eV-1
		t/=2;		// only data taking half the time
		effD = 0.7;
		effO = 0.35;
	}
	else if( detector == 1 ) {		// IAXO
		L = 20/m2eV;								// eV-1
		A = 2.3/(m2eV*m2eV);						// eV-2
		area = 8*0.15/(1e4*m2eV*m2eV);				// eV-2
		bkg = 1e-8*(1e4*1e-3*(m2eV*m2eV*s2eV));		// eV2
		if(gas) { t = days*24*3600/s2eV; }			// eV-1
		else{ t = 3*365*24*3600/s2eV; }				// eV-1
		t/=2;		// only data taking half the time
		effD = 0.8;
		effO = 0.7;
	}
	else if( detector == 2 ) {		// IAXO+
		L = 22/m2eV;								// eV-1
		A = 3.9/(m2eV*m2eV);						// eV-2
		area = 8*0.15/(1e4*m2eV*m2eV);				// eV-2
		bkg = 1e-9*(1e4*1e-3*(m2eV*m2eV*s2eV));		// eV2
		if(gas) { t = days*24*3600/s2eV; }			// eV-1
		else{ t = 5*365*24*3600/s2eV; }				// eV-1
		t/=2;		// only data taking half the time
		effD = 0.8;
		effO = 0.7;
	}
	
	// get bkg count from flux
	double bgCount = bkg*area*t*dE;
	cout << bgCount << endl;
	
	// get "observed" count from poisson
	default_random_engine generator;
	poisson_distribution<int> distribution(bgCount);

	// get for many DP mass values
	//cout << "alive!" << endl;
	for ( int c = 0; c < mtab.size(); c++ ) {
		//if (c<70){continue;}
		double total = 0;
		for( int i = 0; i < samplesize; i++ ) {
			int n = distribution(generator);
			total += CL95( mtab[c], n, bgCount, flux[c], wtab, L, A, t, effO, effD );
			//cout << i << endl;
			//cout << total << endl;
		}
		double chi = total/samplesize;
		chiIAXO.push_back(chi);
		massIAXO.push_back(mtab[c]) ;
		cout << name << ":	m = " << mtab[c] << "	chi = " << chi << endl;
	}
	if(gas) { ext = "-gas.dat"; }
	// write out
	write2D( path + name + ext, massIAXO, chiIAXO );
}





// argc counts number of command line inputs (including the script itself)
// argv is a vector containing cmd args separated by space (including script itself)
int main( int argc, char** argv ) {
	
	string suffix = "-1eV";
	
	// read csv files to vectors
	vector<double> wtab = read("data/wtab.dat");
	vector<double> mtab = read("data/mtab.dat");
	vector<vector<double>> flux = readAtlas("data/output_T.dat");	// cm-2 s-1

	// babyIAXO
	string nameBaby = "babyIAXO" + suffix;
	thread t2( statsAtlas, flux, mtab, wtab, nameBaby, 0 );	

	
	//IAXO baseline
	string nameBaseline = "baselineIAXO" + suffix;
	thread t4(  statsAtlas, flux, mtab, wtab, nameBaseline, 1 );


	/// IAXO upgraded
	string nameUpgraded = "upgradedIAXO" + suffix;
	thread t6(  statsAtlas, flux, mtab, wtab, nameUpgraded, 2 );


	// wait until all threads are finished
	t2.join();
	t4.join();
	t6.join();
	
	cout << "\nProcess completed!\n" << endl;
	
	return(0);
}



	
		
		
		
