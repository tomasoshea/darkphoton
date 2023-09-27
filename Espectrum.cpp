// Tom O'Shea 2022

// get values to plot flux against E for m values for t-DPs

#include "darkphoton.h"	// home of the real code
#include <thread>

using namespace std;

int main(){	

	//for( double m = 1e-5 ; m < 2e5 ; m = m*10 ) { spectrum( m ); }

	vector<double> wp = read("data/wp.dat");	// plasma frequency [eV]
	double max = *max_element( wp.begin(), wp.end() );
	spectrum(max);
	
	return(0);

}
