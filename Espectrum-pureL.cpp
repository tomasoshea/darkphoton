// Tom O'Shea 2022

// get values to plot flux against E for m values for l-DPs

#include "darkphoton.h"	// home of the real code
#include <thread>

using namespace std;

int main(){

	for( double m = 1e-6 ; m < 3e2 ; m = m*10 ) { spectrumPureL( m ); }
	
	return(0);

}
