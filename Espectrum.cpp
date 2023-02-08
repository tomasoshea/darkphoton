// Tom O'Shea 2022

// get values to plot flux against E for m values

// based on Javier Redondo's paper https://arxiv.org/abs/0801.1527

#include "darkphoton.h"	// home of the real code
#include <thread>

using namespace std;

int main(){
	
	// parellel run
	//for_each( execution::par, masses.begin(), masses.end(), spectrum );

	// non-parellel run
	for( double m = 1e-5 ; m < 2e5 ; m = m*10 ) { spectrum( m ); }
	
	return(0);

}
