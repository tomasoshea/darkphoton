// Tom O'Shea 2023

// output .dat file of solar dark photon flux spectrum
// for use in REST axionlib

#include "darkphoton.h"	// home of the real code

using namespace std;
#include <thread>

int main(){

    // m [eV] and chi as args
    thread t1(fluxREST, 1e-3, 1e-11 );
    thread t2(fluxREST, 1e-2, 1e-11 );
    thread t3(fluxREST, 1e-1, 1e-11 );
    thread t4(fluxREST, 1e0, 1e-11 );
    thread t5(fluxREST, 1e1, 1e-11 );
    thread t6(fluxREST, 1e2, 1e-11 );
    thread t7(fluxREST, 1e3, 1e-11 );

    t1.join();
    t2.join();
    t3.join();
    t4.join();
    t5.join();
    t6.join();
    t7.join();
	
	return(0);

}
