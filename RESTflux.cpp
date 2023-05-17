// Tom O'Shea 2023

// output .dat file of solar dark photon flux spectrum
// for use in REST axionlib

#include "darkphoton.h"	// home of the real code

using namespace std;
#include <thread>

int main(){

    // 0: suppressed, 1: resonant, 2: unsuppressed
    thread t0( fluxREST, 0, 1e2, 1e-11 );
    thread t1( fluxREST, 1, 1e2, 1e-11 );
    thread t2( fluxREST, 2, 1e2, 1e-11 );
//    thread t4( fluxREST, 0, 1e1, 1e-11 );
//    thread t5( fluxREST, 1, 1e1, 1e-11 );
//    thread t6( fluxREST, 2, 1e1, 1e-11 );


    t0.join();
    t1.join();
    t2.join();

   
    // m [eV] and chi as args
//    thread t1(fluxREST, 1e-3, 1e-11 );
//    thread t2(fluxREST, 1e-2, 1e-11 );
//    thread t3(fluxREST, 1e-1, 1e-11 );
//    thread t4(fluxREST, 1e0, 1e-11 );
///   thread t5(fluxREST, 1e1, 1e-11 );
//    thread t6(fluxREST, 1e2, 1e-11 );
//    thread t7(fluxREST, 1e3, 1e-11 );
//
//    t1.join();
//    t2.join();
//    t3.join();
//    t4.join();
//    t5.join();
//    t6.join();
//    t7.join();
//	
//    
//    fluxREST( 1e-3, 1e-11);
//    fluxREST( 1e-3, 1e-12);
//    fluxREST( 1e-2, 1e-11);
//    fluxREST( 1e-2, 1e-12);
//    fluxREST( 1e-1, 1e-11);
//    fluxREST( 1e-1, 1e-12);
//    fluxREST( 1e0, 1e-11);
//    fluxREST( 1e0, 1e-12);
//    fluxREST( 1e1, 1e-11);
//    fluxREST( 1e1, 1e-12);
//    fluxREST( 1e2, 1e-11);
//    fluxREST( 1e2, 1e-12);
//    fluxREST( 1e3, 1e-11);
//    fluxREST( 1e3, 1e-12);


	return(0);

}
