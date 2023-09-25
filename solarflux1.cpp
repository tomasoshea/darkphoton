// Tom O'Shea 2023

// get solar flux on earth w/o conversion prob

#include "darkphoton.h"

int main() {
    //for( double m = 1e-3; m < 2e4; m*=10 )
    thread t1(fluxmx,1e-3);
    thread t2(fluxmx,1e-2);
    thread t3(fluxmx,1e-1);
    thread t4(fluxmx,1e0);
    thread t5(fluxmx,1e1);
    thread t6(fluxmx,1e2);
    thread t7(fluxmx,1e3);
    thread t8(fluxmx,1e4);

    t1.join();
    t2.join();
    t3.join();
    t4.join();
    t5.join();
    t6.join();
    t7.join();
    t8.join();
    
    return 0;
}
