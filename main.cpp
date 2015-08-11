/* 
 * File:   main.cpp
 * Author: jnortiz
 *
 * Created on March 10, 2015, 2:12 PM
 */

#include <cstdlib>
#include <math.h>
#include <NTL/RR.h>
#include "HIBE.h"

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp() {
    struct timespec now;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &now);
    return now.tv_nsec + (timestamp_t)now.tv_sec * 1000000000.0;
}

int main(void) {
        
    int action = 1;
    
    timestamp_t ts_start, ts_end;
    int h, k;
    double lambda;
    int m1, m2, q, r;

    /* 128-bit security level */
//    h = 10;
//    k = 7; // n = 2^k is the degree of polynomials in R and R_0
//    q = 2083; //q must be prime and congruent to 3 mod 8
//    m1 = 13;
//    m2 = 170; //m2 >= lambda*m1, such as lambda is the security parameter and lambda = ceil(1 + lg(q))

    /* Intermediate parameter set */
//        h = 10;
//        k = 6; // n = 2^k is the degree of polynomials in R and R_0
//        q = 1019; //q must be prime and congruent to 3 mod 8
//        m1 = 11;
//        m2 = 122; //m2 >= lambda*m1, such as lambda is the security parameter and lambda = ceil(1 + lg(q))
    
    /* Toy parameter set */
    h = 10;
    k = 2; // n = 2^k is the degree of polynomials in R and R_0
    q = 11; //q must be prime and congruent to 3 mod 8
    m1 = 5;
    m2 = 30; //m2 >= lambda*m1, such as lambda is the security parameter and lambda = ceil(1 + lg(q))

    r = ceil((double)(1 + (log(q)/log(3))));
    lambda = ceil(1 + (log(q)/log(2)));

    if(q % 8 != 3) {
        cout << "q must be congruent to 3 mod 8.\n";
        return -1;
    }

    if(/*m1 < sigma ||*/ m1 < r) {
        cout << "m1 must be greater or equal to sigma and r.\n";
        return -1;        
    }

    if(m1 < lambda) {
        cout << "m1 must be greater or equal to lambda.\n";
        return -1;
    }

    if(m2 < lambda) {
        cout << "m2 must be greater or equal to lambda.\n";
        return -1;
    }

    if(m2 < lambda*m1) {
        cout << "m2 must be greater or equal to lambda times m1.\n";
        return -1;
    }
    
    HIBE *hibe = new HIBE(q, m1, m2, k); // Parameters required only in Gaussian sampling from lattices  
        
    timestamp_t avgSetup;
            
    avgSetup = 0.0;
            
    /* Basis generation phase - (Mochetti, and Dahab, 2014) and (Stehlé et al., 2009) */
    ts_start = get_timestamp();
    hibe->Setup(h);
    ts_end = get_timestamp();            

    avgSetup += (ts_end - ts_start);                
    cout << "\n[!] Setup running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;

    delete(hibe);
    
    return 0;
    
}