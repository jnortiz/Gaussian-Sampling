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

    timestamp_t ts_start, ts_end;
    int h, k;
    double lambda;
    int m1, m2, q, r;

    int parameter_set_id = 4;
    
    switch(parameter_set_id) {
        case 0: {
            /* 128-bit security level */
            h = 10;
            k = 7;
            q = 2083; // q must be prime and congruent to 3 mod 8
            m1 = 13;
            m2 = 170; // m2 >= lambda*m1
            break;
        }
        case 1: {
            /* Intermediate parameter set */
            h = 10;
            k = 6; 
            q = 1019;
            m1 = 11;
            m2 = 122;         
            break;
        }
        case 2: {
            /* Small parameter set */
            h = 10;
            k = 4;
            q = 499;
            m1 = 10;
            m2 = 101;                 
            break;
        }
        case 3: {
            /* Tiny parameter set */
            h = 10;
            k = 2;
            q = 11;
            m1 = 5;
            m2 = 30;
            break;
        }
        case 4: {
            h = 10;
            k = 1;
            q = 3;
            m1 = 3;
            m2 = 9;
            break;
        }
        default: {
            cout << "Please, select a parameter set." << endl;
            return -1;
            break;
        }
    }//end-switch    

    r = ceil((double)(1 + (log(q)/log(3))));
    lambda = ceil(1 + (log(q)/log(2)));

    if(q % 8 != 3) {
        cout << "q must be congruent to 3 mod 8.\n";
        return -1;
    }//end-if

    if(m1 < r) {
        cout << "m1 must be greater or equal to r.\n";
        return -1;        
    }//end-if

    if(m1 < lambda) {
        cout << "m1 must be greater or equal to lambda.\n";
        return -1;
    }//end-if

    if(m2 < lambda) {
        cout << "m2 must be greater or equal to lambda.\n";
        return -1;
    }//end-if

    if(m2 < lambda*m1) {
        cout << "m2 must be greater or equal to lambda times m1.\n";
        return -1;
    }//end-if
    
    HIBE *hibe = new HIBE(q, m1, m2, k); // Parameters required only in Gaussian sampling from lattices  
    
    timestamp_t avgSetup = 0.0, avgSampleD = 0.0;            

    int nIterations, tailcut;
    long precision;

    nIterations = 1;
    precision = 107;
    tailcut = 13;

    /* Basis generation phase - (Mochetti, and Dahab, 2014) and (Stehlé et al., 2009) */
    for(int it = 0; it < nIterations; it++) {
        ts_start = get_timestamp();
        hibe->Setup(h);
        ts_end = get_timestamp();            
        
        avgSetup += (ts_end - ts_start);
        cout << "[!] Setup running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;
    }//end-for
    
    cout << endl;
    if(nIterations > 1)
        cout << "Setup average running time: " << (float)(avgSetup/((float)(nIterations)*1000000000.0)) << " s." << endl;   
        
    mat_RR B, B2, Sigma;
    mat_ZZ S, Z;
    vec_ZZ x2;
    RR normOfB, R, s;
    
    int m = hibe->GetM();
    int n = hibe->GetN();
    
    /* Short basis expansion */
    hibe->GetSampler()->RotBasis(S, hibe->GetMsk(), hibe->GetN());
    
    /* Getting the norm of S */
    NTL::conv(B, S);    
    normOfB = hibe->GetSampler()->NormOfBasis(B);
    B.kill();
    
    /* Computing parameters r and s of Peikert's offline phase */
    R = log(m*n)/log(2);
    s = R*(2*normOfB + 1) + 1;    
    NTL::mul(s, s, s);
    
    Sigma.SetDims(m*n, m*n);
    for(int i = 0; i < Sigma.NumRows(); i++)
        Sigma[i][i] = s;
        
    /* Computing the Peikert's algorithm offline phase */
    hibe->GetSampler()->OfflineSampleD(Z, B2, S, q, R, Sigma, m*n);    
    Sigma.kill();
    x2 = hibe->GetSampler()->RefreshSampleD(B2, R, m*n);

    vec_ZZ center;    
    center.SetLength(S.NumRows());
    for(int i = 0; i < center.length(); i++)
        center[i] = RandomBnd(q);    
    
    /* Getting a sample from the lattice using the Peikert's algorithm */
    vec_ZZ sample;
    
    for(int it = 0; it < nIterations; it++) {
        ts_start = get_timestamp();    
        sample = hibe->GetSampler()->SampleD(S, Z, center, x2, q, R);
        ts_end = get_timestamp();    

        avgSampleD += (ts_end - ts_start);
        cout << "[!] SampleD running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;
        cout << "[>] Sample from lattice: " << sample << endl;
    }//end-for
    
    if(nIterations > 1)
        cout << "SampleD average running time: " << (float)(avgSampleD/((float)(nIterations)*1000000000.0)) << " s.\n" << endl;
    
    center.kill();
    S.kill();
    Z.kill();
    x2.kill();
    
    delete(hibe);
    
    return 0;
    
}//end-main() 