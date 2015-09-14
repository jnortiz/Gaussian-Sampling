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

    int parameter_set_id = 2;
    
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
    
    timestamp_t avgSetup = 0.0, avgGSO = 0.0, avgGaussianSampler = 0.0;            

    RR sigmaRR;
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
        cout << "[!] Setup running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s.\n" << endl;
    }//end-for
    
    cout << "Setup average running time: " << (float)(avgSetup/((float)(nIterations)*1000000000.0)) << " s." << endl;   
    
    /* Short basis expansion */
    mat_ZZ S;
    mat_RR T, TTilde;
    vec_RR center, sample;    
    RR normTTilde;

    int i, j, length, n;
    n = hibe->GetN();

    hibe->GetSampler()->RotBasis(S, hibe->GetMsk(), n);

    length = S.NumRows(); // Square matrix
    T.SetDims(length, length);

    for(i = 0; i < T.NumRows(); i++)
        for(j = 0; j < T.NumCols(); j++)
            T[i][j] = to_RR(S[i][j]);
        
    /* Orthogonalization of short basis S - Usual procedure of (Lyubashevsky, and Prest, 2015), Algorithm 1 */
    for(int it = 0; it < nIterations; it++) {
        
        ts_start = get_timestamp();
        normTTilde = hibe->GetSampler()->GramSchmidtProcess(TTilde, T, precision);
        ts_end = get_timestamp();            

        avgGSO += (ts_end - ts_start);
        
        cout << "[!] Norm of Gram-Schmidt reduced basis: " << normTTilde << endl;            
        cout << "[!] Gram-Schmidt orthogonalization process running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s.\n" << endl;
        
    }//end-for
    
    cout << "[!] Gram-Schmidt orthogonalization average time: " << (float)(avgGSO/((float)(nIterations)*1000000000.0)) << " s.\n" << endl;               
    
    T.kill();
    
    /* Parameters of distribution over the lattice span by A and T_A */
    sigmaRR = normTTilde*(log(n)/log(2)); // sigma >= norm(ATilde)*omega(sqrt(log(n)))
    
    center.SetLength(length);
    for(i = 0; i < center.length(); i++)
        center[i] = to_RR(0);    

    for(int it = 0; it < nIterations; it++) {

        /* Sampling from the discrete Gaussian distribution D_{\Lambda, \sigma, c} - (Lyubashevsky, and Prest, 2015), Algorithm 8 */
        ts_start = get_timestamp();    
        sample = hibe->GetSampler()->GaussianSamplerFromLattice(S, TTilde, sigmaRR, precision, tailcut, center);            
        ts_end = get_timestamp();    

        avgGaussianSampler += (ts_end - ts_start);
        cout << "\n[!] GaussianSampler running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s.\n" << endl;               
 
        cout << "\nSample from the lattice: " << sample << endl;

    }//end-for

    cout << "Gaussian-Sampler average running time: " << (float)(avgGaussianSampler/((float)(nIterations)*1000000000.0)) << " s.\n" << endl;

    delete(hibe);
    
    return 0;
    
}//end-main() 