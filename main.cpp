/* 
 * File:   main.cpp
 * Author: jnortiz
 *
 * Created on March 10, 2015, 2:12 PM
 * 
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
}//end-get_timestamp()

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
    
    timestamp_t avgSetup = 0.0, avgRefreshing = 0.0, avgPeikert = 0.0;            
    timestamp_t avgKlein = 0.0, avgUsual = 0.0;
    
    int nIterations, security_level, tailcut;
    long precision;

    security_level = 128;    
    /* (Saarinen, 2015) argues that half of the desired security level is almost always sufficient. */
    precision = security_level/2; 
    tailcut = 13;
    nIterations = 1;

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
        
    mat_ZZ S;
    RR BTilde_norm, normOfB, R, s;
    
    int m = hibe->GetM();
    int n = hibe->GetN();
    int length = m*n;
    
    /* Short basis expansion using the Rot_f operator */
    hibe->GetSampler()->RotBasis(S, hibe->GetMsk(), hibe->GetN());
    
    /* Getting the norm of S */
    mat_RR B;
    NTL::conv(B, S);    
    normOfB = hibe->GetSampler()->NormOfBasis(B);
        
    /* Getting the [norm of the] Gram-Schmidt orthogonalization of S */
    mat_RR BTilde;
    ts_start = get_timestamp();
    //BTilde_norm = hibe->GetSampler()->GramSchmidtProcess(BTilde, B, precision);
    ts_end = get_timestamp();    
    B.kill();
    
    cout << "[!] GramSchmidt Process running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;    
    cout << "[!] Norm of the orthogonal basis: " << BTilde_norm << endl;
    
    vec_RR center;
    center.SetLength(length);
    for(int i = 0; i < length; i++)
        center[i] = to_RR(NTL::RandomBnd(q));
    
    RR sigmaRR = log(length)*BTilde_norm;
    vec_RR sample;
    
    nIterations = 0;    
    for(int it = 0; it < nIterations; it++) {

        ts_start = get_timestamp();    
        sample = hibe->GetSampler()->Klein(S, BTilde, sigmaRR, precision, tailcut, center);
        ts_end = get_timestamp();

        avgKlein += (ts_end - ts_start);

        cout << "[!] Klein's algorithm running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;    
        cout << "[>] Sample from the lattice: " << sample << endl;

        ts_start = get_timestamp();    
        sample = hibe->GetSampler()->GaussianSamplerFromLattice(S, BTilde, sigmaRR, precision, tailcut, center);
        ts_end = get_timestamp();

        avgUsual += (ts_end - ts_start);
        
        cout << "[!] Usual Gaussian sampler running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;    
        cout << "[>] Sample from the lattice: " << sample << endl;

    }//end-for
    
    BTilde.kill();
    center.kill();
    sample.kill();
    
    cout << endl;
    if(nIterations > 1) {
        cout << "[!] Klein's average running time: " << (float)(avgKlein/((float)(nIterations)*1000000000.0)) << " s." << endl;
        cout << "[!] Usual Gaussian sampler average running time: " << (float)(avgUsual/((float)(nIterations)*1000000000.0)) << " s.\n" << endl;
    }//end-if    
    
    /* Computing parameters r and s of Peikert's offline phase */
    R = log(length)/log(2);
    s = R*(R*normOfB + 1);    
    NTL::mul(s, s, s);

    RR factor;
    NTL::div(factor, to_RR(1), sqrt(2*NTL::ComputePi_RR()));
    
    mat_RR Sigma;
    Sigma.SetDims(length, length);
    for(int i = 0; i < Sigma.NumRows(); i++)
        NTL::mul(Sigma[i][i], s, factor);
            
    /* Computing the Peikert's algorithm offline phase */
    mat_RR B2;
    mat_ZZ Z;
    ts_start = get_timestamp();    
    int outputOfflinePeikert = hibe->GetSampler()->OfflinePeikert(Z, B2, S, q, R, Sigma, m*n, precision);    
    ts_end = get_timestamp();
    
    mat_ZZ_p Z_p;
    NTL::conv(Z_p, Z);
    Z.kill();
    Sigma.kill();
    
    cout << "[!] Offline phase of Peikert's algorithm running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;    
    
    RR v = hibe->GetSampler()->ZCreatePartition(64, factor, precision, to_RR(tailcut));
    
    mat_ZZ_p S_p;
    vec_ZZ x2;
    NTL::conv(S_p, S);    
    S.kill();    
    
    nIterations = 10;
    if(outputOfflinePeikert == 0) {        
        
        vec_ZZ_p c_p, x2_p;
        vec_ZZ center, sample;    
        center.SetLength(length);
        
        for(int it = 0; it < nIterations; it++) {
            
            // Getting a fresh vector x2
            ts_start = get_timestamp();        
            x2 = hibe->GetSampler()->RefreshPeikert(B2, R, v, m*n, precision);
            ts_end = get_timestamp();    
            
            cout << "[!] Refreshing phase of Peikert's algorithm running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;    

            avgRefreshing += (ts_end - ts_start);
            
            for(int i = 0; i < center.length(); i++)
                center[i] = RandomBnd(q);    

            NTL::conv(c_p, center);
            NTL::conv(x2_p, x2);
            
            /* Getting a sample from the lattice using the Peikert's algorithm */
            ts_start = get_timestamp();    
            sample = hibe->GetSampler()->Peikert(S_p, Z_p, c_p, x2_p, (long)q, R, precision);
            ts_end = get_timestamp();    

            avgPeikert += (ts_end - ts_start);
            
            cout << "[!] Peikert running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;
            cout << "[>] Sample from the lattice: " << sample << endl;
            
        }//end-for

        cout << endl;
        if(nIterations > 1) {
            cout << "Refreshing phase of Peikert's algorithm average running time: " << (float)(avgRefreshing/((float)(nIterations)*1000000000.0)) << " s." << endl;
            cout << "Peikert average running time: " << (float)(avgPeikert/((float)(nIterations)*1000000000.0)) << " s.\n" << endl;
        }//end-if
        
        c_p.kill();
        x2_p.kill();
        center.kill();
        sample.kill();
        
    }//end-if
    
    B2.kill();
    x2.kill();
    S_p.kill();
    Z_p.kill();
    
    delete(hibe);
    
    return 0;
    
}//end-main() 