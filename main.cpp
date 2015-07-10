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
    h = 10;
    k = 7; // n = 2^k is the degree of polynomials in R and R_0
    q = 2083; //q must be prime and congruent to 3 mod 8
    m1 = 13;
    m2 = 170; //m2 >= lambda*m1, such as lambda is the security parameter and lambda = ceil(1 + lg(q))

    /* Intermediate parameter set */
//    h = 10;
//    k = 6; // n = 2^k is the degree of polynomials in R and R_0
//    q = 1019; //q must be prime and congruent to 3 mod 8
//    m1 = 11;
//    m2 = 122; //m2 >= lambda*m1, such as lambda is the security parameter and lambda = ceil(1 + lg(q))
    
    // Toy parameter set:
//    h = 10;
//    k = 2; // n = 2^k is the degree of polynomials in R and R_0
//    q = 11; //q must be prime and congruent to 3 mod 8
//    m1 = 5;
//    m2 = 30; //m2 >= lambda*m1, such as lambda is the security parameter and lambda = ceil(1 + lg(q))

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
    
    /* Parameter set from (Roy et al., 2013). That depends on the cryptographic system requirements */
    RR sigmaRR = to_RR(3.195); // Standard deviation
    RR c = to_RR(0); // Center of the distribution            
    int tailcut = 13;
    int precision = 107;
    int nRectangles = 63; // Parameter of Ziggurat algorithm
    int omega = precision; // Parameter of Ziggurat algorithm
    
    switch(action) {
        case 1: {

            timestamp_t averageZiggurat, averageKnuthYao;
            
            averageZiggurat = 0.0;
            averageKnuthYao = 0.0;
            
            int i, nIterations = 10;
            
            Vec<int> ZigguratPoly, KnuthPoly;
            int nSamples = 1024; // #coefficients in the polynomial
            
            for(i = 0; i < nIterations; i++) {

                HIBE *hibe = new HIBE(q, m1, m2, k); // Parameters required only in Gaussian sampling from lattices  
                
                cout << endl;
                
                ts_start = get_timestamp();            
                ZigguratPoly = hibe->GetSampler()->PolyGeneratorZiggurat(nSamples, nRectangles, sigmaRR, omega, precision, tailcut); // Coefficients, rectangles, sigma, omega and precision
                ts_end = get_timestamp();            

                averageZiggurat += (ts_end - ts_start);
                
                cout << ZigguratPoly << endl;

                ts_start = get_timestamp();                        
                KnuthPoly = hibe->GetSampler()->PolyGeneratorKnuthYao(nSamples, precision, tailcut, sigmaRR, c); // Coefficients, precision, tailcut, and sigma
                ts_end = get_timestamp();            

                averageKnuthYao += (ts_end - ts_start);
                
                cout << KnuthPoly << endl;

                delete(hibe);
            
            }//end-for
            
            cout << "\nZiggurat average running time for " << nIterations << " iterations: " << (float)(averageZiggurat/((float)(nIterations)*1000000000.0)) << endl;
            cout << "Knuth-Yao average running time for " << nIterations << " iterations: " << (float)(averageKnuthYao/((float)(nIterations)*1000000000.0)) << endl;
            break;
        }
        case 2: {
            
            timestamp_t avgSetup, avgBlockGSO, avgGaussianSampler;
            
            avgSetup = 0.0;
            avgBlockGSO = 0.0;
            avgGaussianSampler = 0.0;
            
            Vec<ZZX> BTilde;
            Vec<ZZ_pX> out;                        
            Vec<ZZ_pX> B;            
            ZZX c, sample;
            RR normBTilde;
            ZZ zero;            
            int i, index, j, it, nIterations;
            
            zero = to_ZZ(0);            
            nIterations = 100;
            
            for(it = 0; it < nIterations; it++) {
                
                HIBE *hibe = new HIBE(q, m1, m2, k);    

                ts_start = get_timestamp();                        
                hibe->Setup(h);
                ts_end = get_timestamp();   

                avgSetup += (ts_end - ts_start);
                cout << "\n[!] Setup running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;
                
                B.SetLength(hibe->GetM()*hibe->GetN());
                
                /* Expansion of vector \tilde{a} into basis A */
                index = 0;
                for(i = 0; i < hibe->GetM(); i++) {
                    hibe->GetSampler()->rot(out, hibe->GetA()[i], hibe->GetN());
                    for(j = 0; j < hibe->GetN(); j++)
                        B[index++] = out[j];                
                }//end-for

                ts_start = get_timestamp();                        
                normBTilde = hibe->GetSampler()->BlockGSO(BTilde, B, hibe->GetM(), hibe->GetN());
                ts_end = get_timestamp();            

                cout << "\n[!] BlockGSO running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;
                avgBlockGSO += (ts_end - ts_start);
                
                cout << "[!] Norm of Gram-Schmidt reduced basis: " << normBTilde << endl;

                c.SetMaxLength(hibe->GetN());
                for(i = 0; i < hibe->GetN(); i++) // Lattice centered in zero
                    c[i] = zero;            

                sigmaRR = normBTilde*(log(hibe->GetN())/log(2)) + 1;

                ts_start = get_timestamp();    
                sample = hibe->GetSampler()->GaussianSamplerFromLattice(B, BTilde, sigmaRR, precision, tailcut, c, hibe->GetN());            
                ts_end = get_timestamp();    

                cout << "\n[!] GaussianSampler running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;
                avgGaussianSampler += (ts_end - ts_start);

                cout << "\nSample from the lattice: " << sample << endl;
                cout << endl;
                
                delete(hibe);
                
            }//end-for
            
            cout << "\n avgSetup: " << avgSetup << endl;
            cout << "\n avgBlockGSO: " << avgBlockGSO << endl;
            cout << "\n avgGaussianSampler: " << avgGaussianSampler << endl;

            cout << "Setup running time: " << (float)(avgSetup/1000000000.0) << " s."  << endl;
            cout << "Block-GSO running time: " << (float)(avgBlockGSO/1000000000.0) << " s." << endl;
            cout << "Gaussian-Sampler running time: " << (float)(avgGaussianSampler/1000000000.0) << " s. " << endl;
            
            cout << "Setup average running time: " << (float)(avgSetup/((float)(nIterations)*1000000000.0)) << " s."  << endl;
            cout << "Block-GSO average running time: " << (float)(avgBlockGSO/((float)(nIterations)*1000000000.0)) << " s." << endl;
            cout << "Gaussian-Sampler average running time: " << (float)(avgGaussianSampler/((float)(nIterations)*1000000000.0)) << " s. " << endl;
                        
            break;
        }
        case 3: {
            
            HIBE *hibe = new HIBE(q, m1, m2, k);    
            
            // Creating the basis needed to sample from the lattice
            ts_start = get_timestamp();                        
            hibe->Setup(h); // Setup algorithm with h = 10
            ts_end = get_timestamp();   
            
            cout << "\n[!] Setup running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;

            delete(hibe);
            
            break;
        }
        default:
            break;            
    }//end-switch
    
    return 0;
    
}