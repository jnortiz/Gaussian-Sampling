/* 
 * File:   main.cpp
 * Author: jnortiz
 *
 * Created on March 10, 2015, 2:12 PM
 */

#include <cstdlib>
#include <cstdint>
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

    /* 128-bit security level */
/*
    h = 10;
    k = 7; // n = 2^k is the degree of polynomials in R and R_0
    q = 2083; //q must be prime and congruent to 3 mod 8
    m1 = 13;
    m2 = 170; //m2 >= lambda*m1, such as lambda is the security parameter and lambda = ceil(1 + lg(q))
*/
    
    /* Intermediate parameter set */
/*
    h = 10;
    k = 6; // n = 2^k is the degree of polynomials in R and R_0
    q = 1019; //q must be prime and congruent to 3 mod 8
    m1 = 11;
    m2 = 122; //m2 >= lambda*m1, such as lambda is the security parameter and lambda = ceil(1 + lg(q))
*/
    
    /* Toy parameter set */

    h = 10;
    k = 2; // n = 2^k is the degree of polynomials in R and R_0
    q = 11; //q must be prime and congruent to 3 mod 8
    m1 = 5;
    m2 = 30; //m2 >= lambda*m1, such as lambda is the security parameter and lambda = ceil(1 + lg(q))
/**/
    
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
    
    int op = 2;
    
    switch(op) {
        case 0: { // Continuous Gaussian sampling using Ziggurat method
            
            /*
             * TODO:
             * Change ZCreatePartition to constructor.
             */
            
            RR output = hibe->GetSampler()->Ziggurat(128, to_RR(1)/sqrt(2*NTL::ComputePi_RR()), 107, to_RR(11));
            cout << "\nContinuous sample: " << output << endl;
            
            break;
        }
        case 1: { // Cholesky decomposition
//                        
//            
//            /* Basis generation phase - (Mochetti, and Dahab, 2014) and (Stehlé et al., 2009) */
//            hibe->Setup(h);
//            
//            Vec< Vec<double> > B, Sigma;
//            int i, j, m = hibe->GetM();
//            double r = sqrt(log(2*m));
//            r *= r;
//            
//            /*
//             * TODO:
//             * First, run Rot(S) in order to S be a square matrix of integers
//             */
//            
//            cout << "\n/** Input matrix **/" << endl;
//            Sigma.SetLength(m);
//            for(i = 0; i < m; i++) {
//                Sigma[i].SetLength(m);
//                for(j = 0; j < m; j++) {
//                    Sigma[i][j] = abs((double)(hibe->GetMsk()[i][j]));
//                    cout << Sigma[i][j] << " ";
//                }
//                cout << endl;
//            }//end-for                       
//            
//            hibe->GetSampler()->CholeskyDecomposition(B, Sigma, m);
//            
//            cout << "\n\n /** Decomposition **/" << endl;
//            for(i = 0; i < m; i++) {
//                for(j = 0; j < m; j++)
//                    cout << B[i][j] << " ";
//                cout << endl;
//            }//end-for            
//            
//            break;
//            
        }//end-case        
        case 2: { // Orthogonalization for the special case of block isometric basis and Gaussian sampling over lattices
            
            timestamp_t avgSetup, avgGSO, avgGaussianSampler;
            
            avgSetup = 0.0;
            avgGSO = 0.0;
            avgGaussianSampler = 0.0;            
            
            RR sigmaRR;
            RR zero;  
            int it, nIterations, tailcut;
            long precision;
            
            sigmaRR = to_RR(3.195);
            precision = 107;
            tailcut = 13;

            zero = to_RR(0);           
            
            /* Basis generation phase - (Mochetti, and Dahab, 2014) and (Stehlé et al., 2009) */
            ts_start = get_timestamp();
            hibe->Setup(h);
            ts_end = get_timestamp();            

            avgSetup += (ts_end - ts_start);                
            cout << "[!] Setup running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;

            /*
             * TODO:
             * 1. Expand the short basis;
             * 2. Run the orthogonalization procedure;
             *  2.1 Adapt the algorithm/basis. The block's disposition has changed.
             * 3. Sampling using hybrid Gaussian sampler of Ducas and Prest.
             */

             /* Short basis expansion */
            mat_ZZ S;
            mat_RR T, TTilde;
            RR normTTilde;

            int i, j, length, n;
            n = hibe->GetN();

            hibe->GetSampler()->RotBasis(S, hibe->GetMsk(), n);
            
            length = S.NumRows();
            T.SetDims(length, length);

            for(i = 0; i < T.NumRows(); i++)
                for(j = 0; j < T.NumCols(); j++)
                    T[i][j] = to_RR(S[i][j]);
                       
#ifdef DEBUG
            cout << "/** Master secret key **/" << endl;                
            cout << hibe->GetMsk() << endl;

            cout << "\n/** S as an integer matrix **/" << endl;
            for(i = 0; i < S.length(); i++) {
                cout << S[i] << endl;
                if((i+1) % n == 0)
                    cout << endl;
            }//end-for

            cout << "\n/** S as a R-module matrix **/" << endl;
            cout << outRotR << endl;
#endif

            /* Orthogonalization of short basis S - Usual procedure of (Lyubashevsky, and Prest, 2015), Algorithm 1 */
            ts_start = get_timestamp();
            normTTilde = hibe->GetSampler()->GramSchmidtProcess(TTilde, T, precision);
            ts_end = get_timestamp();            

            cout << "[!] Gram-Schmidt orthogonalization process running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;
            cout << "[!] Norm of Gram-Schmidt reduced basis: " << normTTilde << endl;            

#ifdef DEBUG                
            for(i = 0; i < TTilde.NumRows(); i++) {
                cout << TTilde[i];
                if((i+1) % TTilde.NumCols())
                    cout << endl;
            }//end-for
#endif

//            return 0;

            sigmaRR = normTTilde*(log(2*hibe->GetN())/log(2)) + 1;

            vec_RR center, sample;                
            hibe->GetSampler()->SetCenter(center, S);
//            cout << "[!] Center of distribution: " << center << endl;

            nIterations = 1;
            
            for(it = 0; it < nIterations; it++) {
            
                /* Sampling from the discrete Gaussian distribution D_{\Lambda, \sigma, c} - (Lyubashevsky, and Prest, 2015), Algorithm 8 */
                ts_start = get_timestamp();    
                sample = hibe->GetSampler()->GaussianSamplerFromLattice(S, TTilde, sigmaRR, precision, tailcut, center);            
                ts_end = get_timestamp();    

                avgGaussianSampler += (ts_end - ts_start);
                cout << "\n[!] GaussianSampler running time: " << (float)((ts_end - ts_start)/1000000.0) << " ms." << endl;                
                cout << "\nSample from the lattice: " << sample << endl;

            }//end-for

            cout << "Gaussian-Sampler average running time: " << (float)(avgGaussianSampler/((float)(nIterations)*1000000000.0)) << " s. " << endl;

            break;
            
        }//end-case
        
    }//end-switch
    
    delete(hibe);
    
    return 0;
    
}//end-main()