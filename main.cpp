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
    
    HIBE *hibe = new HIBE(q, m1, m2, k);  
    
    int m = hibe->GetM();
    int n = hibe->GetN();
    
    ZZ_pX f;
    f.SetLength(n+1);
    SetCoeff(f, 0, 1);
    SetCoeff(f, n, 1);    
    
    Samplers *samplers = new Samplers(q, f);
    
    timestamp_t avgSetup = 0.0, avgRefreshing = 0.0, avgPeikert = 0.0;            
    timestamp_t avgKlein = 0.0, avgUsual = 0.0;
    
    int nIterations, tailcut;
    long precision;

    // (Roy, Vercauteren and Verbauwhede, 2013): precision and sigma values for statistical distance less than 2^{-90}
    precision = 107;
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
    RR BTilde_norm, normOfB;
    
    int length = m*n;
    
    /* Short basis expansion using the Rot_f operator */
    samplers->RotBasis(S, hibe->GetMsk(), n);
    
    /* Getting the norm of S */
    mat_RR B;
    NTL::conv(B, S);    
    normOfB = samplers->NormOfBasis(B);
    
    int procedure_id = 4;
    
    switch(procedure_id) {
        
        case 0: {
            
            /*
             * Compact Gaussian sampler of (Lyubashevsky and Prest, 2015)
             */
            
            mat_RR OrthoBasis;
            vec_RR D;
                            
            ts_start = get_timestamp();    
            BTilde_norm = samplers->GramSchmidtProcess(OrthoBasis, D, B, precision);
            ts_end = get_timestamp();

            cout << "[!] Gram-Schmidt process running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;    

            vec_RR I, v;            
            ts_start = get_timestamp();    
            samplers->PrepareToSampleCGS(v, I, OrthoBasis, D, B[0], precision);
            ts_end = get_timestamp();

            cout << "[!] PrepareToSampleCGS running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;    

            vec_RR Bn, center, sample;
            Bn = OrthoBasis[OrthoBasis.NumRows()-1];                        
            OrthoBasis.kill();
            
            center.SetLength(length);
            for(int i = 0; i < center.length(); i++)
                center[i] = to_RR(0);                
            
            RR sigma = log(length)*BTilde_norm;

            ts_start = get_timestamp();    
            sample = samplers->CompactGaussianSampler(B, center, Bn, v, I, D, sigma, precision);
            ts_end = get_timestamp();

            cout << "[!] Compact-Gaussian-Sampler running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;    
            cout << "\n[>] Sample from the lattice: " << sample << endl;

            B.kill();
            S.kill();
            center.kill();
            v.kill();
            I.kill();
            D.kill();
            Bn.kill();
            
            break;
            
        }//end-case-0
    
        case 1: {
            
            /*
             * Klein's Gaussian sampler
             */
                       
            /* Getting the [norm of the] Gram-Schmidt orthogonalization of S */
            mat_RR BTilde;
            ts_start = get_timestamp();
            BTilde_norm = samplers->GramSchmidtProcess(BTilde, B, precision);
            ts_end = get_timestamp();
            
            B.kill();

            cout << "[!] Gram-Schmidt process running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;    
            cout << "[!] Norm of the orthogonal basis: " << BTilde_norm << endl;
                        
            vec_RR center;
            center.SetLength(length);
            for(int i = 0; i < length; i++)
                center[i] = to_RR(NTL::RandomBnd(q));

            RR sigmaRR = log(length)*BTilde_norm;
            vec_RR sample;

            nIterations = 1;    
            for(int it = 0; it < nIterations; it++) {

                ts_start = get_timestamp();    
                sample = samplers->Klein(S, BTilde, sigmaRR, precision, tailcut, center);
                ts_end = get_timestamp();

                avgKlein += (ts_end - ts_start);

                cout << "[!] Klein's algorithm running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;    
                cout << "[>] Sample from the lattice: " << sample << endl;
                
            }//end-for

            BTilde.kill();
            center.kill();
            sample.kill();

            cout << endl;
            if(nIterations > 1)
                cout << "[!] Klein's average running time: " << (float)(avgKlein/((float)(nIterations)*1000000000.0)) << " s." << endl;
                        
            break;
            
        }//end-case-1
        
        case 2: {
            
            /*
             * Usual Gaussian sampler
             */
                       
            /* Getting the [norm of the] Gram-Schmidt orthogonalization of S */
            mat_RR BTilde;
            ts_start = get_timestamp();
            BTilde_norm = samplers->GramSchmidtProcess(BTilde, B, precision);
            ts_end = get_timestamp();
            
            B.kill();

            cout << "[!] Gram-Schmidt process running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;    
            cout << "[!] Norm of the orthogonal basis: " << BTilde_norm << endl;
                        
            vec_RR center;
            center.SetLength(length);
            for(int i = 0; i < length; i++)
                center[i] = to_RR(NTL::RandomBnd(q));

            RR sigmaRR = log(length)*BTilde_norm;
            vec_RR sample;

            nIterations = 1;    
            for(int it = 0; it < nIterations; it++) {

                ts_start = get_timestamp();    
                sample = samplers->GaussianSamplerFromLattice(S, BTilde, sigmaRR, precision, tailcut, center);
                ts_end = get_timestamp();

                avgUsual += (ts_end - ts_start);

                cout << "[!] Usual Gaussian sampler running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;    
                cout << "[>] Sample from the lattice: " << sample << endl;

            }//end-for

            BTilde.kill();
            center.kill();
            sample.kill();

            cout << endl;
            if(nIterations > 1)
                cout << "[!] Usual Gaussian sampler average running time: " << (float)(avgUsual/((float)(nIterations)*1000000000.0)) << " s.\n" << endl;
            
            break;
            
        }//end-case-2
        
        case 3: {
            
            /*
             * Peikert's method for Gaussian sampling from q-ary lattices 
             */
                        
            RR R, s;
            
            R = log(length)/log(2);
            s = R*(2*length*normOfB + 1);  
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
            int outputOfflinePeikert = samplers->OfflinePeikert(Z, B2, S, q, R, Sigma, length, precision);    
            ts_end = get_timestamp();

            mat_ZZ_p Z_p;
            NTL::conv(Z_p, Z);
            Z.kill();
            Sigma.kill();

            cout << "[!] Offline phase of Peikert's algorithm running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;    

            RR v = samplers->ZCreatePartition(64, factor, precision, to_RR(tailcut));

            mat_ZZ_p S_p;
            vec_ZZ x2;
            NTL::conv(S_p, S);    
            S.kill();    

            nIterations = 10;
            if(outputOfflinePeikert == 0) {        

                vec_ZZ_p c_p, x2_p;
                vec_ZZ center, sample;    
                center.SetLength(S_p.NumRows());

                for(int it = 0; it < nIterations; it++) {

                    // Getting a fresh vector x2
                    ts_start = get_timestamp();        
                    x2 = samplers->RefreshPeikert(B2, R, v, length, precision);
                    ts_end = get_timestamp();    

                    cout << "\n[!] Refreshing phase of Peikert's algorithm running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;    

                    avgRefreshing += (ts_end - ts_start);

                    for(int i = 0; i < center.length(); i++)
                        center[i] = RandomBnd(q);    

                    NTL::conv(c_p, center);
                    NTL::conv(x2_p, x2);

                    /* Getting a sample from the lattice using the Peikert's algorithm */
                    ts_start = get_timestamp();    
                    sample = samplers->Peikert(S_p, Z_p, c_p, x2_p, (long)q, R, precision);
                    ts_end = get_timestamp();    

                    avgPeikert += (ts_end - ts_start);

                    cout << "\n[!] Peikert running time: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;
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
            
            break;
            
        }//end-case-3
        
        case 4: {
            
            vec_RR samples;
            timestamp_t avgZiggurat = 0.0;
            int dimension = 8194, m = 128;            
            RR factor, v;
            
            NTL::div(factor, to_RR(1), sqrt(2*NTL::ComputePi_RR()));
            
            samples.SetLength(dimension);            
            
            nIterations = 1000;
            
            for(int it = 0; it < nIterations; it++) {
                cout << endl;
                ts_start = get_timestamp();
                v = samplers->ZCreatePartition(m, factor, precision, to_RR(tailcut));            
                for(int j = 0; j < dimension; j++)
                    samples[j] = samplers->Ziggurat(m, factor, precision, v);                
                ts_end = get_timestamp();
                avgZiggurat += (ts_end - ts_start);
                cout << "[>] Continuous samples:" << samples << endl;
            }//end-for
            
            if(nIterations > 1)
                cout << "\n[!] (Continuous) Ziggurat algorithm average running time: " << (float)(avgZiggurat/((float)(nIterations)*1000000000.0)) << " s." << endl;
            
            samples.kill();            
            break;
            
        }//end-case-4        
        default: {
            break;
        }
            
    }//end-switch
    
    delete(hibe);
    delete(samplers);
    
    return 0;
    
}//end-main() 
