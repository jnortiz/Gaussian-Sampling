/* 
 * File:   main.cpp
 * Author: jnortiz
 *
 * Created on March 10, 2015, 2:12 PM
 */

#include <cstdlib>
#include <math.h>
#include "HIBE.h"

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp() {
    struct timespec now;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &now);
    return now.tv_nsec + (timestamp_t)now.tv_sec * 1000000000.0;
}

int main(void) {
        
    int action = 2;
    
    int h, k;
    double lambda;
    int m1, m2, q, r;

    h = 10;
    k = 3; // n = 2^k is the degree of polynomials in R and R_0
    q = 19; //q must be prime and congruent to 3 mod 8
    m1 = 6;
    m2 = 49; //m2 >= lambda*m1, such as lambda is the security parameter and lambda = ceil(1 + lg(q))
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

    HIBE *hibe = new HIBE(q, m1, m2, k);    
    
    RR sigmaRR = to_RR(3.195); // Standard deviation
    RR c = to_RR(0); // Center of the distribution            
    int tailcut = 13;
    int precision = 107;
    int nRectangles = 63; // Parameter of Ziggurat algorithm
    int omega = precision; // Parameter of Ziggurat algorithm
    
    switch(action) {
        case 1: {
            /* Parameter set from (Roy et al., 2013). That depends on the cryptographic system requirements */
            
            if(sigmaRR*to_RR(tailcut) > power2_RR(sizeof(int)*8+1)-1) {
                cout << "Error! This distribution can not be simulated. Aborting..." << endl;
                return -1;
            }//end-if                
            
            Vec<int> ZigguratPoly, KnuthPoly;
            int nSamples = 25000; // #coefficients in the polynomial
            
            timestamp_t ts_start, ts_end;
            
            ts_start = get_timestamp();            
            ZigguratPoly = hibe->GetSampler()->PolyGeneratorZiggurat(nSamples, nRectangles, sigmaRR, omega, precision, tailcut); // Coefficients, rectangles, sigma, omega and precision
            ts_end = get_timestamp();            
                        
            cout << "[!] Ziggurat running time for " << nSamples << " samples: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;
            cout << ZigguratPoly << endl;
            
            ts_start = get_timestamp();                        
            KnuthPoly = hibe->GetSampler()->PolyGeneratorKnuthYao(nSamples, precision, tailcut, sigmaRR, c); // Coefficients, precision, tailcut, and sigma
            ts_end = get_timestamp();            
            
            cout << "[!] Knuth-Yao running time for " << nSamples << " samples: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;
            cout << KnuthPoly << endl;
            
            break;
        }
        case 2: { // Case still under analysis
            
            // Creating the basis needed to sample from the lattice
            hibe->Setup(h); // Setup algorithm with h = 10
            
            Vec<ZZX> BTilde;
            Vec<ZZ> C;
            Vec<double> D;
            ZZX c; // Center of the lattice
            ZZ zero = to_ZZ(0);
            c.SetMaxLength(hibe->GetN());
            
            for(int i = 0; i < hibe->GetN(); i++) // Lattice centered in zero
                c[i] = zero;            
            
            // Lattice with random center
//            random(c, hibe->GetN()); // n = 2^{k-1} = \phi(m)           
            
            hibe->GetSampler()->FasterIsometricGSO(BTilde, C, D, hibe->GetA(), k);
            cout << "\n/* Basis A */" << endl;
            cout << hibe->GetA() << endl;
            
            cout << "\n/* Gram-Schmidt reduced basis */" << endl;
            cout << BTilde << endl;
            
            cout << "\nSample from the lattice: \n" << hibe->GetSampler()->GaussianSamplerFromLattice(hibe->GetA(), BTilde, sigmaRR, precision, tailcut, c, k) << endl;            
                        
            break;
        }
        case 3: {
            
            RR area, percentage;
            
            hibe->GetSampler()->DZCreatePartition(nRectangles, sigmaRR, precision, tailcut);
            
            area = hibe->GetSampler()->CoverageAreaZiggurat(sigmaRR);
            percentage = 100*(area/to_RR(0.5));
            
            cout << "Coverage area: " << area << endl;
            cout << "Percentage: " << percentage << endl; 
            
            break;
        }
        default:
            break;            
    }//end-switch
    
    delete(hibe);
    return 0;
    
}