/* 
 * File:   main.cpp
 * Author: jnortiz
 *
 * Created on March 10, 2015, 2:12 PM
 */

#include <cstdlib>
#include <math.h>
#include "HIBE.h"
#include "Samplers.h"

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp() {
    struct timespec now;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &now);
    return now.tv_nsec + (timestamp_t)now.tv_sec * 1000000000.0;
}

int main(void) {
        
    int action = 3;
    int k = 8; // n = 2^k is the degree of polynomials in R and R_0
    
    switch(action) {
        case 1: { // Setup algorithm from HIBE scheme
            double lambda;
            int q, m1, m2, r, sigma;

            q = 587; //q must be prime and congruent to 3 mod 8
            m1 = 13;
            m2 = 150; //m2 >= lambda*m1, such as lambda is the security parameter and lambda = ceil(1 + lg(q))

            r = ceil((double)(1 + (log(q)/log(3))));
            sigma = 1;

            lambda = ceil(1 + (log(q)/log(2)));

            if(q % 8 != 3) {
                cout << "q must be congruent to 3 mod 8.\n";
                return -1;
            }

            if(m1 < sigma || m1 < r) {
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

            HIBE hibe((double)q, m1, m2, k, sigma);
            hibe.Setup(10); // Setup algorithm with h = 10            
            
            break;
        }
        
        case 2: {
            /* Parameter set from (Roy et al., 2013). That depends on the cryptographic system requirements */
            RR sigma = to_RR(3.195); // Standard deviation
            int tailcut = 13;
            
            if(sigma*to_RR(tailcut) > power2_RR(sizeof(int)*8+1)-1) {
                cout << "Error! This distribution can not be simulated. Aborting..." << endl;
                return -1;
            }//end-if                
            
            Vec<int> ZigguratPoly, KnuthPoly;
            Samplers sampler(k);
            int nSamples = 1024; // #coefficients in the polynomial
            int nRectangles = 63; // Parameter of Ziggurat algorithm
            int omega = 107; // Parameter of Ziggurat algorithm
            int precision = 107;
            
            timestamp_t ts_start, ts_end;
            
            ts_start = get_timestamp();            
            ZigguratPoly = sampler.PolyGeneratorZiggurat(nSamples, nRectangles, sigma, omega, precision, tailcut); // Coefficients, rectangles, sigma, omega and precision
            ts_end = get_timestamp();            
                        
            cout << "[!] Ziggurat running time for " << nSamples << " samples: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;
            cout << ZigguratPoly << endl;
            
            ts_start = get_timestamp();                        
            KnuthPoly = sampler.PolyGeneratorKnuthYao(nSamples, precision, tailcut, sigma); // Coefficients, precision, tailcut, and sigma
            ts_end = get_timestamp();            
            
            cout << "[!] Knuth-Yao running time for " << nSamples << " samples: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;
            cout << KnuthPoly << endl;
            
            break;
        }
        case 3: {
            
            k = 3;//For testing
            Samplers sampler(k);

            int m = 8;

            double outD;
            int out, q;
            Vec< int > a, b;
            a.SetLength(m);
            b.SetLength(m);
            q = 4;
            
            for(int i = 0; i < m; i++) {
                a[i] = NTL::RandomBnd(q);
                b[i] = NTL::RandomBnd(q);
            }//end-for
            
            cout << "Polynomial a: ";
            for(int i = 0; i < m; i++)
                cout << a[i] << " ";
            cout << endl;
            cout << "Polynomial b: ";            
            for(int i = 0; i < m; i++)
                cout << b[i] << " ";
            cout << endl;
            
            sampler.InnerProduct(out, a, b);            
            cout << "<a, b>: " << out << endl;
            
            sampler.Norm(outD, a);
            cout << "Norm of a: " << outD << endl;
            sampler.Norm(outD, b);
            cout << "Norm of b: " << outD << endl;
            
            break;
        }
            
    }//end-switch
    
    return 0;
    
}

