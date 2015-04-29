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

int main(void) {
        
    int action = 2;
    
    switch(action) {
        case 1: { // Setup algorithm from HIBE scheme
            double lambda;
            int q, m1, m2, k, r, sigma;

            q = 587; //q must be prime and congruent to 3 mod 8
            m1 = 13;
            m2 = 150; //m2 >= lambda*m1, such as lambda is the security parameter and lambda = ceil(1 + lg(q))
            k = 8; // n = 2^k is the degree of polynomials in R and R_0

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
            Samplers sampler;
            int nSamples = 512; // #coefficients in the polynomial
            RR nRectangles = to_RR(63); // Parameter of Ziggurat algorithm
            RR sigma = to_RR(3.195); // Standard deviation
            ZZ omega = to_ZZ(107); // Parameter of Ziggurat algorithm
            RR precision = to_RR(107);
            int tailcut = 13;    
            
            cout << sampler.PolyGeneratorZiggurat(nSamples, nRectangles, sigma, omega, precision, to_RR(tailcut)); // Coefficients, rectangles, sigma, omega and precision
            cout << endl;
            cout << sampler.PolyGeneratorKnuthYao(nSamples, to_int(precision), tailcut, sigma); // Coefficients, precision, tailcut, and sigma
            break;
        }            
    }//end-switch
    
    return 0;
    
}

