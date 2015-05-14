/* 
 * File:   main.cpp
 * Author: jnortiz
 *
 * Created on March 10, 2015, 2:12 PM
 */

#include <cstdlib>
#include <math.h>
#include "HIBE.h"

using namespace std;

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp() {
    struct timespec now;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &now);
    return now.tv_nsec + (timestamp_t)now.tv_sec * 1000000000.0;
}

int main(int argc, char** argv) {
    
    int q, m1, m2, k, sigma;

    q = 587; //q must be prime and congruent to 3 mod 8
    m1 = 13;
    m2 = 150; //m2 >= lambda*m1, such as lambda is the security parameter and lambda = ceil(1 + lg(q))
    k = 8; // n = 2^k is the degree of polynomials in R and R_0
    sigma = 1;

    HIBE hibe((double)q, m1, m2, k, sigma);
    
    int option = 2;
    
    switch(option) {
        
        case 1: {

            int r;
            double lambda;
    
            r = ceil((double)(1 + (log(q)/log(3))));
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

            hibe.Setup(10); // Setup algorithm with h = 10
            
            break;
        }
        
        case 2: {

            Vec<int> ZigguratPoly;
            int nSamples = 1024; // #coefficients in the polynomial
            RR nRectangles = to_RR(63); // Parameter of Ziggurat algorithm
            RR sigma = to_RR(3.195); // Standard deviation
            ZZ omega = to_ZZ(107); // Parameter of Ziggurat algorithm
            RR precision = to_RR(107);
            int tailcut = 13;    
            
            timestamp_t ts_start, ts_end;
            
            ts_start = get_timestamp();            
            ZigguratPoly = hibe.PolyGeneratorZigguratO(nSamples, nRectangles, sigma, omega, precision, to_RR(tailcut)); // Coefficients, rectangles, sigma, omega and precision
            ts_end = get_timestamp();            
                        
            cout << "[!] Ziggurat running time for " << nSamples << " samples: " << (float)((ts_end - ts_start)/1000000000.0) << " s." << endl;
//            cout << ZigguratPoly << endl;            
            
            break;
        }
    
    }
    return 0;
}

