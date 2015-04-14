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

int main(int argc, char** argv) {
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
//    hibe.Setup(10); // Setup algorithm with h = 10
    
    RR m, sigmaRR, n;
    ZZ omega;
    m = to_RR(63);
    sigmaRR = to_RR(32);
    n = to_RR(106);
    omega = 106;
    
    int nPositive = 0, nNegative = 0, nZero = 0;
    ZZ sample;
    
    // Following the reference implementation, omega = bit precision (n)    
    // Getting samples from the distribution    
    cout << "[*] Ziggurat sampling" << endl;
    for(int i = 0; i < 1000000; i++) {
        sample = hibe.ZiggutatO(m, sigmaRR, omega, n);
        if(sample > 0)
            nPositive++;
        else if(sample < 0)
            nNegative++;
        else
            nZero++;
    }//end-for
    
    cout << "\nPositive numbers: " << nPositive << endl;
    cout << "Negative numbers: " << nNegative << endl;
    cout << "Zero numbers: " << nZero << endl;
    
    return 0;
}

