/* 
 * File:   main.cpp
 * Author: jnortiz
 *
 * Created on March 10, 2015, 2:12 PM
 */

#include <cstdlib>
#include <math.h>
#include <NTL/ZZ.h>

#include "HIBE.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    double lambda, q;
    int m1, m2, k;
    
    q = 43; //q must be prime and congruent to 3 mod 8
    m1 = 7;
    m2 = 50; //m2 >= lambda*m1, such as lambda is the security parameter and lambda = ceil(1 + lg(q))
    k = 3; // n = 2^k is the degree of polynomials in R and R_0

    lambda = ceil(1 + (log(q)/log(2)));
    
    if((int)(q) % 8 != 3) {
        cout << "q must be congruent to 3 mod 8.\n";
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
    
    HIBE hibe(q, m1, m2, k);
    hibe.IdealTrapGen();
    
    return 0;
}

