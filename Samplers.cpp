/* 
 * File:   Samplers.cpp
 * Author: jnortiz
 * 
 * Created on April 24, 2015, 3:51 PM
 */

#include "Samplers.h"

#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include <NTL/RR.h>

#define CONSTANT_TIME 1

using namespace NTL;
using namespace std;

Samplers::Samplers() {
}

Samplers::Samplers(const Samplers& orig) {
}

Samplers::~Samplers() {
}

// The most expensive function in Ziggurat algorithm. It calculates the probability associated to x
RR Samplers::Rho(RR sigma, RR x) {
    
    if(x == 0)
        return to_RR(1);
    
    return exp(-(power(x/sigma, 2))/2.0);
    
}//end-Rho()

// Function used in ZigguratO algorithm to decide if x is up or down to the curve (Rho function)
ZZ Samplers::sLine(ZZ x0, ZZ x1, ZZ y0, ZZ y1, ZZ x, long int i) {
    
    // Consider x0 = x_{i-1}, x1 =s x_i, y0 = \bar{y}_{i-1} and y1 = \bar{y}_i
    if(x1 == x0)
        return to_ZZ(-1);
    
    ZZ y0Hat, y1Hat;    
    y0Hat = 0; // Just in case i is neither greater than 1 nor equal to 1
    y1Hat = to_int(y1);
    
    if(i > 1)
        y0Hat = y0;
    if(i == 1)
        y0Hat = to_ZZ(1);
    
    return (y1Hat - y0Hat)/(x1 - x0)*(x - x1);
    
}//end-sLine()

// Sampling algorithm with optimization to avoid Rho function computation
ZZ Samplers::ZiggutatO(RR m, RR sigma, ZZ omega) {
        
    ZZ curve, sigma_ZZ, s, sample, x, powerOmega, yPrime, yBar;
    long int i, mInt;
        
    powerOmega = power2_ZZ(to_int(omega)); // 2^{\omega}
    yPrime = RandomBnd(powerOmega - 1);                
    
    mInt = to_int(m);
    sigma_ZZ = to_ZZ(sigma);
    unsigned bit = 0;    
    sample = 14*sigma_ZZ;
    ZZ b;
    
#ifdef CONSTANT_TIME
    
    int index;
    for(index = 0; index < 3; index++) { // Estimated number of iterations required to obtain all requested samples
        i = RandomBnd(mInt) + 1; // Sample a random value in {0, ..., m-1} and add one to the result, e.g. select a rectangle in {1, ..., m}
        s = 1 - 2 * RandomBits_long(1); // Sample a random signal s
        x = RandomBnd(this->X_ZZ[i] + 1); // Sample a x value between 0 and floor(x_i)
        b = RandomBits_long(1);
        yBar = yPrime * (this->Y_ZZ[i-1] - this->Y_ZZ[i]);
        curve = powerOmega * to_ZZ(Rho(sigma, to_RR(x)) - Y[i]);
        
        // The sampled x is in the left side of the i-th rectangle
        bit = (x > to_ZZ(0) && x <= this->X_ZZ[i-1])? 1 : 0;
        
        // If x = 0, define s*x as a sample with probability of 50%
        bit = (x == to_ZZ(0) && b == to_ZZ(0))? 1 : 0;
                
        // If x is in the right size of the i-th rectangle
        bit = (this->X_ZZ[i] + 1 <= sigma_ZZ // In concave-down case
            && (yBar <= powerOmega * sLine(this->X_ZZ[i-1], this->X_ZZ[i], this->Y_ZZ[i-1],this->Y_ZZ[i], x, i)
            ||  yBar <= curve ))? 1 : 0;        
        bit = (sigma_ZZ <= this->X_ZZ[i-1] // In concave-up case
            && yBar < powerOmega * sLine(this->X_ZZ[i-1], this->X_ZZ[i], this->Y_ZZ[i-1], this->Y_ZZ[i], (x-1), i) 
            && yBar < curve)? 1 : 0;
        bit = (yBar <= to_ZZ(to_RR(powerOmega) * (Rho(sigma, to_RR(x)) - this->Y[i])))? 1 : 0;
        
        if(bit)            
            sample = s*x;
        
    }//end-for
    
#else
    
    while(!bit) {        
        i = RandomBnd(mInt) + 1; // Sample a random value in {0, ..., m-1} and add one to the result, e.g. select a rectangle in {1, ..., m}
        s = 1 - 2 * RandomBits_long(1); // Sample a random signal s
        x = RandomBnd(this->X_ZZ[i] + 1); // Sample a x value between 0 and floor(x_i)
        
        if(x > to_ZZ(0) && x <= this->X_ZZ[i-1]) { // The sampled x is in the left side of the i-th rectangle
            sample = s*x;
            bit = 1;
        }//end-if
        else
            if(x == to_ZZ(0)) {
                b = RandomBits_long(1); // It selects between 0 or 1
                if(b == to_ZZ(0)) {
                    sample = s*x;
                    bit = 1;
                }//end-if
            }//end-if
            else {
                yPrime = RandomBnd(powerOmega - 1);                
                yBar = yPrime * (this->Y_ZZ[i-1] - this->Y_ZZ[i]);
                curve = (powerOmega * to_ZZ(Rho(sigma, to_RR(x)) - Y[i]));
                
                if(this->X_ZZ[i] + 1 <= sigma_ZZ) // In concave-down case
                    if(yBar <= powerOmega * sLine(this->X_ZZ[i-1], this->X_ZZ[i], this->Y_ZZ[i-1],this->Y_ZZ[i], x, i) || yBar <= curve) {
                        sample = s*x;
                        bit = 1;
                    }//end-if
                else
                    if(sigma_ZZ <= this->X_ZZ[i-1]) // In concave-up case
                        if(yBar < powerOmega * sLine(this->X_ZZ[i-1], this->X_ZZ[i], this->Y_ZZ[i-1], this->Y_ZZ[i], (x-1), i) && yBar < curve) {
                            sample = s*x;
                            bit = 1;
                        }//end-if
                    else
                        if(yBar <= to_ZZ(to_RR(powerOmega) * (Rho(sigma, to_RR(x)) - this->Y[i]))) {
                           sample = s*x; 
                           bit = 1;
                        }//end-if
            }//end-else
    }//end-while
    
#endif    
    
    return sample;    
    
}//end-ZigguratO()

// DZCreatePartition defines the x and y axes of rectangles in the Gaussian distribution
void Samplers::DZCreatePartition(RR m, RR sigma, RR n) {
    /*
     * Parameters description:
     * m: number of rectangles
     * sigma: Gaussian distribution parameter
     * n: bit precision
     */
    
    // The approximation error must be less than 2^{-n}
    RR statDistance = power2_RR(-to_long(n));
    
    RR c, tailcut, y0;
    long nRect = to_long(m);
    bool validPartition = 0;
    
    this->X.SetLength(nRect + 1);
    this->Y.SetLength(nRect + 1);
            
    tailcut = to_RR(13) * sigma; // We assume that tailcut = 13
    
    while(tailcut < to_RR(14)*sigma) {
        
        c = to_RR(1);        
        y0 = to_RR(0);
        this->X[nRect] = tailcut;        
        
        while(y0 < 1 || abs(y0)-1 > statDistance) {
            
            y0 = DZRecursion(m, c, sigma);
            
            if(y0 == -1) // Error in "inv" calculus
                break;
            
            if(y0 >= 1) { // The found partitioning is valid
                validPartition = 1;
                break;
            }//end-if
            
            c++; // If no partitioning is found, repeat the process with incremented constant
            
        }//end-while
        
        if(validPartition) // We return the first valid partition found
            break;
        
        if(y0 < 1 && abs(y0)-1 > statDistance)
            tailcut++;
        else
            break;
        
    }//end-while
    
    if(validPartition) {
        cout << "[*] DZCreatePartition status: Pass!" << endl;
        
        this->X_ZZ.SetLength(this->X.length());
        this->Y_ZZ.SetLength(this->Y.length());

        long int i;
        // Computing the floor and ZZ format of values in X and Y vectors
        for(i = 0; i < this->X_ZZ.length(); i++) { // X and Y vectors have the same length m
            this->X_ZZ[i] = to_ZZ(this->X[i]);
            this->Y_ZZ[i] = to_ZZ(this->Y[i]);
        }//end-for
                
        cout << "Final distance: " << y0 << endl;
    }
    else // No valid partition was found
        cout << "[*] DZCreatePartition status: Error!" << endl;
        
    
}//end-DZCreatePartition()

// Used in DZCreatePartition to define the distance y0
RR Samplers::DZRecursion(RR m, RR c, RR sigma) {
    
    RR inv, minus2, overM, over2, S;
    int nRect;
    
    nRect = to_int(m);
    minus2 = to_RR(-2);
    overM = to_RR(1)/m;
    over2 = to_RR(1)/to_RR(2);
    
    if(inv < to_RR(0))
        return to_RR(-1);
    
    // "Size" of each rectangle
    S = sigma * overM * sqrt(ComputePi_RR() * over2) * to_RR(c);
    
    /* 
     * In the reference code, this statement is always executed. Although, it overwrites the 
     * "this->X[nRect] = tailcut;" statement in DZCreatePartition function.
     */
    //    X[nRect] = to_RR(13*sigma);
    this->Y[nRect] = Rho(sigma, this->X[nRect]);
    
    inv = minus2 * log(S/to_RR(TruncToZZ(this->X[nRect]) + to_ZZ(1)));
    
    this->X[nRect-1] = sigma * sqrt(inv);
    this->Y[nRect-1] = Rho(sigma, this->X[nRect-1]);
    
    for(int i = nRect-2; i > 0; i--) {
        inv = minus2 * log(S/to_RR(TruncToZZ(this->X[i+1]) + to_ZZ(1))) + Rho(sigma, this->X[i]);
        
        if(inv < to_RR(0))
            return to_RR(-1);
        
        this->X[i] = sigma * sqrt(inv); // Rho(sigma, x_i) = y_i
        this->Y[i] = exp(-over2 * power(this->X[i]/sigma, 2)); // Rho^{-1}(sigma, Rho(sigma, x_i)) = x_i               
    }//end-for
    
    this->Y[0] = (S/(to_RR(1) + this->X[1])) + Rho(sigma, this->X[1]);
        
    return this->Y[0];
    
}//end-DZCreatePartition()

Vec<ZZ> Samplers::PolyGeneratorZiggurat(int dimension, RR m, RR sigma, ZZ omega, RR n) {
    
    cout << "[*] Ziggurat Gaussian sampling" << endl;

    /* Output precision setup */
    RR::SetOutputPrecision(to_long(n));
    
    Vec<ZZ> polynomial;    
    polynomial.SetLength((long)dimension);

    // Only for statistical purposes
    int nPositive, nNegative, nZero;
    nPositive = nNegative = nZero = 0;
    
    // Creating the rectangles partitioning
    this->DZCreatePartition(m, sigma, n);
    
    // Following the reference implementation, omega = bit precision (n)    
    // Getting samples from the distribution    
    for(int i = 0; i < dimension; i++) {
        polynomial[i] = this->ZiggutatO(m, sigma, omega);
        if(polynomial[i] <= 13*to_ZZ(sigma)) { // Tailcut t = 13
            if(polynomial[i] > 0)
                nPositive++;
            else if(polynomial[i] < 0)
                nNegative++;
            else
                nZero++;
        }//end-if
    }//end-for

    cout << "\nPositive numbers: " << nPositive << endl;
    cout << "Negative numbers: " << nNegative << endl;
    cout << "Zero numbers: " << nZero << endl;

    int samplesGen = nPositive + nNegative + nZero;

    if(samplesGen == dimension)
        cout << "All samples were successfully generated." << endl;
    else
        cout << "Correctness percentage: " << (float)(100*samplesGen)/dimension << endl;
    
    return polynomial;
    
}//end-PolyGeneratorZiggurat()
