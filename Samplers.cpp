/* 
 * File:   Samplers.cpp
 * Author: jnortiz
 * 
 * Created on April 24, 2015, 3:51 PM
 */

#include "Samplers.h"
#include <NTL/RR.h>
#include <NTL/matrix.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>

using namespace NTL;
using namespace std;

Samplers::Samplers() {
    
    int randint;
    int bytes_read;
    int fd = open("/dev/urandom", O_RDONLY);
    
    if (fd != -1) {
        bytes_read = read(fd, &randint, sizeof(randint));
        if (bytes_read != sizeof(randint))            
            fprintf(stderr, "read() failed (%d bytes read)\n", bytes_read);
    } else        
        fprintf(stderr, "open() failed\n");
    
    close(fd);
          
    NTL::SetSeed(to_ZZ(randint));
    
}//end-Samplers()

Samplers::~Samplers() {
}

/* It calculates the probability associated to x */
RR Samplers::Rho(RR sigma, RR x) {
    
    if(x == 0)
        return to_RR(1);
    
    return exp(-(power(x/sigma, 2))/2.0);
    
}//end-Rho()

/* It computes the probability associated with a sample x */
RR Samplers::Probability(RR x, RR sigma) {
    RR S = sigma*sqrt(2*ComputePi_RR());
    RR overS = 1/S;
    
    if(x == to_RR(0))
        return overS;
    
    return overS*exp(-(power(x/sigma, 2))/2.0);
    
}//end-Probability()

/* It selects between two given values depending on bit value. If the bit value 
 * is zero, the output becomes "a" */
int Select(int a, int b, unsigned bit) {
    unsigned mask;
    int output;
    mask = -bit;
    output = mask & (a ^ b);
    output = output ^ a;
    return output;
}//end-Select()

/* It produces a n-dimension polynomial with coefficient sampled from a 
 * Gaussian distribution using Ziggurat algorithm */
Vec<int> Samplers::PolyGeneratorZiggurat(int dimension, RR m, RR sigma, ZZ omega, RR n, RR tail) {
    
    cout << "[*] Ziggurat Gaussian sampling" << endl;

    /* Output precision setup */
    RR::SetOutputPrecision(to_long(n));
    
    Vec<int> polynomial;    
    polynomial.SetLength((long)dimension);

    // Only for statistical purposes
    int nIterations, nPositive, nNegative, nZero, samplesGen;
    nIterations = 0;
    nPositive = nNegative = nZero = 0;
    
    // Creating the rectangles partitioning
    this->DZCreatePartition(m, sigma, n, tail);
    
    // Following the reference implementation, omega = bit precision (n)    
    // Getting samples from the discrete distribution    
    do {
        samplesGen = 0;
        for(int i = 0; i < dimension; i++) {
            polynomial[i] = this->Ziggurat(m, sigma, omega);
            if(polynomial[i] <= to_int(tail*sigma))
                samplesGen++;
        }//end-for
        nIterations++;
    }while(samplesGen < dimension);

    cout << "[*] All samples were successfully generated in " << nIterations << " iteration(s)." << endl;
    
    return polynomial;
    
}//end-PolyGeneratorZiggurat()

/* Algorithm for sampling from discrete distribution using rejection method.
 * It does not avoid Rho function computation */
int Samplers::Ziggurat(RR m, RR sigma, ZZ omega) {
        
    ZZ b, s, x, powerOmega, yPrime, yBar;
    int aux, bit, curve, i, invalidSample, isZero, j, mInt, S;
    unsigned less, lessEqual, equal, greater;
        
    powerOmega = power2_ZZ(to_int(omega)); // 2^{\omega}
    mInt = to_int(m); // Number of rectangles
    invalidSample = 14*to_int(sigma);
    bit = 0;
    
    /* Estimated number of iterations required to obtain all requested samples 
     * in one iteration of PolyGeneratorZiggurat algorithm */
    for(int index = 0; index < 6; index++) { 
        i = RandomBnd(mInt) + 1; // Sample a random value in {0, ..., m-1} and add one to the result, e.g. select a rectangle in {1, ..., m}
        s = 1 - 2*RandomBits_long(1); // Sample a random signal s
        x = RandomBnd(this->X_ZZ[i] + 1); // Sample a x value between 0 and floor(x_i)
        b = RandomBits_long(1);
        yPrime = RandomBnd(powerOmega - 1);                 
        yBar = yPrime * (this->Y_ZZ[i-1] - this->Y_ZZ[i]);
        curve = to_int(to_RR(powerOmega) * (Rho(sigma, to_RR(x)) - this->Y[i]));
        
        isZero = to_int(-x);
    
        for(j = 0; j < (sizeof(int)*8)-1; j++) {
          aux = isZero >> 1;
          isZero = isZero | aux;
        }//end-for
        
        greater = (isZero+1)%2; //Certainly x > 0 if it's not zero
        
        less = to_int(x) - this->X_ZZ[i-1]; //If x > X_ZZ[i] then "less" is a negative number
        equal = less - 1;
        
        for(j = 0; j < sizeof(int)*8-1; j++)
            less >>= 1;
        
        for(j = 0; j < sizeof(int)*8-1; j++) {
            aux = equal >> 1;            
            equal = equal & aux;
        }//end-for
        
        lessEqual = less | equal;
        
        // First case: The sampled x is in the left side of the i-th rectangle; e.g., 0 < x <= this->X_ZZ[i-1]        
        bit = (bit | (greater & lessEqual));
        
        // Second case: If x = 0, define s*x as a sample with probability of 50%        
        bit = (bit | (isZero & (b+1)%2));
                        
        less = to_int(x) - curve;
        equal = less - 1;
        
        for(j = 0; j < sizeof(int)*8-1; j++)
            less >>= 1;
        
        for(j = 0; j < sizeof(int)*8-1; j++) {
            aux = equal >> 1;            
            equal = equal & aux; //If "equal" becomes 1, so x = 0
        }//end-for        
           
        // Third case: the sampled x is below to the curve
        bit = (bit | (less | equal));
        
        /* If the bit becomes 1, the valid sample s*x is assigned to S. 
         * The bit is an or operation between the bit value of the last two iterations. 
         * It prevents a valid sample to be overwritten. */
        S = Select(invalidSample, to_int(s*x), bit); 

    }//end-for
    
    return S;    
    
}//end-ZigguratO()

/* DZCreatePartition defines the x and y axes of rectangles in the Gaussian distribution */
void Samplers::DZCreatePartition(RR m, RR sigma, RR n, RR tail) {
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
            
    tailcut = tail*sigma;
    
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
            this->X_ZZ[i] = to_int(this->X[i]);
            this->Y_ZZ[i] = to_int(this->Y[i]);
        }//end-for
                
    }
    else // No valid partition was found
        cout << "[*] DZCreatePartition status: Error!" << endl;
        
    
}//end-DZCreatePartition()

/* It is used in DZCreatePartition to define the distance y0 */
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
    //    X[nRect] = to_RR(tail*sigma);
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

/* Using the Knuth-Yao algorithm, it produces a n-dimension polynomial 
 * with coefficients from the Gaussian distribution */
Vec<int> Samplers::PolyGeneratorKnuthYao(int dimension, int precision, int tailcut, RR sigma) {
    
    cout << "\n[*] Knuth-Yao Gaussian sampling" << endl;
    
    /* Output precision setup */
    RR::SetOutputPrecision(to_long(precision));
    
    Vec<int> polynomial;
    int bound, samplesGen, iterations;
    
    polynomial.SetLength((long)dimension);
    bound = tailcut*to_int(sigma);
    iterations = 0;
    
    do {
        samplesGen = 0; // It counts the number of successfully generated samples
        for(int i = 0; i < dimension; i++) {
            polynomial.put(i, this->KnuthYao(precision, tailcut, sigma));
            // Samples equal to the bound won't be accepted
            if(polynomial.get(i) < bound && polynomial.get(i) > -bound)
                samplesGen++;
        }//end-for
        iterations++;
    }while(samplesGen < dimension);
    
    if(samplesGen == dimension)
        cout << "[*] All samples were successfully generated in " 
                << iterations << " iteration(s)." << endl;

    return polynomial;
    
}//end-PolyGeneratorKnuthYao()

/* Knuth-Yao algorithm to obtain a sample from the discrete Gaussian */
int Samplers::KnuthYao(int precision, int tailcut, RR sigma) {
    
    if(this->P.NumRows() != precision) {
        this->BuildProbabilityMatrix(precision, tailcut, sigma);
        cout << "[*] Probability matrix building status: Pass!" << endl;
    }//end-if
    
    ZZ r;
    int bound, col, d, i, invalidSample, searchRange, S;
    unsigned aux, enable, hit;
    
    bound = tailcut*to_int(sigma);
    d = 0; //Distance
    hit = 0;
    invalidSample = 3*bound;
    /* Search range required to obtain all samples with only one iteration 
     * in PolyGeneratorKnuthYao() algorithm */
    searchRange = this->P.NumRows()/4;
    S = 0;
    
    for(int row = 0; row < searchRange; row++) {
        
        r = RandomBits_long(1); // Random choice between 0 and 1
        d = 2*d + to_int(r); // Distance calculus
        
        for(col = 0; col < this->P.NumCols(); col++) {
            
            d = d - to_int(this->P[row][col]);
            
            // Enable turns one just in case d is equal to -1
            enable = (unsigned)d;
            for(i = 0; i < (sizeof(int)*8)-1; i++) {
              aux = enable >> 1;
              enable = enable & aux;
            }//end-for
            
            /* When enable&!hit becomes 1, "col" is added to "S";
             * e.g. enable = 1 and hit = 0 */
            S += Select(invalidSample, col, (enable & !hit));
            hit += (enable & !hit);
                            
        }//end-for
    }//end-for
    
    /* Note: the "col" value is in [0, 2*bound]. So, the invalid sample must be 
     * greater than 2*bound. */
    S = S % invalidSample;
    S -= bound;
    
    return S;
    
}//end-Knuth-Yao()

/* This method build the probability matrix for samples in the range 
 * given by [-tailcut*\floor(sigma), +tailcut*\floor(sigma)] */
void Samplers::BuildProbabilityMatrix(int precision, int tailcut, RR sigma) {
    
    // The random variable consists of elements in [-tailcut*sigma, tailcut*sigma]
    mat_ZZ aux_P;
    int i, bound;
    RR probOfX;
    
    i = 0;
    bound = tailcut*to_int(sigma);
    
    aux_P.SetDims(2*bound+1, precision);
    this->P.SetDims(precision, 2*bound+1);
    
    for(int x = -bound; x <= bound; x++, i++) {
        probOfX = Probability(to_RR(x), sigma);
        BinaryExpansion(aux_P, probOfX, precision, i);
    }//end-for
    
    // Changing the elements positioning in P to decrease page fault in future access
    int row_aux_P = 0;
    for(int col = this->P.NumCols()-1; col >= 0; col--) {
        for(int row = 0; row < this->P.NumRows(); row++)
            this->P.put(row, col, aux_P.get(row_aux_P, row));
        row_aux_P++;        
    }//end-for
    
//    this->PrintMatrix("Probability matrix", this->P);
    
}//end-BuildProbabilityMatrix()

/* Method for computing the binary expansion of a given probability in [0, 1] */
void Samplers::BinaryExpansion(mat_ZZ& aux_P, RR probability, int precision, int index) {
        
    RR pow;
    int i, j;
    i = -1;
    j = 0;
    
    while(probability > 0 && j < precision) {
        pow = power2_RR(i--); // 2^{i}
        if(pow <= probability) {
            aux_P.put(index, j, to_ZZ(1));
            probability -= pow;
        } else
            aux_P.put(index, j, to_ZZ(0));
        j++;
    }//end-while
        
}//end-BinaryExpansion()

void Samplers::PrintMatrix(const string& name, const mat_ZZ& matrix) {
    
    cout << "\n/** " << name << " **/" << endl;
    for(int i = 0; i < matrix.NumRows(); i++) {
        for(int j = 0; j < matrix.NumCols(); j++)
            cout << matrix[i][j] << " ";
        cout << endl;
    }//end-for
    
}//end-PrintVectorZZX()