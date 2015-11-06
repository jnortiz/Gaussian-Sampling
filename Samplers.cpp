/* 
 * File:   Samplers.cpp
 * Author: jnortiz
 * 
 * Created on April 24, 2015, 3:51 PM
 */

#include "Samplers.h"
#include <NTL/ZZ_pEX.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>

using namespace NTL;
using namespace std;

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp() {
    struct timespec now;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &now);
    return now.tv_nsec + (timestamp_t)now.tv_sec * 1000000000.0;
}

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

Samplers::~Samplers() {};

/* Probability associated with x */
RR Samplers::Rho(RR sigma, RR x) {
    
    if(x == 0)
        return to_RR(1);
    
    return exp(-(power(x/sigma, 2))/2.0);
    
}//end-Rho()

RR Samplers::Probability(RR x, RR sigma, RR c) {
    RR S = sigma*sqrt(2*ComputePi_RR());
    RR overS = 1/S;
    
    if(x == to_RR(0))
        return overS;
    
    return overS*exp(-(power((x-c)/sigma, 2))/2.0);
    
}//end-Probability()

/* It selects between two given values depending on a bit. 
 * If the bit is zero, the output becomes "a" */
int Select(int a, int b, unsigned bit) {
    unsigned mask;
    int output;
    mask = -bit;
    output = mask & (a ^ b);
    output = output ^ a;
    return output;
}//end-Select()

/* Testing if x = 0 */
int isZero(int x) {
    unsigned zero;
    zero = x;
    zero = 1 ^ ((zero | -zero) >> 31) & 1;
    return zero;    
}//end-isZero()

/* Testing if x is less than y */
unsigned lessThan(int x, int y) {    
    unsigned less;    
    less = x-y;
    less >>= sizeof(int)*8-1;    
    return less;        
}//end-lessThan()

/* Testing if x = y */
unsigned isEqual(int x, int y) {    
    unsigned equal;    
    equal = x-y; // "equal" turns 0 if x = y    
    equal = 1 ^ ((equal | -equal) >> 31) & 1; // "equal" turns 1 iff enable was 0
    return equal;    
}//end-isEqual()

/* It produces a n-dimension polynomial with coefficient sampled from a 
 * Gaussian distribution using Ziggurat algorithm */
Vec<int> Samplers::PolyGeneratorZiggurat(int dimension, int m, RR sigma, int omega, int n, int tail) {
    
    cout << "\n[*] Ziggurat Gaussian sampling" << endl;
    
    Vec<int> polynomial;    
    polynomial.SetLength((long)dimension);

    // Only for statistical purposes
    int bound, samplesGen;
    bound = to_int(tail*sigma);
    
    // Creating the rectangles partitioning
    this->DZCreatePartition(m, sigma, n, tail);
        
    // Following the reference implementation, omega = bit precision (n)    
    // Getting samples from the discrete distribution    
    do {
        samplesGen = 0;
        for(int i = 0; i < dimension; i++) {
            polynomial[i] = this->Ziggurat(m, n, sigma, omega);
            if(polynomial[i] <= bound)
                samplesGen++;
        }//end-for
    }while(samplesGen < dimension);

    return polynomial;
    
}//end-PolyGeneratorZiggurat()

/* Algorithm for sampling from discrete distribution using rejection method.
 * It does not avoid Rho function computation */
int Samplers::Ziggurat(int m, int n, RR sigma, int omega) {
        
    // Precision of floating point operations
    RR::SetPrecision(n);    

    // See below to explanation about this comments
    ZZ powerOmega;
    int curve;
    int bit, i, invalidSample, zero, s, S, x;
    unsigned less, lessEqual, equal, greater;
    long b;
        
    powerOmega = power2_ZZ(omega); // 2^{\omega}
    invalidSample = 14*to_int(sigma); // Value out of the range
    bit = 0;
    
    /* Estimated number of iterations required to obtain all requested samples 
     * in one iteration of PolyGeneratorZiggurat algorithm */
    for(int index = 0; index < 6; index++) { 
        i = RandomBnd(m) + 1; // Sample a random value in {0, ..., m-1} and add one to the result, e.g. select a rectangle in {1, ..., m}
        s = 1 - 2*RandomBits_long(1); // Sample a random signal s
        x = RandomBnd(this->X_ZZ[i] + 1); // Sample a x value between 0 and floor(x_i)
        b = RandomBits_long(1);
        
        zero = isZero(x);               
        greater = (zero+1)%2; //Certainly x > 0 if it's not zero        
        less = lessThan(x, this->X_ZZ[i-1]);
        equal = isEqual(x, this->X_ZZ[i-1]);                
        lessEqual = less | equal;
        
        /* First case: The sampled x is in the left side of the i-th rectangle; 
         * e.g., 0 < x <= this->X_ZZ[i-1] */
        bit = (bit | (greater & lessEqual));
        
        /* Second case: If x = 0, define s*x as a sample with probability of 50% */
        bit = (bit | (zero & (b+1)%2));
                        
        // Suppression of 3rd case
        curve = to_int(to_RR(powerOmega) * (Rho(sigma, to_RR(x)) - this->Y[i]));
        less = lessThan(x, curve);
        equal = isEqual(x, curve);
        
        // Third case: the sampled x is below to the curve and in the left rectangle
        bit = (bit | (less | equal)); 
        
        /* If the bit becomes 1, the valid sample s*x is assigned to S. 
         * The bit is an OR operation between the last and the current bit value. 
         * It prevents a valid sample to be overwritten. */
        S = Select(invalidSample, s*x, bit); 

    }//end-for
    
    return S;    
    
}//end-Ziggurat()

/* DZCreatePartition defines the x and y axes of rectangles in the Gaussian distribution */
void Samplers::DZCreatePartition(int m, RR sigma, int n, int tail) {
    /*
     * Parameters description:
     * m: number of rectangles
     * sigma: Gaussian distribution parameter
     * n: bit precision
     */
    
    /* Output precision setup */
    RR::SetPrecision((long)n);
    
    // The approximation error must be less than 2^{-n}
    RR statDistance = power2_RR(-((long)n));
    
    Vec<RR> bestX, X, Y;
    RR bestdiff, c, cc, cl, cu, lastdiff, mRR, tailcut, y0;
    int i;
    
    bestX.SetLength(m + 1);
    X.SetLength(m + 1);
    Y.SetLength(m + 1);
            
    bestX[m] = -1;
    mRR = to_RR(m);        
    
    tailcut = to_RR(tail)*sigma;
    c = 1 + 1/mRR;
    bestdiff = to_RR(3);
    
    while(tailcut < to_RR(tail+1)*sigma) {
        
        cu = to_RR(0);
        cl = to_RR(1);
        y0 = to_RR(-1);
        lastdiff = to_RR(-2);
        X[m] = tailcut;        
        
        while(y0 < 0 || abs(y0) > statDistance && abs(y0 - lastdiff) > statDistance) {
            
            cc = c;
            lastdiff = y0;
            
            y0 = DZRecursion(X, Y, m, tail, c, sigma) - to_RR(1);            
            
            if(y0 == -2) // Error in "inv" calculus
                break;
            
            if(y0 >= 0) { // The found partitioning is valid

                for(int i = 0; i < m+1; i++)
                    bestX[i] = X[i];
                
                cc = c;
                cu = c;
                bestdiff = y0;                
                
            } else
                cl = c;
            
            if(cu < cl)
                c += 1/mRR;
            else
                c = (cu + cl)/to_RR(2);
            
            if(c >= 11 || y0 == lastdiff)
                break;
            
        }//end-while
                
        if(y0 < 0 || abs(y0) > statDistance && abs(y0 - lastdiff) > statDistance)
            tailcut++;
        else
            break;
        
    }//end-while
    
    if(bestX[m] != -1) {
        cout << "[*] DZCreatePartition status: Pass!" << endl;
        
        this->X.SetLength(m + 1);
        this->Y.SetLength(m + 1);
        this->X_ZZ.SetLength(m + 1);
        this->Y_ZZ.SetLength(m + 1);

        // Computing the floor and ZZ format of values in X and Y vectors
        for(i = 0; i < this->X_ZZ.length(); i++) { // X and Y vectors have the same length m
            this->X[i] = bestX[i];
            this->Y[i] = Y[i];
            this->X_ZZ[i] = to_int(this->X[i]);
            this->Y_ZZ[i] = to_int(this->Y[i]);
        }//end-for        
        
    }//end-if
    else // No valid partition was found
        cout << "[*] DZCreatePartition status: Error!" << endl;
            
}//end-DZCreatePartition()

/* It is used in DZCreatePartition to define the distance y0 */
RR Samplers::DZRecursion(Vec<RR>& X, Vec<RR>& Y, int m, int tail, RR c, RR sigma) {
    
    RR inv, minus2, overM, over2, S;
    
    minus2 = to_RR(-2);
    overM = to_RR(1)/to_RR(m);
    over2 = to_RR(1)/to_RR(2);
        
    // "Size" of each rectangle
    S = sigma * overM * sqrt(ComputePi_RR()/to_RR(2)) * c;
    
    X[m] = to_RR(tail)*sigma;
    Y[m] = this->Rho(sigma, to_RR(TruncToZZ(X[m]) + to_ZZ(1)));
    
    inv = minus2 * log(S/to_RR(TruncToZZ(X[m]) + to_ZZ(1)));
    
    if(inv < to_RR(0))
        return to_RR(-1);
    
    X[m-1] = sigma * sqrt(inv);
    Y[m-1] = this->Rho(sigma, X[m-1]);
    
    for(int i = m-2; i > 0; i--) {
        
        inv = minus2 * log(S/to_RR(TruncToZZ(X[i+1]) + to_ZZ(1)) + Rho(sigma, X[i+1]));
        
        if(inv < to_RR(0))
            return to_RR(-1);
        
        X[i] = sigma * sqrt(inv); // Rho(sigma, x_i) = y_i
        Y[i] = exp(-over2 * power(X[i]/sigma, 2)); // Rho^{-1}(sigma, Rho(sigma, x_i)) = x_i               
        
    }//end-for
    
    Y[0] = (S/to_RR(to_ZZ(1) + TruncToZZ(X[1]))) + this->Rho(sigma, X[1]);
    
    return Y[0];
    
}//end-DZCreatePartition()

/* Using the Knuth-Yao algorithm, it produces a n-dimension polynomial 
 * with coefficients from the Gaussian distribution */
Vec<int> Samplers::PolyGeneratorKnuthYao(int dimension, int precision, float tailcut, RR sigma, RR c) {
        
    cout << "\n[*] Knuth-Yao Gaussian sampling" << endl;
        
    this->BuildProbabilityMatrix(precision, tailcut, sigma, c);
    
    Vec<int> polynomial;
    int bound, center, sample;
    
    polynomial.SetLength((long)dimension);
    bound = ((int)tailcut)*to_int(sigma);
    center = to_int(c);
    
    for(int i = 0; i < dimension; i++) {
        do{
            sample = this->KnuthYao(tailcut, sigma, c);
        } while(sample >= (center + bound) || sample <= (center - bound));                
        polynomial.put(i, sample);
    }//end-for
    
    return polynomial;
    
}//end-PolyGeneratorKnuthYao()

/* Knuth-Yao algorithm to obtain a sample from the discrete Gaussian */
int Samplers::KnuthYao(float tailcut, RR sigma, RR c) {

    int bound, center, col, d, invalidSample, pNumRows, pNumCols, S, signal;
    unsigned enable, hit;
    unsigned long r;
    
    bound = ((int)tailcut)*to_int(sigma);
    center = to_int(c);
    d = 0; //Distance
    hit = 0;
    signal = 1 - 2*RandomBits_long(1); // Sample a random signal s    
    invalidSample = bound+1;
    pNumRows = this->P.length(); // Precision
    pNumCols = this->P[0].length();    
    
    Vec<int> randomBits;
    randomBits.SetLength(pNumRows);
    
    int i, index, j, length;
    length = sizeof(unsigned long)*8; // 64 bits 
    
    index = 0;
    for(i = 0; i < (pNumRows/length+1); i++) {
        r = RandomWord(); // It returns a word filled with pseudo-random bits
        for(j = 0; j < length, index < pNumRows; j++, r >>= 1)
            randomBits[index++] = (r & 1); // Getting the least significant bit
    }//end-for
    
    S = 0;
    
    for(int row = 0; row < pNumRows; row++) {
        
        d = 2*d + randomBits[row]; // Distance calculus
        
        for(col = this->begin[row]; col < pNumCols; col++) {
            
            d = d - this->P[row][col];
            
            enable = (unsigned)(d + 1); // "enable" turns 0 iff d = -1
            enable = 1 ^ ((enable | -enable) >> 31) & 1; // "enable" turns 1 iff "enable" was 0
             
            /* When enable&!hit becomes 1, "col" is added to "S";
             * e.g. enable = 1 and hit = 0 */
            S += Select(invalidSample, col, (enable & !hit));
            hit += (enable & !hit);
                            
        }//end-for
        
    }//end-for
    
    /* Note: the "col" value is in [0, bound]. So, the invalid sample must be 
     * greater than bound. */
    S %= invalidSample;
    S = S - bound + center;
    S *= signal;
    
    return S;
    
}//end-Knuth-Yao()

/* This method build the probability matrix for samples in the range 
 * [-tailcut*\floor(sigma), +tailcut*\floor(sigma)] */
void Samplers::BuildProbabilityMatrix(int precision, float tailcut, RR sigma, RR c) {
    
    RR::SetPrecision(to_long(precision));

    Vec< Vec<int> > auxP;
    Vec<int> auxBegin;
    
    // The random variable consists of elements in [c-tailcut*sigma, c+tailcut*sigma]
    int i, j, bound, pNumCols, pNumRows, x;
    vec_RR probOfX;
    RR pow;
    
    bound = ((int)tailcut)*to_int(sigma);
    
    probOfX.SetLength(bound+1);
       
    auxP.SetLength(precision);
    for(i = 0; i < auxP.length(); i++)
        auxP[i].SetLength(bound+1);

    for(x = bound; x > 0; x--)
        probOfX[bound-x] = Probability(to_RR(x) + c, sigma, c);
    div(probOfX[bound], Probability(to_RR(0) + c, sigma, c), to_RR(2));
    
    i = -1;
    for(j = 0; j < precision; j++) {
        pow = power2_RR(i--); // 2^{i}
        for(x = bound; x >= 0; x--) {
            auxP[j][bound-x] = 0;                
            if(probOfX[bound-x] >= pow) {
                auxP[j][bound-x] = 1;
                probOfX[bound-x] -= pow;
            }//end-if
        }//end-for
    }//end-while
    
    this->P = auxP;
    
    // Uncomment this line if you want to preview the probability matrix P
//    this->PrintMatrix("Probability matrix", this->P);
    
    pNumCols = this->P[0].length();
    pNumRows = this->P.length();
    
    auxBegin.SetLength(pNumRows);
    
    // Computing in which position the non-zero values in P start and end 
    for(i = 0; i < pNumRows; i++) {
        
        auxBegin[i] = pNumCols-1;
        
        for(j = 0; j < pNumCols; j++)
            if(this->P[i][j] == 1) {
                auxBegin[i] = j;
                break;
            }//end-if
        
    }//end-for
    
    this->begin = auxBegin;
                
}//end-BuildProbabilityMatrix()

/* Method for computing the binary expansion of a given probability in [0, 1] */
void Samplers::BinaryExpansion(Vec< Vec<int> >& auxP, RR probability, int precision, int index) {
    
    RR pow;
    int i, j;
    i = -1;
    j = 0;
    
    while(probability > 0 && j < precision) {
        pow = power2_RR(i--); // 2^{i}
        if(pow <= probability) {
            auxP[j][index] = 1;
            probability -= pow;
        } else
            auxP[j][index] = 0;
        j++;
    }//end-while
            
}//end-BinaryExpansion()

void Samplers::PrintMatrix(const string& name, const Vec< Vec<int> >& matrix) {
    
    cout << "\n/** " << name << " **/" << endl;
    for(int i = 0; i < matrix.length(); i++) {
        for(int j = 0; j < matrix[0].length(); j++)
            cout << matrix[i][j] << " ";
        cout << endl;
    }//end-for
    
}//end-PrintVectorZZX()