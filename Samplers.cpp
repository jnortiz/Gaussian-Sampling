/* 
 * File:   Samplers.cpp
 * Author: jnortiz
 * 
 * Created on April 24, 2015, 3:51 PM
 */

#include "Samplers.h"
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/RR.h>
#include <NTL/matrix.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <NTL/c_lip.h>
#include <complex>
#include <NTL/ZZX.h>

using namespace NTL;
using namespace std;

Samplers::Samplers(int k, int q, const ZZ_pX& f) {
    
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
    
    ZZ_p::init(conv<ZZ>(q)); // Coefficients modulo
    this->f = f;
    ZZ_pE::init(f); // Ring elements modulo
    
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
    
    cout << "[*] Ziggurat Gaussian sampling" << endl;

    /* Output precision setup */
    RR::SetOutputPrecision((long)n);
    
    Vec<int> polynomial;    
    polynomial.SetLength((long)dimension);

    // Only for statistical purposes
    int bound, nIterations, nPositive, nNegative, nZero, samplesGen;
    bound = to_int(tail*sigma);
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
            if(polynomial[i] <= bound)
                samplesGen++;
        }//end-for
        nIterations++;
    }while(samplesGen < dimension);

    cout << "[*] All samples were successfully generated in " << nIterations << " iteration(s)." << endl;
    
    return polynomial;
    
}//end-PolyGeneratorZiggurat()

/* Algorithm for sampling from discrete distribution using rejection method.
 * It does not avoid Rho function computation */
int Samplers::Ziggurat(int m, RR sigma, int omega) {
        
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
    
}//end-ZigguratO()

/* DZCreatePartition defines the x and y axes of rectangles in the Gaussian distribution */
void Samplers::DZCreatePartition(int m, RR sigma, int n, int tail) {
    /*
     * Parameters description:
     * m: number of rectangles
     * sigma: Gaussian distribution parameter
     * n: bit precision
     */
    
    /* Output precision setup */
    RR::SetOutputPrecision((long)n);
    
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
Vec<int> Samplers::PolyGeneratorKnuthYao(int dimension, int precision, int tailcut, RR sigma, RR c) {
    
    cout << "\n[*] Knuth-Yao Gaussian sampling" << endl;
    
    /* Output precision setup */
    RR::SetOutputPrecision(to_long(precision));
    
    this->BuildProbabilityMatrix(precision, tailcut, sigma, c);
    cout << "[*] Probability matrix building status: Pass!" << endl;
    
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
    
    int bound, col, d, invalidSample, pNumRows, pNumCols, r, searchRange, S;
    unsigned enable, hit;
    
    bound = tailcut*to_int(sigma);
    d = 0; //Distance
    hit = 0;
    invalidSample = 3*bound;
    pNumRows = precision;
    pNumCols = 2*bound+1;    
    
    /* Approximated search range required to obtain all samples with only one iteration 
     * in PolyGeneratorKnuthYao() algorithm */
    searchRange = pNumRows/4;
    S = 0;
    
    for(int row = 0; row < searchRange; row++) {
        
        r = RandomBits_long(1); // Random choice between 0 and 1
        d = 2*d + r; // Distance calculus
        
        for(col = 0; col < pNumCols; col++) {
            
            d = d - this->P[row][col];
            
            enable = (unsigned)(d + 1); // "enable" turns 0 iff d = -1
            enable = 1 ^ ((enable | -enable) >> 31) & 1; // "enable" turns 1 iff "enable" was 0
             
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
 * [-tailcut*\floor(sigma), +tailcut*\floor(sigma)] */
void Samplers::BuildProbabilityMatrix(int precision, int tailcut, RR sigma, RR c) {
        
    Vec< Vec<int> > auxP;
    // The random variable consists of elements in [-tailcut*sigma, tailcut*sigma]
    int bound, center, i, upperBound;
    RR probOfX;
    
    bound = tailcut*to_int(sigma);
    center = to_int(c);
    
    auxP.SetLength(precision);
    for(i = 0; i < auxP.length(); i++)
        auxP[i].SetLength(2*bound+1);

    upperBound = center + bound;
    i = 2*bound;    
    for(int x = (center - bound); x <= upperBound, i >= 0; x++, i--) {
        probOfX = Probability(to_RR(x), sigma, c);
        BinaryExpansion(auxP, probOfX, precision, i);
    }//end-for
       
    this->P = auxP;
    // Uncomment this line if you want to preview the probability matrix P
//    this->PrintMatrix("Probability matrix", this->P);
    
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

/* Generation of a length-degree polynomial modulo q */
void Samplers::PolyGenerator(ZZX& b, int length, int q) {
    b.SetLength(length);
    for(int i = 0; i < length; i++)
        b[i] = NTL::RandomBnd(q);
}//end-PolyGenerator()

/* Giving a polynomial g, out contains (b*x)%phi(x) */
ZZX Samplers::Isometry(ZZX& b) {
    
    if(IsZero(b))
        return b;
    
    return MulByXMod(b, to_ZZX(this->f));
    
}//end-Isometry()

/* Norm of a polynomial */
RR Samplers::Norm(const ZZX& b, int n) {
    
    ZZ norm = to_ZZ(0);

    for(int i = 0; i < n; i++)
        norm += b[i]*b[i];
    
    return SqrRoot(to_RR(norm));
        
}//end-Norm()

/* Computation of the inner product as matrix multiplication */
ZZ Samplers::InnerProduct(const ZZX& a, const ZZX& b, int n) {
    
    ZZ innerp = to_ZZ(0);

    for(int i = 0; i < n; i++)
        innerp += a[i]*b[i];
    
    return innerp;    
    
}//end-InnerProduct()

/* Computation of Euler's phi function for m = 2^k, p = 2 */
int Samplers::EulerPhiPowerOfTwo(int k) {
    return pow(2.0, k-1);    
}//end-EulerPhiPowerOfTwo()

void Samplers::PrintMatrix(const string& name, const Vec< Vec<int> >& matrix) {
    
    cout << "\n/** " << name << " **/" << endl;
    for(int i = 0; i < matrix.length(); i++) {
        for(int j = 0; j < matrix[0].length(); j++)
            cout << matrix[i][j] << " ";
        cout << endl;
    }//end-for
    
}//end-PrintVectorZZX()

RR Samplers::GramSchmidtProcess(Vec<ZZX>& BTilde, const Vec<ZZ_pX>& B, int n) {
    
    ZZX mult;
    ZZ innerpr;
    RR mu, norm;
    int i, j, m;
    m = B.length();
    
    // Gram-Schmidt reduced basis
    BTilde.SetLength(m);
    mult.SetLength(n);

    for(i = 0; i < m; i++)
        BTilde[i].SetLength(n); // n coefficients        
    
    for(i = 0; i < m; i++) {
        
        BTilde[i] = to_ZZX(B[i]);
        
        for(j = 0; j < i; j++) {            
            
            innerpr = this->InnerProduct(to_ZZX(B[i]), BTilde[j], n);
            norm = to_RR(this->InnerProduct(BTilde[j], BTilde[j], n));
            
            if(norm == to_RR(0)) {
                cout << "Division by zero." << endl;
                return to_RR(-1);
            }//end-if
                
            div(mu, to_RR(innerpr), norm);            
            this->Mult(mult, BTilde[j], mu, n);
            sub(BTilde[i], BTilde[i], mult);            
            
        }//end-for
        
    }//end-for    
    
    RR norm1, norm2, normB, normBTilde;    
    normB = to_RR(0);
    normBTilde = to_RR(0);
    cout << endl;
    for(i = 0; i < m; i++) {
        norm1 = this->Norm(to_ZZX(B[i]), n);
        norm2 = this->Norm(BTilde[i], n);
        if(norm1 > normB)
            normB = norm1;
        if(norm2 > normBTilde)
            normBTilde = norm2;
    }//end-for
    cout << "Norm of B: " << normB << endl;
    cout << "Norm of GS reduced basis: " << normBTilde << endl;
    
    return normBTilde;
    
}//end-GramSchmidtProcess() 

// TODO:
void Samplers::SampleD(const Vec<ZZ_pX>& B, const Vec<ZZX>& BTilde, RR sigma, RR c, RR norm) {
    
    int m;
    m = B.length();
    
    if(sigma <= norm * log(m)/log(2)) {
        cout << "\n[!] Sigma must be greater or equal to the norm of GSO B times log(m)." << endl;
        return;
    }//end-if        
    
};

///* Sampling a lattice point from a Gaussian centered in zero */
//ZZX Samplers::GaussianSamplerFromLattice(const Vec<ZZ_pX>& B, const Vec<ZZX>& BTilde, RR sigma, int precision, int tailcut, ZZX c, int k) {
//        
//    RR::SetOutputPrecision(to_long(precision));
//
//    Vec<ZZX> auxB, C;
//    Vec<int> Z;
//    ZZX C0, mult, sample;
//    ZZ norm;
//    RR d, sigma_i;
//    int i, m, n, phi;
//    
//    m = B.length(); // m = m1 + m2 ~ 55
//    n = pow(2, k); // n = 2^k
//    phi = this->EulerPhiPowerOfTwo(k);
//    
//    Z.SetLength(m);
//    C0.SetLength(n);
//    mult.SetLength(n);
//    sample.SetLength(n);
//            
//    auxB.SetLength(m);
//    C.SetLength(m);
//    for(i = 0; i < m; i++) {
//        auxB[i].SetLength(n);
//        C[i].SetLength(n);
//        auxB[i] = to_ZZX(B[i]); // Basis conversion from ZZ_pX to ZZX (no changes are made in coefficients)
//    }//end-for
//    
//    C[m-1] = c; // Center of the lattice        
//    
//    /* The inner product and norm operations are taken 
//     * in the inner product space H */
//    for(i = m-1; i > 0; i--) {
//        
//        norm = this->Norm(BTilde[i], n);
//        
//        if(norm == to_ZZ(0)) {
//            cout << "\n[!] Error: division by zero ahead. Norm is zero. Aborting..." << endl;
//            return to_ZZX(-1);
//        }//end-if
//                    
//        // The new center for the discrete Gaussian
//        div(d, to_RR(this->InnerProduct(C[i], BTilde[i], n)), to_RR(this->InnerProduct(BTilde[i], BTilde[i], n)));        
//
//        // And the new standard deviation
//        div(sigma_i, sigma, to_RR(norm));
//        
//        // Problem: what if sigma_i < 1? Should be possible to sample from a discrete Gaussian distribution?
//        // It can be a problem in GSO algorithm...
//        if(sigma_i < 1) // If sigma is smaller than 1, there is no integer to be sampled
//            sigma_i += 1;
//
//        this->BuildProbabilityMatrix(precision, tailcut, sigma_i, d); // For Knuth-Yao algorithm                   
//        Z[i] = this->KnuthYao(precision, tailcut, sigma_i); // Sampling from a different Gaussian each iteration
//
//        mul(mult, auxB[i], (long)(Z[i])); // Multiply Z[i] with all coefficients of auxB        
//        sub(C[i-1], C[i], mult);
//        
//    }//end-for
//    
//    norm = this->Norm(BTilde[i], n);
//
//    if(norm == to_ZZ(0)) {
//        cout << "\n[!] Error: division by zero ahead. Norm is zero. Aborting..." << endl;
//        return to_ZZX(-1);
//    }//end-if
//
//    div(d, to_RR(this->InnerProduct(C[i], BTilde[i], n)), to_RR(this->InnerProduct(BTilde[i], BTilde[i], n)));        
//    div(sigma_i, sigma, to_RR(norm));
//    
//    if(sigma_i < 1)
//        sigma_i += 1;
//
//    this->BuildProbabilityMatrix(precision, tailcut, sigma_i, d);                   
//    Z[i] = this->KnuthYao(precision, tailcut, sigma_i);
//
//    mul(mult, auxB[i], (long)(Z[i]));  
//    sub(C0, C[i], mult);    
//    sub(sample, c, C0);
//    
//    return sample;
//    
//}//end-GaussianSamplerFromLattice()

///* Gram-Schmidt reduced basis generation */
//void Samplers::FasterIsometricGSO(Vec<ZZX>& BTilde, Vec<ZZ>& C, Vec<ZZ>& D, const Vec<ZZ_pX>& B, int n) {
//    
//    cout << "\n[*] Running FasterIsometricGSO: ";
//    
//    /* This implementation follows the Algorithm 3 
//     * of (Lyubashevsky, and Prest, 2015) */
//    
//    ZZX isometry, mult;
//    ZZX B1; // B[0] reduced modulo Phi_m(x)
//    ZZX V; // Updated at each iteration
//    ZZ CD, norm;
//    int i, m;
//    
//    // Lengths:
//    m = B.length();
//    
//    // Vectors:
//    BTilde.SetLength(m);
//    C.SetLength(m);
//    D.SetLength(m);
//    
//    for(i = 0; i < m; i++)
//        BTilde[i].SetMaxLength(n);
//    
//    // Polynomials:
//    isometry.SetMaxLength(n);
//    mult.SetMaxLength(n);
//    B1.SetMaxLength(n);
//    V.SetMaxLength(n);
//    
//    B1 = to_ZZX(B[0]);
//    BTilde[0] = B1;
//    V = B1;
//    
//    C[0] = this->InnerProduct(B1, to_ZZX(this->Isometry(B1)), n);
//    norm = this->Norm(B1, n);
//    mul(D[0], norm, norm);
//    
//    for(i = 0; i < m-1; i++) {
//        
//        if(D[i] == to_ZZ(0)) {
//            cout << "\n[!] Error: division by zero ahead. Norm is zero. Aborting..." << endl;
//            return;
//        }//end-if
//        
//        div(CD, C[i], D[i]);
//        isometry = this->Isometry(BTilde[i]);
//        mul(mult, V, CD);
//        sub(BTilde[i+1], isometry, mult); 
//        mul(mult, isometry, CD);
//        sub(V, V, mult);        
//        C[i+1] = this->InnerProduct(B1, this->Isometry(BTilde[i+1]), n);        
//        mul(D[i+1], CD, C[i]);
//        sub(D[i+1], D[i], D[i+1]);
//        
//    }//end-for
//    
//    cout << "Pass!" << endl;
//        
//}//end-FasterIsometricGSO()

void Samplers::Mult(ZZX& out, const ZZX& V, RR c, int n) {    
    
    RR v;    
    
    if(!IsZero(V)) {
        for(int i = 0; i < n; i++) {
            v = to_RR(V[i]);
            out[i] = to_ZZ(v*c);
        }//end-for
    } else
        out = to_ZZX(0);
    
}//end-Mult()

RR Samplers::CoverageAreaZiggurat(RR sigma) {
    
    RR area, one;
    int m;
    
    one = to_RR(1);    
    m = this->X.length(); // Actually, m is the numbers of rectangles plus one
    
    area = to_RR(0);
    for(int i = 1; i < m-1; i++) // Riemann sum
        area += this->Probability(this->X[i], sigma, to_RR(0)) * (this->X[i] - this->X[i-1]);
    
    return area;
    
}//end-CoverageAreaZiggurat()