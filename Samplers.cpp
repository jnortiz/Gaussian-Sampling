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
#include <NTL/quad_float.h>

using namespace NTL;
using namespace std;

Samplers::Samplers(int q, const ZZ_pX& f) {
    
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
            polynomial[i] = this->Ziggurat(m, n, sigma, omega);
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
Vec<int> Samplers::PolyGeneratorKnuthYao(int dimension, int precision, int tailcut, RR sigma, RR c) {
    
    
    cout << "\n[*] Knuth-Yao Gaussian sampling" << endl;
        
    this->BuildProbabilityMatrix(precision, tailcut, sigma, c);
    cout << "[*] Probability matrix building status: Pass!" << endl;
    
    Vec<int> polynomial;
    int bound, center, sample;
    
    polynomial.SetLength((long)dimension);
    bound = tailcut*to_int(sigma);
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
int Samplers::KnuthYao(int tailcut, RR sigma, RR c) {

    int bound, center, col, d, invalidSample, pNumRows, pNumCols, r, S, signal;
    unsigned enable, hit;
    
    bound = tailcut*to_int(sigma);
    center = to_int(c);
    d = 0; //Distance
    hit = 0;
    signal = 1 - 2*RandomBits_long(1); // Sample a random signal s    
    invalidSample = bound+1;
    pNumRows = this->P.length(); // Precision
    pNumCols = this->P[0].length();    
    
    S = 0;
    
    for(int row = 0; row < pNumRows; row++) {
        
        r = RandomBits_long(1); // Random choice between 0 and 1
        d = 2*d + r; // Distance calculus
        
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
void Samplers::BuildProbabilityMatrix(int precision, int tailcut, RR sigma, RR c) {
    
    RR::SetPrecision(to_long(precision));

    Vec< Vec<int> > auxP;
    Vec<int> auxBegin;
    
    // The random variable consists of elements in [c-tailcut*sigma, c+tailcut*sigma]
    int i, j, bound, pNumCols, pNumRows;
    RR probOfX;    
    
    bound = tailcut*to_int(sigma);
       
    auxP.SetLength(precision);
    for(i = 0; i < auxP.length(); i++)
        auxP[i].SetLength(bound+1);

    div(probOfX, Probability(to_RR(0) + c, sigma, c), to_RR(2));
    BinaryExpansion(auxP, probOfX, precision, bound);
    
    i = bound-1;
    for(int x = 1; x <= bound, i >= 0; x++, i--) {
        probOfX = Probability(to_RR(x) + c, sigma, c);
        BinaryExpansion(auxP, probOfX, precision, i);
    }//end-for
    
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

ZZX Samplers::GaussianSamplerFromLattice(const Vec<ZZX>& B, const mat_RR& BTilde, RR sigma, int precision, int tailcut, ZZX center, int n) {

    // Precision of floating point operations
    RR::SetPrecision(precision);    

    ZZX C, mult, sample;
    int Z;
    
    vec_RR auxC;
    RR d, norm, sigma_i;
    int i, j, mn;
    
    mn = B.length(); // mn = (m1 + m2) * n
    
    auxC.SetLength(n);
    C.SetLength(n);
    mult.SetLength(n);
    sample.SetLength(n);
                
    C = center; // Center of the lattice        
    
    /* The inner product and norm operations are taken 
     * in the inner product space H */
    for(i = mn-1; i >= 0; i--) {
        
        norm = this->Norm(BTilde[i], n);
        
        for(j = 0; j < n; j++)
            auxC[j] = to_RR(C[j]);
        
        // The new center for the discrete Gaussian
        div(d, this->InnerProduct(auxC, BTilde[i], n), this->InnerProduct(BTilde[i], BTilde[i], n));        
        
        // And the new standard deviation
        div(sigma_i, sigma, norm);
        
        this->BuildProbabilityMatrix(precision, tailcut, sigma_i, d);             
        
        Z = this->KnuthYao(tailcut, sigma_i, d);

        mul(mult, B[i], Z);       
        
        sub(C, C, mult);
        
    }//end-for
    
    sub(sample, center, C);
    
    return sample;
    
}//end-GaussianSamplerFromLattice()

RR Samplers::BlockGSO(mat_RR& BTilde, const Vec<ZZX>& B, int m, int n, int precision) {
    
    // Precision of floating point operations
    RR::SetPrecision(precision);    
    
    cout << "\n[!] Norm of basis B: " << this->NormOfBasis(B, m*n, n);    
    cout << "\n[*] Block_GSO status: ";
    
    BTilde.SetDims(m*n, n);
    
    mat_RR outRot, outGSO; // n x n matrices
    vec_RR conv, ortho, mult;
    RR mu, innerpr, norm;
    int i, j, k;
    
    conv.SetLength(n);
    ortho.SetLength(n);
    mult.SetLength(n);
    
    for(i = 0; i < m; i++) { // For each isometric block i in [1, ..., m]

        for(k = 0; k < n; k++) {
            conv[k] = to_RR(B[i*n][k]); // Copy the vector whose expansion is an isometric basis
            ortho[k] = to_RR(0);
        }//end-for

         /* For each expansion of a vector makes B[i*n] orthogonal 
          * to the previous vectors */
        for(j = 0; j < i*n; j++) {
            innerpr = this->InnerProduct(conv, BTilde[j], n);            
            norm = this->InnerProduct(BTilde[j], BTilde[j], n);            
            div(mu, innerpr, norm);           
            mul(mult, BTilde[j], mu);
            add(ortho, ortho, mult);
        }//end-for
        
        sub(BTilde[i*n], conv, ortho);            
                
        // Expansion of the vector that is already orthogonal to its predecessors
        // In (Lyubashevsky, and Prest, 2015), it uses B[i*n] instead
        this->rot(outRot, BTilde[i*n], n);
                
        // Computes its orthogonal basis
        this->FasterIsometricGSO(outGSO, outRot);
        
        // Copying the orthogonal basis of B[i*n] to the output
        for(j = 0; j < n; j++)
            BTilde[i*n + j] = outGSO[j];
        
    }//end-for    
    
    cout << "Pass!" << endl;
    
    norm = this->NormOfBasis(BTilde, m*n, n);
    
    return norm;
    
}//end-BlockGSO()

void Samplers::FasterIsometricGSO(mat_RR& BTilde, const mat_RR& B) {
    
    /* This implementation follows the Algorithm 3 
     * of (Lyubashevsky, and Prest, 2015) */
    
    vec_RR C, D, isometry, mult, V;
    RR CD;
    int n;
    
    n = B.NumCols(); // B is a square matrix
    
    BTilde.SetDims(n, n);
    C.SetLength(n);
    D.SetLength(n);
        
    isometry.SetMaxLength(n);
    mult.SetMaxLength(n);
    V.SetMaxLength(n);
    
    BTilde[0] = B[0];
    V = B[0];    
    C[0] = this->InnerProduct(V, this->Isometry(BTilde[0], n), n);
    D[0] = this->InnerProduct(BTilde[0], BTilde[0], n);
    
    for(int i = 0; i < n-1; i++) {        
        div(CD, C[i], D[i]);        
        isometry = this->Isometry(BTilde[i], n);         
        mul(mult, V, CD);                
        sub(BTilde[i+1], isometry, mult);            
        mul(mult, isometry, CD);        
        sub(V, V, mult);                        
        C[i+1] = this->InnerProduct(B[0], this->Isometry(BTilde[i+1], n), n);         
        sub(D[i+1], D[i], CD*C[i]);                
    }//end-for
    
}//end-FasterIsometricGSO()

void Samplers::Rot(Vec<ZZX>& A, const Vec<ZZ_pX>& a, int m, int n) {
    
    Vec<ZZX> out;                            
    int i, index, j;
    
    A.SetLength(m*n);        
    
    for(i = 0; i < m*n; i++)
        A[i].SetMaxLength(n);
    
    for(i = 0, index = 0; i < m; i++) {
        this->rot(out, a[i], n);
        for(j = 0; j < n; j++)
            A[index++] = out[j];                
    }//end-for
    
}//end-Rot()

/* Rot_f(b) operation: generates a circulant matrix which contains the (n-1) rotations of vector b */
void Samplers::rot(Vec<ZZX>& out, const ZZ_pX& b, int n) {
    
    ZZX auxB, isometry;
    
    out.SetLength(n);
    auxB.SetLength(n);
    
    auxB = to_ZZX(b);        
    out[0] = auxB;    
    isometry = auxB;
    
    for(int i = 1; i < n; i++) {
        isometry = this->Isometry(isometry, n);
        out[i] = isometry;
    }//end-for
   
}//end-Rot()

void Samplers::rot(mat_RR& out, const vec_RR& b, int n) {
        
    vec_RR isometry;
    
    out.SetDims(n, n);    
    
    out[0] = b;    
    isometry = b;    
    
    for(int i = 1; i < n; i++) {
        isometry = this->Isometry(isometry, n);
        out[i] = isometry;
    }//end-for
    
}//end-Rot()

vec_RR Samplers::Isometry(vec_RR& b, int n) {
    
    if(IsZero(b))
        return b;
    
    vec_RR out;
    
    out.SetLength(n);
        
    out[0] = -b[n-1];
    
    for(int i = 1; i < n; i++)
        out[i] = b[i-1];
    
    return out;
    
}//end-Isometry()

/* Giving a polynomial g, out contains (b*x)%f(x) */
ZZX Samplers::Isometry(ZZX& b, int n) {
    
    if(IsZero(b))
        return b;
    
    ZZX out;
    
    out.SetLength(n);
        
    out[0] = -b[n-1];    
    
    for(int i = 1; i < n; i++)
        out[i] = b[i-1];
    
    return out;
    
}//end-Isometry()

RR Samplers::NormOfBasis(const mat_RR& B, int m, int n) {
    
    RR norm, normB;    
    
    normB = to_RR(0);
    
    for(int i = 0; i < m; i++) {
        norm = this->Norm(B[i], n);
        if(norm > normB)
            normB = norm;
    }//end-for
    
    return normB;
    
}//end-NormOfBasis()

RR Samplers::NormOfBasis(const Vec<ZZX>& B, int m, int n) {
    
    RR norm, normB;    
    
    normB = to_RR(0);
    
    for(int i = 0; i < m; i++) {
        norm = this->Norm(B[i], n);
        if(norm > normB)
            normB = norm;
    }//end-for
    
    return normB;
    
}//end-NormOfBasis()

/* Norm of a polynomial */
RR Samplers::Norm(const ZZX& b, int n) {
    
    ZZ norm = to_ZZ(0);

    for(int i = 0; i < n; i++)
        norm += b[i]*b[i];
    
    return SqrRoot(to_RR(norm));
        
}//end-Norm()

RR Samplers::Norm(const vec_RR& b, int n) {
    
    RR norm, mult;
    norm = to_RR(0);

    for(int i = 0; i < n; i++) {
        mul(mult, b[i], b[i]);
        norm += mult;
    }//end-for
    
    return SqrRoot(norm);
        
}//end-Norm()

RR Samplers::InnerProduct(const vec_RR& a, const vec_RR& b, int n) {
    
    RR innerp, mult;
    innerp = to_RR(0);

    for(int i = 0; i < n; i++) {
        mul(mult, a[i], b[i]);
        innerp += mult;
    }//end-for
    
    return innerp;
    
}//end-InnerProduct()

void Samplers::PrintMatrix(const string& name, const Vec< Vec<int> >& matrix) {
    
    cout << "\n/** " << name << " **/" << endl;
    for(int i = 0; i < matrix.length(); i++) {
        for(int j = 0; j < matrix[0].length(); j++)
            cout << matrix[i][j] << " ";
        cout << endl;
    }//end-for
    
}//end-PrintVectorZZX()