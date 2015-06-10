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

Samplers::Samplers(int k) {
    
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
    
    int phi;
    phi = this->EulerPhiPowerOfTwo(k); //It computes \phi(m), with m = 2^k
    this->phi.SetLength(phi+1);
    SetCoeff(this->phi, 0, 1);
    SetCoeff(this->phi, phi, 1);
    
    //It's used for sampling from lattice
    this->BuildVandermondeMatrix(k);
    
}//end-Samplers()

Samplers::~Samplers() {};

/* Probability associated with x */
RR Samplers::Rho(RR sigma, RR x) {
    
    if(x == 0)
        return to_RR(1);
    
    return exp(-(power(x/sigma, 2))/2.0);
    
}//end-Rho()

RR Samplers::Probability(RR x, RR sigma) {
    RR S = sigma*sqrt(2*ComputePi_RR());
    RR overS = 1/S;
    
    if(x == to_RR(0))
        return overS;
    
    return overS*exp(-(power(x/sigma, 2))/2.0);
    
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
        
    ZZ powerOmega, yPrime, yBar;
    int bit, curve, i, invalidSample, zero, s, S, x;
    unsigned less, lessEqual, equal, greater;
    long b;
        
    powerOmega = power2_ZZ(omega); // 2^{\omega}
    invalidSample = 14*to_int(sigma);
    bit = 0;
    
    /* Estimated number of iterations required to obtain all requested samples 
     * in one iteration of PolyGeneratorZiggurat algorithm */
    for(int index = 0; index < 6; index++) { 
        i = RandomBnd(m) + 1; // Sample a random value in {0, ..., m-1} and add one to the result, e.g. select a rectangle in {1, ..., m}
        s = 1 - 2*RandomBits_long(1); // Sample a random signal s
        x = RandomBnd(this->X_ZZ[i] + 1); // Sample a x value between 0 and floor(x_i)
        b = RandomBits_long(1);
        yPrime = RandomBnd(powerOmega - 1);
        yBar = yPrime * (this->Y_ZZ[i-1] - this->Y_ZZ[i]);
        curve = to_int(to_RR(powerOmega) * (Rho(sigma, to_RR(x)) - this->Y[i]));
        
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
                        
        less = lessThan(x, curve);
        equal = isEqual(x, curve);
        
        /* Third case: the sampled x is below to the curve */
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
    
    // The approximation error must be less than 2^{-n}
    RR statDistance = power2_RR(-((long)n));
    
    RR tailcut, y0;
    int c;
    bool validPartition = 0;
    
    this->X.SetLength(m + 1);
    this->Y.SetLength(m + 1);
            
    tailcut = to_RR(tail)*sigma;
    
    while(tailcut < to_RR(14)*sigma) {
        
        c = 1;        
        y0 = to_RR(0);
        this->X[m] = tailcut;        
        
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
RR Samplers::DZRecursion(int m, int c, RR sigma) {
    
    RR inv, minus2, overM, over2, S;
    
    minus2 = to_RR(-2);
    overM = to_RR(1)/to_RR(m);
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
    this->Y[m] = Rho(sigma, this->X[m]);
    
    inv = minus2 * log(S/to_RR(TruncToZZ(this->X[m]) + to_ZZ(1)));
    
    this->X[m-1] = sigma * sqrt(inv);
    this->Y[m-1] = Rho(sigma, this->X[m-1]);
    
    for(int i = m-2; i > 0; i--) {
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
    
    this->BuildProbabilityMatrix(precision, tailcut, sigma);
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
void Samplers::BuildProbabilityMatrix(int precision, int tailcut, RR sigma) {
    
    // The random variable consists of elements in [-tailcut*sigma, tailcut*sigma]
    int i, bound;
    RR probOfX;
    
    bound = tailcut*to_int(sigma);    
    
    this->P.SetLength(precision);    
    for(i = 0; i < this->P.length(); i++)
        this->P[i].SetLength(2*bound+1);

    i = 2*bound;    
    for(int x = -bound; x <= bound, i>= 0; x++, i--) {
        probOfX = Probability(to_RR(x), sigma);
        BinaryExpansion(probOfX, precision, i);
    }//end-for
       
    // Uncomment this line if you want to preview the probability matrix P
//    this->PrintMatrix("Probability matrix", this->P);
    
}//end-BuildProbabilityMatrix()

/* Method for computing the binary expansion of a given probability in [0, 1] */
void Samplers::BinaryExpansion(RR probability, int precision, int index) {
    
    RR pow;
    int i, j;
    i = -1;
    j = 0;
    
    while(probability > 0 && j < precision) {
        pow = power2_RR(i--); // 2^{i}
        if(pow <= probability) {
            this->P[j][index] = 1;
            probability -= pow;
        } else
            this->P[j][index] = 0;
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
    b %= this->phi;
    return MulByXMod(b, this->phi);
}//end-Isometry()

/* Norm of a polynomial */
double Samplers::Norm(const ZZX& b) {
    
    Vec< complex<double> > mult;
    complex<double> sum;
    int aux, colsV, i, j, rowsV;
    
    colsV = this->V[0].length();
    rowsV = this->V.length();
    mult.SetLength(rowsV);
    
    /* Each row contains b(\psi_m^{i}) */
    for(i = 0; i < rowsV; i++) {
        sum = 0;
        for(j = 0; j < colsV; j++) {
            aux = to_int(b[j]);
            sum += this->V[i][j]*(complex<double>)aux;
        }//end-for
        mult[i] = sum;
    }//end-for
    
    sum = 0; 
    for(i = 0; i < rowsV; i++)
        sum += pow(mult[i].real(), 2.0) + pow(mult[i].imag(), 2.0);        
    
    return sqrt(sum).real();
    
}//end-Norm()

/* Computation of the inner product as matrix multiplication */
ZZ Samplers::InnerProduct(const ZZX& a, const ZZX& b) {
    // <a, b> = a^{T} * \bar{V^{T}} * V * b
    
    int colsV, i, j, rowsV;
    colsV = this->V[0].length();
    rowsV = this->V.length();
    
    Vec< Vec< complex<double> > > transpV;
    Vec< ZZX > C;
    transpV.SetLength(colsV);
    C.SetLength(colsV);
    
    for(i = 0; i < colsV; i++) {
        C[i].SetLength(colsV);
        transpV[i].SetLength(rowsV);
    }//end-for
    
    /* Transposition of Vandermonde matrix V */
    for(i = 0; i < colsV; i++)
        for(j = 0; j < rowsV; j++)
            transpV[i][j] = this->V[j][i];
    
    /* Computation of conjugate of V^{T} */
    this->ConjugateOfMatrix(transpV);
    
    /* Integer matrix multiplication \bar{V^{T}} * V */
    this->ComplexMatrixMult(C, transpV, this->V);
    
    /* D = a^{T} * \bar{V^{T}} * V */
    ZZX D;
    ZZ sum;
    D.SetLength(colsV);
    
    for(i = 0; i < colsV; i++) {
        sum = to_ZZ(0);
        for(j = 0; j < colsV; j++)
            sum = sum + a[j]*C[j][i];
        D[i] = sum;
    }//end-for
    
    /* Final computation. The output is an integer */
    sum = 0;
    for(i = 0; i < colsV; i++)
        sum += D[i]*b[i];
    
    return sum;
    
}//end-InnerProduct()

/* The Vandermonde matrix computes the i-th m-th roots of unity 
 * and their exponentiations */
void Samplers::BuildVandermondeMatrix(int k) {
    
    complex<double> rootOfUnity;
    double pi;
    int i, index, j, m, phi;
    
    pi = to_double(ComputePi_RR());
    m = pow(2, k);
    phi = this->EulerPhiPowerOfTwo(k);
    
    this->V.SetLength(phi);
    
    index = 0;
    for(i = 0; i < m; i++) {        
        this->V[index].SetLength(phi);
        if(GCD(i, m) == 1) {
            rootOfUnity = std::polar(1.0, (double)((2*pi*i)/(double)m));        
            for(j = 0; j < phi; j++)
                this->V[index][j] = pow(rootOfUnity, j);
            index++;
        }//end-if        
    }//end-for
    
}//end-BuildVandermondMatrix()

/* Computation of Euler's phi function for m = 2^k, p = 2 */
int Samplers::EulerPhiPowerOfTwo(int k) {
    return pow(2.0, k-1);    
}//end-EulerPhiPowerOfTwo()

/* It computer the conjugate of each element in the matrix */
void Samplers::ConjugateOfMatrix(Vec< Vec< complex<double> > >& M) {
    
    int cols, rows;
    cols = M[0].length();
    rows = M.length();
    
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            M[i][j] = std::conj(M[i][j]);
    
}//end-ConjugateOfMatrix()

/* Matrix multiplication of complex numbers generating an integer matrix */
void Samplers::ComplexMatrixMult(Vec< ZZX >& c, const Vec< Vec< complex<double> > >& a, const Vec< Vec< complex<double> > >& b) {       

    int colsA, colsB, i, j, k, rowsA, rowsB;
    colsA = a[0].length();
    colsB = b[0].length();
    rowsA = a.length();
    rowsB = b.length();
    
    if(colsA == rowsB) {
        
        complex<double> sum;
        
        c.SetLength(rowsA);
        for(i = 0; i < rowsA; i++)
            c[i].SetLength(colsB);        
        
        for (k = 0; k < rowsA; k++) {
          for (j = 0; j < colsB; j++) {
            sum = 0;
            for (i = 0; i < rowsB; i++)
              sum = sum + a[k][i]*b[i][j];
            c[k][j] = to_ZZ(sum.real());
          }//end-for
        }//end-for
        
    }//end-if
    
}//end-Mult()

void Samplers::PrintMatrix(const string& name, const Vec< Vec<int> >& matrix) {
    
    cout << "\n/** " << name << " **/" << endl;
    for(int i = 0; i < matrix.length(); i++) {
        for(int j = 0; j < matrix[0].length(); j++)
            cout << matrix[i][j] << " ";
        cout << endl;
    }//end-for
    
}//end-PrintVectorZZX()

/* Sampling a lattice point from a Gaussian centered in zero */
ZZX Samplers::GaussianSamplerFromLattice(const Vec<ZZX>& B, const Vec<ZZX>& BTilde, RR sigma, int precision, int tailcut) {
    
    ZZX c, mult;
    RR sigmaI;
    ZZ zero;
    double d, norm;
    int i, m, sample;
    
    m = B.length();
    c.SetLength(m);
    zero = to_ZZ(0);
   
    for(i = 0; i < m; i++) //Centered in zero
        c[i] = zero;
    
    for(i = m-1; i >= 0; i--) {
        norm = this->Norm(BTilde[i]);
        d = to_double(this->InnerProduct(c, BTilde[i]))/(pow(norm, 2.0));
        sigmaI = sigma/to_RR(norm);
        sample = this->KnuthYao(precision, tailcut, sigmaI);
        sample = (int)(sample + d);
        mul(mult, B[i], (long)(sample));        
        sub(c, c, mult);
    }//end-for
    
    for(i = 0; i < m; i++)
        c[i] = -c[i];
    
    return c;
    
}//end-GaussianSamplerFromLattice()

/* Gram-Schmidt reduced basis generation */
void Samplers::FasterIsometricGSO(Vec<ZZX>& BTilde, Vec<ZZ>& C, Vec<double>& D, const Vec<ZZ_pX>& B, int k) {
    
    ZZX isometry;
    ZZX B1; // B[0] reduced modulo Phi_m(x)
    ZZX V; // Updated at each iteration
    ZZX V1; // It stores the first value of V
    double div, Cdouble;
    int i, m, phi;
    
    m = B.length();
    phi = this->EulerPhiPowerOfTwo(k);
    
    B1.SetLength(phi);
    V.SetLength(phi);
    V1.SetLength(phi);
    C.SetLength(m);
    D.SetLength(m);
    
    BTilde.SetLength(m);
    for(i = 0; i < BTilde.length(); i++)
        BTilde[i].SetLength(phi);
    
    B1 = to_ZZX(B[0]) % this->phi;
    
    BTilde[0] = B1; V = B1; V1 = B1;
    
    C[0] = this->InnerProduct(V1, this->Isometry(B1));
    D[0] = pow(this->Norm(B1), 2.0);
    
    for(i = 0; i < m-1; i++) {
        
        Cdouble = to_double(C[i]);
        div = Cdouble/D[i];
        isometry = this->Isometry(BTilde[i]);
        
        BTilde[i+1] = isometry - this->Mult(V, div, phi); // V_k*C_k/D_k
        V = V - this->Mult(isometry, div, phi);
        C[i+1] = this->InnerProduct(V1, this->Isometry(BTilde[i+1]));
        D[i+1] = D[i] - div*Cdouble;
        
    }//end-for
    
}//end-FasterIsometricGSO()

ZZX Samplers::Mult(ZZX V, double c, int phi) {
    for(int i = 0; i < phi; i++)
        V[i] = to_ZZ(to_double(V[i]) * c);
    return V;
}//end-Mult()