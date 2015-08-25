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

/* Knuth-Yao algorithm to obtain a sample from the discrete Gaussian */
int Samplers::KnuthYao(int tailcut, RR sigma, RR c) {

    int bound, center, col, d, invalidSample, pNumRows, pNumCols, S, signal;
    unsigned enable, hit;
    unsigned long r;
    
    bound = tailcut*to_int(sigma);
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
void Samplers::BuildProbabilityMatrix(int precision, int tailcut, RR sigma, RR c) {
    
    RR::SetPrecision(to_long(precision));

    Vec< Vec<int> > auxP;
    Vec<int> auxBegin;
    
    // The random variable consists of elements in [c-tailcut*sigma, c+tailcut*sigma]
    int i, j, bound, pNumCols, pNumRows, x;
    vec_RR probOfX;
    RR pow;
    
    bound = tailcut*to_int(sigma);
    
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

RR Samplers::Ziggurat(int m, RR sigma, int precision, RR tail) {
    
    cout << "[*] Ziggurat status: ";
    RR::SetPrecision((long)precision);
    RR::SetOutputPrecision((long)precision);
    
    /* Create partitioning first */
    RR v;    
    this->DZCreatePartition(m, sigma, precision, tail, v);
    
    /* Sampling phase */    
    RR c, U, z;
    int i;
    double ulong_size = (double)(sizeof(unsigned long)*8);
    
    c = to_RR(0);
    
    for(;;) {
    
        i = RandomBnd(m); // Sample a random value in {0, ..., m-1}
        U = to_RR((RandomBnd(ulong_size-1) + 1)/ulong_size);
        
        if(i > 0)
            z = U*this->X[i];
        else
            z = U*v/this->Probability(X[i], sigma, c);
        
        if(z > this->X[i+1])
            return z;
        
        if(i == 0)
            return this->NewMarsagliaTailMethod(this->X[m-1]);
        
        U = to_RR((RandomBnd(ulong_size-1) + 1)/ulong_size);
        
        if(i > 0 && U*(this->Probability(X[i], sigma, c) - this->Probability(this->X[i+1], sigma, c)) < (this->Probability(z, sigma, c) - this->Probability(this->X[i+1], sigma, c)))
            return z;
            
    }//end-for
    
    cout << "Pass!\n";
    
}//end-Ziggurat()

/* DZCreatePartition defines the x and y axes of rectangles in the Gaussian distribution */
void Samplers::DZCreatePartition(int m, RR sigma, int n, RR tail, RR& v) {
    
    cout << "\n[*] DZCreatePartition status: ";
    /* The Ziggurat algorithm was designed for centered Gaussian distributions; 
     i.e., c = 0. */
    
    /*
     * Parameters description:
     * m: number of rectangles
     * sigma: Gaussian distribution parameter
     * n: bit precision
     */
    
    Vec<RR> bestX, X;
    RR bestdiff, r, z;
    
    RR statDistance = power2_RR(-((long)n));
    RR overM = to_RR(1)/to_RR(m);
    RR minusOne = to_RR(-1);
    RR zero = to_RR(0);
    
    int i;
    int first = 1;
    
    X.SetLength(m);
    bestX.SetLength(m);            
    
    bestX[m-1] = minusOne;    
    z = minusOne;    
    r = (tail*sigma)/2; //r = x_m, the mth x-value
    bestdiff = r; // Arbitrary initial value
    
    while(z != zero && r > zero) { // If r = 0 then the area v is also 0 (what doesn't make sense...)

        z = this->DZRecursion(X, m, r, sigma, v);        

        if(z == minusOne && first) { // Error in "inv" or square root computation
            
            first = 0;
            add(r, r, 2*overM);
            div(overM, overM, to_RR(m));
            
            while(z != zero && r > zero) { // If r = 0 then the area v is also 0 (what doesn't make sense...)
                
                z = this->DZRecursion(X, m, r, sigma, v);

                if(z == minusOne)            
                    break;

                if(abs(z) < abs(bestdiff)) { // If the actual z is closest to zero, then that's the better partitioning until now
                    for(int i = 0; i < m; i++)
                        bestX[i] = X[i];                
                    bestdiff = z;
                }//end-if

                sub(r, r, overM);
                
            }//end-while            
            
        }//end-if
        
        if(z == minusOne && !first)
            break;
        
        if(abs(z) < abs(bestdiff)) { // If the actual z is closest to zero, then that's the better partitioning until now
            for(int i = 0; i < m; i++)
                bestX[i] = X[i];                
            bestdiff = z;
        }//end-if
        
        sub(r, r, overM);
        
    }//end-while
    
    if(z == zero)
        for(int i = 0; i < m; i++)
            bestX[i] = X[i];                    
               
    if(bestX[m-1] != -1) { // Some partitioning was found
        
        cout << "Pass!" << endl;        
        cout << "\nFinal z value: " << bestdiff << endl;
        cout << "Final r value: " << X[m-1] << endl;
        cout << "Statistical difference: " << to_RR(statDistance - bestdiff) << endl;
        
        cout << "\n[!] Partitioning for " << m << " rectangles: \n";
        this->X.SetLength(m);
        for(i = 0; i < this->X.length(); i++) {
            this->X[i] = bestX[i];
            cout << this->X[i] << endl;
        }//end-for
        
    } else // No valid partition was found
        cout << "Error!" << endl;
            
}//end-DZCreatePartition()

/* It is used in DZCreatePartition to define the distance y0 */
RR Samplers::DZRecursion(Vec<RR>& X, int m, RR r, RR sigma, RR& v) {
    
    RR zero = to_RR(0);
    RR c = zero; // Center of distribution
    RR interm, overPi;
    RR minusOne = to_RR(-1);
    div(overPi, minusOne, ComputePi_RR());
    
    X[m-1] = r;
    
    v = r*this->Probability(r, sigma, c) /*+ integrate r to infinity Probability(x) */;
    
    for(int i = (m-2); i >= 0; i--) {
        
        if(X[i+1] == zero)
            return minusOne;
        
        // TODO: general formula for rho^{-1}(x)
        // This inversion of rho(x) works only when variance is equal to 1/2*pi
        interm = overPi * log(v/X[i+1] + this->Probability(X[i+1], sigma, c));
        
        if(interm < zero)
            return minusOne;
            
        X[i] = sqrt(interm);
       
    }//end-for
    
    return (v - X[0] + X[0]*this->Probability(X[0], sigma, c));
        
}//end-DZCreatePartition()

/* Given a threshold r, this function returns a value in the tail region (Thomas et al., 2007) */
RR Samplers::NewMarsagliaTailMethod(RR r) {
    
    /* r is the right most x-coordinate in Ziggurat partitioning */
    RR a, b;
    RR x, y;
    double ulong_size = (double)(sizeof(unsigned long)*8);
    
    do {                        
        // Two variables with values in (0, 1)
        a = to_RR((RandomBnd(ulong_size-1) + 1)/ulong_size);
        b = to_RR((RandomBnd(ulong_size-1) + 1)/ulong_size);
        div(x, -log(a), r);
        y = -log(b);
    } while(y + y > x*x);
    
    return (r+x);
    
}//end-NewMarsagliaTailMethod()


void Samplers::GramSchmidtProcess(Vec< Vec<double> >& T_ATilde, const Vec< Vec<int> >& T_A, int n) {
    
    cout << "\n[*] Gram-Schmidt process status: ";
        
    int i, j, k, m;
    double mu;
    
    m = T_A.length();
    
    T_ATilde.SetLength(m);
    for(i = 0; i < m; i++)
        T_ATilde[i].SetLength(m*n);
                
    for(i = 0; i < m; i++) {
        this->CopyIntToDoubleVec(T_ATilde[i], T_A[i]);
        for(j = 0; j < i; j++) {
            mu = this->InnerProduct(T_A[i], T_ATilde[j])/this->InnerProduct(T_ATilde[j], T_ATilde[j]);
            for(k = 0; k < (m*n); k++)
                T_ATilde[i][k] = T_ATilde[i][k] - mu*T_ATilde[j][k];            
        }//end-for
    }//end-for
    
    cout << "Pass!";
    
}//end-GramSchmidtProcess() 

void Samplers::CopyIntToDoubleVec(Vec<double>& B, const Vec<int>& A) {
    
    B.SetLength(A.length());
    
    for(int i = 0; i < B.length(); i++)
        B[i] = (double)A[i];
    
}

/* Generic method for Gaussian Sampling over lattices */
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

RR Samplers::Norm(const vec_RR& b, int n) {
    
    RR norm, mult;
    norm = to_RR(0);

    for(int i = 0; i < n; i++) {
        mul(mult, b[i], b[i]);
        norm += mult;
    }//end-for
    
    return SqrRoot(norm);
        
}//end-Norm()

double Samplers::InnerProduct(const Vec<int>& a, const Vec<int>& b) {
    
    double innerp, mult;
    innerp = 0.0;

    for(int i = 0; i < a.length(); i++) {
        mult = ((double)a[i])*((double)b[i]);
        innerp += mult;
    }//end-for
    
    return innerp;
    
}//end-InnerProduct()

double Samplers::InnerProduct(const Vec<double>& a, const Vec<double>& b) {
    
    double innerp = 0.0;

    for(int i = 0; i < a.length(); i++)
        innerp += a[i]*b[i];
    
    return innerp;
    
}//end-InnerProduct()

double Samplers::InnerProduct(const Vec<int>& a, const Vec<double>& b) {
    
    double innerp, mult;
    innerp = 0.0;

    for(int i = 0; i < a.length(); i++) {
        mult = ((double)a[i])*b[i];
        innerp += mult;
    }//end-for
    
    return innerp;
    
}//end-InnerProduct()

RR Samplers::InnerProduct(const vec_RR& a, const vec_RR& b, int n) {
    
    RR innerp, mult;
    innerp = to_RR(0);

    for(int i = 0; i < n; i++) {
        mul(mult, a[i], b[i]);
        innerp += mult;
    }//end-for
    
    return innerp;
    
}//end-InnerProduct()

double Samplers::NormOfBasis(const Vec< Vec<double> >& T_ATilde) {
    
    double norm, normT_ATilde = 0.0;
    
    for(int i = 0; i < T_ATilde.length(); i++) {
        norm = sqrt(this->InnerProduct(T_ATilde[i], T_ATilde[i]));
        if(norm > normT_ATilde)
            normT_ATilde = norm;
    }//end-for
    
    return normT_ATilde;
    
}//end-NormOfBasis()

double Samplers::NormOfBasis(const Vec< Vec<int> >& T_A) {
    
    double norm, normT_A = 0.0;
    
    for(int i = 0; i < T_A.length(); i++) {
        norm = sqrt(this->InnerProduct(T_A[i], T_A[i]));
        if(norm > normT_A)
            normT_A = norm;
    }//end-for
    
    return normT_A;
    
}//end-NormOfBasis()

void Samplers::PrintMatrix(const string& name, const Vec< Vec<int> >& matrix) {
    
    cout << "\n/** " << name << " **/" << endl;
    for(int i = 0; i < matrix.length(); i++) {
        for(int j = 0; j < matrix[0].length(); j++)
            cout << matrix[i][j] << " ";
        cout << endl;
    }//end-for
    
}//end-PrintVectorZZX()