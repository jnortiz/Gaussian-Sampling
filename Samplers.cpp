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