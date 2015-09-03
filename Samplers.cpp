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
    
    /* Create the partitioning first */
    RR v; // Area of each rectangle
    this->ZCreatePartition(m, sigma, precision, tail, v);
    
    /* Sampling phase */    
    RR c, U, z;
    int i;
    double ulong_size = (double)(sizeof(unsigned long)*8);
    
    c = to_RR(0);
    
    for(;;) {
    
        i = RandomBnd(m); // Sample a random value in {0, ..., m-1}
        U = to_RR((RandomBnd(ulong_size-1) + 1)/ulong_size); // Uniform sample in (0, 1)
        
        if(i > 0)
            z = U*this->X[i];
        else
            z = U*v/this->Probability(X[i], sigma, c);
        
        if(z > this->X[i+1])
            return z;
        
        if(i == 0) // Sampling from the tail
            return this->NewMarsagliaTailMethod(this->X[m-1]);
        
        U = to_RR((RandomBnd(ulong_size-1) + 1)/ulong_size);
        
        if(i > 0 && U*(this->Probability(X[i], sigma, c) - this->Probability(this->X[i+1], sigma, c)) < (this->Probability(z, sigma, c) - this->Probability(this->X[i+1], sigma, c)))
            return z;
            
    }//end-for
    
    cout << "Pass!\n";
    
}//end-Ziggurat()

/* DZCreatePartition defines the x and y axes of rectangles in the Gaussian distribution */
void Samplers::ZCreatePartition(int m, RR sigma, int n, RR tail, RR& v) {
    
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
    r = (tail*sigma)/2; //r = x_m, the m-th x-value
    bestdiff = r; // Arbitrary initial value
    
    while(z != zero && r > zero) { // If r = 0 then the area v is also 0 (what doesn't make sense...)

        z = this->ZRecursion(X, m, r, sigma, v);        

        if(z == minusOne && first) { // Error in "inv" or square root computation
            
            first = 0;
            add(r, r, 2*overM);
            div(overM, overM, to_RR(m));
            
            while(z != zero && r > zero) { // If r = 0 then the area v is also 0 (what doesn't make sense...)
                
                z = this->ZRecursion(X, m, r, sigma, v);

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
        
        this->X.SetLength(m);
        for(i = 0; i < this->X.length(); i++) {
            this->X[i] = bestX[i];
        }//end-for
        
    } else // No valid partition was found
        cout << "Error!" << endl;
            
}//end-DZCreatePartition()

/* It is used in DZCreatePartition to define the distance y0 */
RR Samplers::ZRecursion(Vec<RR>& X, int m, RR r, RR sigma, RR& v) {
    
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

RR Samplers::GramSchmidtProcess(mat_RR& T_ATilde, const mat_RR& T_A, long precision) {
    
    RR::SetPrecision(precision);
    
    cout << "\n[!] Norm of the short basis: " << this->NormOfBasis(T_A) << endl;    
    cout << "[*] Gram-Schmidt process status: ";
    
    vec_RR mult;
    RR mu, norm;
    int i, j;
    int cols, rows;
    
    cols = T_A.NumCols();
    rows = T_A.NumRows();
    
    T_ATilde.SetDims(rows, cols);
                
    for(i = 0; i < rows; i++) { // For each row vector
        
        T_ATilde[i] = T_A[i];
        
        for(j = 0; j < i; j++) {
            div(mu, this->InnerProduct(T_A[i], T_ATilde[j]), this->InnerProduct(T_ATilde[j], T_ATilde[j]));
            cout << mu << " ";
            mul(mult, T_ATilde[j], mu);
            sub(T_ATilde[i], T_ATilde[i], mult);
        }//end-for
        
    }//end-for
    
    norm = this->NormOfBasis(T_ATilde);
    this->TTilde = T_ATilde;
    cout << "Pass!" << endl;
    
    return norm;
    
}//end-GramSchmidtProcess() 

double Samplers::GramSchmidtProcess(Vec< Vec<double> >& T_ATilde, const Vec< Vec<int> >& T_A) {
    
    /*
     * Input: a square matrix T_A (nm x nm)
     */
    
    cout << "\n[!] Norm of the short basis: " << this->NormOfBasis(T_A) << endl;
    
    cout << "[*] Gram-Schmidt process status: ";
    
    int i, j, k;
    int length;
    double mu, norm;
    
    length = T_A.length();
    
    T_ATilde.SetLength(length);
    for(i = 0; i < length; i++)
        T_ATilde[i].SetLength(length);
                
    for(i = 0; i < length; i++) { // For each vector
        
        this->CopyIntToDoubleVec(T_ATilde[i], T_A[i]);
        
        for(j = 0; j < i; j++) {
            
            mu = this->InnerProduct(T_A[i], T_ATilde[j])/this->InnerProduct(T_ATilde[j], T_ATilde[j]);
            
            for(k = 0; k < length; k++)
                T_ATilde[i][k] = T_ATilde[i][k] - mu*T_ATilde[j][k];            
            
        }//end-for
    }//end-for
    
    norm = this->NormOfBasis(T_ATilde);
    
    cout << "Pass!" << endl;
    
    return norm;
    
}//end-GramSchmidtProcess() 

void Samplers::CopyIntToDoubleVec(Vec<double>& B, const Vec<int>& A) {
    
    B.SetLength(A.length());
    
    for(int i = 0; i < B.length(); i++)
        B[i] = (double)A[i];
    
}//end-CopyIntToDoubleVec()

/* Computes the decomposition of A into B*B^t and returns the lower-triangular matrix B */
void Samplers::CholeskyDecomposition(Vec< Vec<double> >& B, const Vec< Vec<double> >& A, int n) {
    
    int i, j;
    
    B.SetLength(n);
    for(i = 0; i < n; i++) {
        B[i].SetLength(n);
        for(j = 0; j < n; j++) {
            B[i][j] = 0.0;
            cout << A[i][j] << " ";
        }
        cout << endl;
    }//end-for
    
    int k;
    double s;
    
    for(i = 0; i < n; i++) { // Every line
        for(j = 0; j <= i; j++) { // Only columns before diagonal
            
            s = 0.0;
            
            for(k = 0; k < j; k++) // Sum
                s += B[i][k]*B[j][k];
            
            cout << "s: " << s << endl;
            if(i == j) {
                cout << "Caso i == j: " << sqrt(A[i][i] - s) << endl;
                B[i][j] = sqrt(A[i][i] - s);                    
            }
            else
                B[i][j] = 1.0/B[j][j]*(A[i][j] - s);
            cout << B[i][j] << " ";
        }//end-for
        cout << endl;
    }
    
}//end-CholeskyDecomposition()

RR Samplers::BlockGSO(mat_RR& BTilde, const Vec<ZZX>& B, int n, int precision) {
    
    // Precision of floating point operations
    RR::SetPrecision(precision);    
    
    cout << "[!] Norm of basis B: " << this->NormOfBasis(B, B.length(), n) << endl; 
    cout << "\n[*] Block_GSO status: ";

    mat_RR outRot, outGSO; // (n x n)-dimension matrices
    vec_RR conv, ortho, mult;
    RR mu, innerpr, norm;
    int i, j, k;
    
    BTilde.SetDims(B.length(), n);
    
    k = B.length()/n; // Number of isometric blocks
    
    conv.SetLength(n);
    ortho.SetLength(n);
    mult.SetLength(n);
    
    for(i = 0; i < k; i++) { // For each isometric block i in [1, ..., k]

        for(j = 0; j < n; j++) {
            conv[j] = to_RR(B[i*n][j]); // Copy the vector whose expansion is an isometric basis
            ortho[j] = to_RR(0);
        }//end-for

         /* For each expansion of a vector makes B[i*n] orthogonal 
          * to the previous vectors */
        for(j = 0; j < i*n; j++) {
            if(!IsZero(BTilde[j])) {
                innerpr = this->InnerProduct(conv, BTilde[j]);            
                norm = this->InnerProduct(BTilde[j], BTilde[j]);            
                div(mu, innerpr, norm);           
                mul(mult, BTilde[j], mu);
                add(ortho, ortho, mult);
            }//end-if
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
    
    norm = this->NormOfBasis(BTilde);
    
    return norm;
    
}//end-BlockGSO()

void Samplers::FasterIsometricGSO(mat_RR& BTilde, const mat_RR& B) {
    
    /* This implementation follows the Algorithm 3 
     * of (Lyubashevsky, and Prest, 2015) */
    if(IsZero(B)) {
//        cout << "[!] Warning: the input basis is null." << endl;
        BTilde = B;
        return;
    }//end-if
    
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
    C[0] = this->InnerProduct(V, this->Isometry(BTilde[0], n));
    D[0] = this->InnerProduct(BTilde[0], BTilde[0]);
    
    for(int i = 0; i < n-1; i++) {        
        div(CD, C[i], D[i]);        
        isometry = this->Isometry(BTilde[i], n);         
        mul(mult, V, CD);                
        sub(BTilde[i+1], isometry, mult);            
        mul(mult, isometry, CD);        
        sub(V, V, mult);                        
        C[i+1] = this->InnerProduct(B[0], this->Isometry(BTilde[i+1], n));         
        sub(D[i+1], D[i], CD*C[i]);                
    }//end-for
    
}//end-FasterIsometricGSO()

/* Giving a polynomial g, out contains (b*x)%f(x) */
ZZX Samplers::Isometry(ZZX& b, int n) {
    
    ZZX out;
    out.SetLength(n);
    
    if(IsZero(b))
        out = to_ZZX(0);        
    else {
        out[0] = -b[n-1];    
        for(int i = 1; i < n; i++)
            out[i] = b[i-1];
    }//end-else
    
    return out;
    
}//end-Isometry()

vec_RR Samplers::Isometry(vec_RR& b, int n) {
    
    vec_RR out;    
    out.SetLength(n);
        
    out[0] = -b[n-1];    
    for(int i = 1; i < n; i++)
        out[i] = b[i-1];
    
    return out;
    
}//end-Isometry()

void Samplers::RotBasis(Vec<ZZX>& T, const Vec< Vec<ZZX> >& S, int n) {
    
    Vec< Vec<ZZX> > outInnerRot;
    int i, j, k, l;
    int m = S.length();
    
    T.SetLength(m*n);
    for(i = 0; i < (m*n); i++)
        T[i].SetLength(m*n);
    
    for(i = 0; i < m; i++) {
        this->Rot(outInnerRot, S[i], m, n);        
        for(j = 0; j < n; j++) // For each row
            T[i*n + j] = outInnerRot[j];
    }//end-for    
    
}//end-RotBasis()

void Samplers::RotBasis(Vec< Vec<int> >& T, const Vec< Vec<ZZX> >& S, int n) {
    
    /*
     * Input:
     * A (m x m)-dimension matrix of polynomials in R_0
     * Output:
     * A (nm x nm)-dimension integer matrix
     */
        
    Vec< Vec<ZZX> > outInnerRot;
    int i, j, k, l;
    int m = S.length();
    
    T.SetLength(m*n);
    for(i = 0; i < (m*n); i++)
        T[i].SetLength(m*n);
    
    for(i = 0; i < m; i++) {
        this->Rot(outInnerRot, S[i], m, n);        
        for(j = 0; j < n; j++) // For each row
            for(k = 0; k < m; k++) // For each column
                for(l = 0; l < n; l++) // For each coefficient in S[i][j]
                    T[i*n+j][k*n+l] = to_int(outInnerRot[j][k][l]);                        
    }//end-for    
    
}//end-RotBasis()

void Samplers::RotBasis(Vec<ZZX>& T, const Vec< Vec<ZZX> >& S, int n) {
    
    Vec<ZZX> outrot;
    int i, j, k;
    int m = S.length();
    
    T.SetLength(m*m*n);
    for(i = 0; i < (m*m*n); i++)
        T[i].SetLength(n);
    
    int index = 0;
    for(i = 0; i < m; i++) {
        for(j = 0; j < m; j++) {
            this->rot(outrot, S[i][j], n);
            for(k = 0; k < n; k++)
                T[index+k] = outrot[k];
            index += n;
        }//end-for
    }//end-for    
        
}//end-RotBasis()

/* Applies the rot operator component-wise in a row vector a */
void Samplers::Rot(Vec< Vec<ZZX> >& A, const Vec<ZZX>& a, int m, int n) {
    
    Vec<ZZX> out;
    int i, j;
    
    A.SetLength(n);    
    for(i = 0; i < n; i++) {
        A[i].SetLength(m);
        for(j = 0; j < m; j++)
            A[i][j].SetLength(n);
    }//end-for
    
    for(j = 0; j < m; j++) {
        this->rot(out, a[j], n);
        for(i = 0; i < n; i++) {
            A[i][j] = out[i];
            A[i][j].normalize();
        }//end-for
    }//end-for
    
}//end-Rot()

void Samplers::rot(Vec<ZZX>& out, const ZZX& b, int n) {
    
    ZZX isometry;   
    out.SetLength(n); 
    
    out[0].SetLength(n);
    isometry.SetLength(n);
    
    out[0] = isometry = b;
    
    for(int i = 1; i < n; i++) {
        isometry = this->Isometry(isometry, n);
        out.SetLength(n);
        out[i] = isometry;
    }//end-for
   
}//end-rot()

void Samplers::rot(mat_RR& out, const vec_RR& b, int n) {
        
    vec_RR isometry;
    isometry.SetLength(n);
    
    out.SetDims(n, n);    
    
    out[0] = b;    
    isometry = b;    
    
    for(int i = 1; i < n; i++) {
        isometry = this->Isometry(isometry, n);
        out[i] = isometry;
    }//end-for
    
}//end-rot()

void Samplers::SetCenter(vec_RR& c, const Vec< Vec<int> >& S) {
    
    int cols, rows;    
    cols = S[0].length();
    rows = S.length();
    
    c.SetLength(rows);
    RR acc = to_RR(0);
    
    for(int i = 0; i < rows; i++) {
        acc = to_RR(0);
        for(int j = 0; j < cols; j++)
            acc += to_RR(S[i][j]*S[i][j]);
        c[i] = sqrt(acc);
    }//end-for
    
}//end-SetCenter()


vec_RR Samplers::CompactGaussianSampler(const Vec< Vec<int> >& B, RR sigma, const vec_RR center, const vec_RR& BTilden, const vec_RR& Vn, const vec_RR& H, const vec_RR& I) {
    
}

vec_RR Samplers::GaussianSamplerFromLattice(const Vec< Vec<int> >& B, const mat_RR& BTilde, RR sigma, int precision, int tailcut, const vec_RR center) {

    cout << "\n[*] Gaussian Sampler status: ";
    
    RR::SetPrecision(precision);    

    vec_RR C, sample;
    int Z;
    
    RR d, innerp, sigma_i;
    int cols, i, j, rows;
    
    cols = BTilde.NumCols();
    rows = BTilde.NumRows();
    
    C.SetLength(cols);
    sample.SetLength(cols);
                
    C = center; // Center of the lattice        
    
    /* The inner product and norm operations are taken 
     * in the inner product space H */
    for(i = rows-1; i >= 0; i--) {
        
        innerp = this->InnerProduct(BTilde[i], BTilde[i]);
        
        // The new center for the discrete Gaussian
        div(d, this->InnerProduct(C, BTilde[i]), innerp);        
        
        // And the new standard deviation
        div(sigma_i, sigma, sqrt(innerp));
        
        cout << "sigma_i: " << sigma_i << endl;
        cout << "t*sigma: " << to_RR(tailcut)*sigma_i << endl;
        
        this->BuildProbabilityMatrix(precision, tailcut, sigma_i, d);             
        
        Z = this->KnuthYao(tailcut, sigma_i, d);
        
        for(j = 0; j < B.length(); j++)
            C[j] = C[j] - to_RR(B[i][j]*Z);
        
    }//end-for
    
    sub(sample, center, C);
    
    cout << "Pass!" << endl;
    
    return sample;
    
}//end-GaussianSamplerFromLattice()

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
        div(d, this->InnerProduct(auxC, BTilde[i]), this->InnerProduct(BTilde[i], BTilde[i]));        
        
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
        add(norm, norm, mult);
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

RR Samplers::InnerProduct(const vec_RR& a, const vec_RR& b) {
    
    RR innerp, mult;
    innerp = to_RR(0);
    
    for(int i = 0; i < a.length(); i++) {
        mul(mult, a[i], b[i]);
        innerp += mult;
    }//end-for
    
    return innerp;
    
}//end-InnerProduct()

RR Samplers::NormOfBasis(const mat_RR& B) {
    
    RR norm, normB;    
    
    normB = to_RR(0);
    
    for(int i = 0; i < B.NumRows(); i++) {
        norm = this->Norm(B[i], B.NumCols());
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

/* Norm of a polynomial */
RR Samplers::Norm(const ZZX& b, int n) {
    
    ZZ norm = to_ZZ(0);

    for(int i = 0; i < n; i++)
        norm += b[i]*b[i];
    
    return SqrRoot(to_RR(norm));
        
}//end-Norm()

void Samplers::PrintMatrix(const string& name, const Vec< Vec<int> >& matrix) {
    
    cout << "\n/** " << name << " **/" << endl;
    for(int i = 0; i < matrix.length(); i++) {
        for(int j = 0; j < matrix[0].length(); j++)
            cout << matrix[i][j] << " ";
        cout << endl;
    }//end-for
    
}//end-PrintVectorZZX()