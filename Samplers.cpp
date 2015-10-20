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
#include <NTL/ZZX.h>

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

Samplers::~Samplers() {
    this->P.kill();
    this->begin.kill();
    this->X.kill();
    this->f.kill();
}

RR Samplers::Probability(RR x, RR sigma, RR c) {
    
    RR S = sigma*sqrt(2*ComputePi_RR());
    RR overS = 1/S;
    
    if(x == to_RR(0))
        return overS;
    
    return overS*exp(-(power((x-c)/sigma, 2))/2.0);
    
}//end-Probability()

// If the bit is zero, the output becomes "a" 
int Select(int a, int b, unsigned bit) {
    
    unsigned mask;
    int output;
    
    mask = -bit;
    output = mask & (a ^ b);
    output = output ^ a;
    
    return output;
    
}//end-Select()

int Samplers::KnuthYao(int tailcut, RR sigma, RR c) {

    int bound, center, col, d, invalidSample, pNumRows, pNumCols, S, signal;
    unsigned enable, hit;
    unsigned long r;
    
    bound = tailcut*to_int(sigma);
    center = to_int(c);
    d = 0;
    hit = 0;
    signal = 1 - 2*RandomBits_long(1);
    invalidSample = bound+1;
    pNumRows = this->P.length();
    pNumCols = this->P[0].length();    
    
    Vec<int> randomBits;
    randomBits.SetLength(pNumRows);
    
    int i, index, j, length;
    length = sizeof(unsigned long)*8; 
    
    index = 0;
    for(i = 0; i < (pNumRows/length+1); i++) {
        r = RandomWord();
        for(j = 0; j < length, index < pNumRows; j++, r >>= 1)
            randomBits[index++] = (r & 1);
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

RR Samplers::Ziggurat(int m, RR sigma, long precision, RR v) {

/*
 * Important: run ZCreatePartition() procedure before calling the Ziggurat sampler.
 */
    
#ifdef DEBUG
    cout << "[*] Ziggurat status: ";
#endif
    
    RR::SetPrecision(precision);
    
    RR c, probX, probX1, probZ, U, z;
    int i;
    double ulong_size = (double)(sizeof(unsigned long)*8);
    
    c = to_RR(0);
    
    for(;;) {
    
        i = NTL::RandomBnd((long)(m-1)) + 1;
        U = to_RR((NTL::RandomBnd((long)(ulong_size-1)) + 1)/ulong_size); // Uniform sample in (0, 1)
        
        if(i < (m-1))
            NTL::mul(z, U, this->X[i]);
        else // Base strip
            NTL::div(z, U*v, this->Probability(X[i], sigma, c));
            
        if(z < this->X[i-1]) // The sample is in the left-most part of the i-th rectangle
            return z;

        if(i == (m-1))
            return this->NewMarsagliaTailMethod(this->X[m-1]);
        
        probX = this->Probability(X[i], sigma, c);
        probX1 = this->Probability(this->X[i+1], sigma, c);
        probZ = this->Probability(z, sigma, c);

        U = to_RR((RandomBnd(ulong_size-1) + 1)/ulong_size);

        // The sample can be in the right-most part of the i-th rectangle        
        if(i > 0 && U*(probX - probX1) < (probZ - probX1))
            return z;
            
    }//end-for
    
#ifdef DEBUG
    cout << "Pass!" << endl;
#endif
    
}//end-Ziggurat()

/* DZCreatePartition defines the x- and y-axes of the rectangles in the Gaussian distribution */
RR Samplers::ZCreatePartition(int m, RR sigma, long n, RR tail) {
    
    RR::SetPrecision(n);
    
#ifdef DEBUG
    cout << "\n[*] ZCreatePartition status: ";
#endif
    
    /* The Ziggurat algorithm was designed for centered Gaussian distributions; 
     i.e., c = 0. */
    
    /*
     * Parameters description:
     * m: number of rectangles
     * sigma: Gaussian distribution parameter
     * n: bit precision
     */
    
    Vec<RR> bestX, X;
    RR bestdiff, r, v, z;
    
    RR overM = to_RR(1)/to_RR(m);
    RR minusOne = to_RR(-1);
    RR zero = to_RR(0);
    
    int i;
    int first = 1;
    
    X.SetLength(m);
    bestX.SetLength(m);            
    
    bestX[m-1] = minusOne;    
    z = minusOne;    
    r = (tail*sigma)/2; // r = x_m, the m-th x-value
    v = 0;
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
        
#ifdef DEBUG
        cout << "Pass!" << endl;        
        cout << "\nFinal z value: " << bestdiff << endl;
        cout << "Final r value: " << X[m-1] << endl;
        cout << "Statistical difference: " << to_RR(power2_RR(-((long)n)) - bestdiff) << endl;
#endif
        
        this->X.SetLength(m);
        for(i = 0; i < this->X.length(); i++) {
            this->X[i] = bestX[i];
        }//end-for
        
    } else // No valid partition was found
#ifdef DEBUG
        cout << "Error!" << endl;
#endif
    
    return v;
    
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

// Usual algorithm for Gram-Schmidt orthogonalization
RR Samplers::GramSchmidtProcess(mat_RR& T_ATilde, const mat_RR& T_A, long precision) {
    
    RR::SetPrecision(precision);
    
    cout << "\n[!] Norm of short basis T: " << this->NormOfBasis(T_A) << endl;    
    cout << "[*] Gram-Schmidt process status: ";
    
    mat_RR mu;
    vec_RR innerpT_ATilde, mult, sum;
    RR inner1;
    int cols, rows;
    
    cols = T_A.NumCols();
    rows = T_A.NumRows();
    
    mu.SetDims(rows, rows);
    T_ATilde.SetDims(rows, cols);
    mult.SetLength(cols);
    sum.SetLength(cols);
    innerpT_ATilde.SetLength(rows);
    
    for(int i = 0; i < rows; i++)            
        mu[i][i] = to_RR(1);
    
    T_ATilde[0] = T_A[0];

    for(int i = 1; i < rows; i++) {
        
        clear(sum);
        
        for(int j = 0; j < (i-1); j++) {
            NTL::InnerProduct(inner1, T_A[i], T_ATilde[j]);
            div(mu[i][j], inner1, innerpT_ATilde[j]);
            mul(mult, T_ATilde[j], mu[i][j]);
            add(sum, sum, mult);
        }//end-for

        NTL::InnerProduct(inner1, T_A[i], T_ATilde[i-1]);
        NTL::InnerProduct(innerpT_ATilde[i-1], T_ATilde[i-1], T_ATilde[i-1]);
        div(mu[i][i-1], inner1, innerpT_ATilde[i-1]);
        mul(mult, T_ATilde[i-1], mu[i][i-1]);
        add(sum, sum, mult);
        
        sub(T_ATilde[i], T_A[i], sum);
        
    }//end-for
    
    mu.kill();
    innerpT_ATilde.kill();
    mult.kill();
    sum.kill();
    
    cout << "Pass!" << endl;
    RR norm = this->NormOfBasis(T_ATilde);
    return norm;
    
}//end-GramSchmidtProcess() 

// Gram-Schmidt orthogonalization procedure for block isometric basis
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
                NTL::InnerProduct(innerpr, conv, BTilde[j]);
                NTL::InnerProduct(norm, BTilde[j], BTilde[j]); // Square of norm is the inner product
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

// Gram-Schmidt orthogonalization procedure for isometric basis
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
    
    NTL::InnerProduct(C[0], V, this->Isometry(BTilde[0]));
    NTL::InnerProduct(D[0], BTilde[0], BTilde[0]);
    
    for(int i = 0; i < n-1; i++) {        
        div(CD, C[i], D[i]);        
        isometry = this->Isometry(BTilde[i]);         
        mul(mult, V, CD);                
        sub(BTilde[i+1], isometry, mult);            
        mul(mult, isometry, CD);        
        sub(V, V, mult);
        NTL::InnerProduct(C[i+1], B[0], this->Isometry(BTilde[i+1]));
        sub(D[i+1], D[i], CD*C[i]);                
    }//end-for
    
}//end-FasterIsometricGSO()

vec_RR Samplers::GaussianSamplerFromLattice(const mat_ZZ& B, const mat_RR& BTilde, RR sigma, int precision, int tailcut, const vec_RR center) {

    cout << "\n[*] Gaussian Sampler status: ";
    
    RR::SetPrecision(precision);    

    double sizeof_RR = pow(2.0, sizeof(RR));
    RR d, innerp, innerp1, sigma_i, Z;
    vec_RR C, sample;    
    int cols, i, j, rows;
    
    cols = BTilde.NumCols();
    rows = BTilde.NumRows();
    
    C.SetLength(cols);
    sample.SetLength(cols);
                
    C = center;
    
    for(i = rows-1; i >= 0; i--) {
        
        NTL::InnerProduct(innerp, BTilde[i], BTilde[i]);
        NTL::InnerProduct(innerp1, C, BTilde[i]);        
        div(d, innerp1, innerp);        
        div(sigma_i, sigma, sqrt(innerp));
        
        if(floor(sigma_i) == 0)
            Z = d;
        else {
            if(sigma_i > sizeof_RR)
                sigma_i = to_RR(2)*sigma;
            this->BuildProbabilityMatrix(precision, tailcut, sigma_i, d);                     
            Z = to_RR(this->KnuthYao(tailcut, sigma_i, d));            
        }//end-else
        
        for(j = 0; j < B.NumCols(); j++)
            C[j] = C[j] - to_RR(B[i][j])*Z;
        
    }//end-for
    
    sub(sample, center, C);    
    cout << "Pass!" << endl;
    
    return sample;
    
}//end-GaussianSamplerFromLattice()

vec_RR Samplers::Klein(const mat_ZZ& B, const mat_RR& BTilde, RR sigma, long precision, int tailcut, const vec_RR center) {

    cout << "\n[*] Klein's sampler status: ";
    
    RR::SetPrecision(precision);    

    double sizeof_RR = pow(2.0, sizeof(RR));
    RR aux_B, c_prime, innerp, innerp1, sigma_prime, Z;
    vec_RR c, v;    
    int cols, i, j, rows;
    
    cols = BTilde.NumCols();
    rows = BTilde.NumRows();
    
    c.SetLength(cols);
    v.SetLength(cols);
                   
    clear(v);
    c = center;
    
//    unsigned isZero, lessThan;
//    RR aux_Z, newSigma;
    
    for(i = rows-1; i >= 0; i--) {
        
        NTL::InnerProduct(innerp1, c, BTilde[i]);        
        NTL::InnerProduct(innerp, BTilde[i], BTilde[i]);
        div(c_prime, innerp1, innerp);        
        div(sigma_prime, sigma, sqrt(innerp));
        
//        isZero = IsZero(to_int(to_ZZ(floor(sigma_prime))));
//        lessThan = lessThan(sizeof_RR, sigma_prime);
//        
//        NTL::mul(newSigma, to_RR(2), sigma);
//        this->BuildProbabilityMatrix(precision, tailcut, sigma_prime, c_prime);                     
//        aux_Z = to_RR(this->KnuthYao(tailcut, sigma_prime, c_prime));                    
//        
//        sigma_prime = Select(sigma_prime, newSigma, lessThan);
//        Z = Select(c_prime, aux_Z, isZero);
        
        if(NTL::floor(sigma_prime) == 0)
            Z = c_prime;
        else {
            if(sigma_prime > sizeof_RR)
                NTL::mul(sigma_prime, to_RR(2), sigma);
            this->BuildProbabilityMatrix(precision, tailcut, sigma_prime, c_prime);                     
            Z = to_RR(this->KnuthYao(tailcut, sigma_prime, c_prime));            
        }//end-else
                
        for(j = 0; j < B.NumCols(); j++) {
            aux_B = to_RR(B[i][j]);
            v[j] = v[j] + aux_B*Z;
            c[j] = c[j] - aux_B*Z;
        }//end-for
        
    }//end-for
    
    cout << "Pass!" << endl;
    
    return v;
    
}//end-Klein()

int Samplers::OfflinePeikert(mat_ZZ& Z, mat_RR& B2, const mat_ZZ& B, int q, RR r, const mat_RR& Sigma, int n, long precision) {
    
    cout << "[*] Offline SampleD status: ";
    
    RR::SetPrecision(precision);
    
    RR factor;
    NTL::div(factor, to_RR(1), sqrt(2*NTL::ComputePi_RR()));
    
    ZZ d;
    mat_ZZ invB;    
    NTL::inv(d, invB, B, 0);
    
    mat_RR aux_invB;
    conv(aux_invB, invB);
    invB.kill();
    
    RR over_d;
    NTL::div(over_d, to_RR(1), to_RR(d));
    RR mult = to_RR(q)*over_d;
    
    mat_RR auxZ;   
    NTL::mul(auxZ, aux_invB, mult);
    aux_invB.kill();
    
    Z.SetDims(auxZ.NumRows(), auxZ.NumCols());
    for(int i = 0; i < Z.NumRows(); i++)
        for(int j = 0; j < Z.NumCols(); j++)
            NTL::TruncToZZ(Z[i][j], auxZ[i][j]);    
    auxZ.kill();
        
    RR r_square;
    NTL::mul(r_square, r, r);    
    
    mat_ZZ transposeB, square;
    mat_RR Sigma1, square_RR;    
    NTL::transpose(transposeB, B);   
    NTL::mul(square, B, transposeB);
    transposeB.kill();
    NTL::conv(square_RR, square);
    square.kill();
    NTL::mul(Sigma1, square_RR, to_RR(2)*r_square*factor);
    square_RR.kill();
    
    mat_RR Sigma2;
    NTL::sub(Sigma2, Sigma, Sigma1);
    Sigma1.kill();    
    
    for(int i = 0; i < Sigma2.NumRows(); i++) 
        for(int j = 0; j < Sigma2.NumCols(); j++) {
            NTL::mul(Sigma2[i][j], Sigma2[i][j], factor);
            NTL::sub(Sigma2[i][j], Sigma2[i][j], r_square);
        }//end-for
    int outputCholesky = this->CholeskyDecomposition(B2, Sigma2, n); // Computes the decomposition of Sigma2 into B1.B1^t
    Sigma2.kill();
    
    if(outputCholesky == -1)
        return -1;

    cout << "Pass!" << endl;
    
    return 0;
    
}//end-OfflineSampleD()

vec_ZZ Samplers::RefreshPeikert(const mat_RR& B2, RR r, RR v, int n, long precision) {
    
    RR::SetPrecision(precision);
    
    RR std_deviation;
    NTL::div(std_deviation, to_RR(1), sqrt(2*NTL::ComputePi_RR()));
    
    vec_ZZ x2;
    x2.SetLength(n);
    
    vec_RR w;
    w.SetLength(n);
    
    for(int i = 0; i < w.length(); i++)
        w[i] = this->Ziggurat(64, std_deviation, precision, v);
    
    vec_RR mult;
    NTL::mul(mult, w, B2);
    w.kill();
    
    for(int i = 0; i < n; i++) {
        this->BuildProbabilityMatrix(precision, 13, r, mult[i]); // Standard deviation r, centered in w[i]
        x2[i] = to_ZZ(this->KnuthYao(13, r, mult[i]));
    }//end-for
    mult.kill();

    return x2;
    
}//end-RefreshSampleD()

vec_ZZ Samplers::Peikert(const mat_ZZ_p& B, const mat_ZZ_p Z, const vec_ZZ_p& c, const vec_ZZ_p& x2, long q, RR r, long precision) {
    
    RR::SetPrecision(precision);
    
    vec_ZZ_p subt, mult;
    NTL::sub(subt, c, x2);    
    NTL::mul(mult, Z, subt);
    subt.kill();
    
    vec_ZZ aux_mult;
    NTL::conv(aux_mult, mult);
    mult.kill();
    
    vec_RR center;    
    center.SetLength(aux_mult.length());    
    for(int i = 0; i < aux_mult.length(); i++)
        NTL::div(center[i], to_RR(aux_mult[i]), to_RR(q));
    aux_mult.kill();
    
    vec_ZZ_p rounding;
    rounding.SetLength(center.length());        
    for(int i = 0; i < center.length(); i++) {
        this->BuildProbabilityMatrix(precision, 13, r, center[i]);
        rounding[i] = to_ZZ_p((long)(this->KnuthYao(13, r, center[i])));
    }//end-for
    center.kill();
    
    vec_ZZ_p mult1;
    NTL::mul(mult1, B, rounding);
    rounding.kill();
    
    vec_ZZ_p x;
    NTL::sub(x, c, mult1);
    mult1.kill();
    
    vec_ZZ aux_x;
    NTL::conv(aux_x, x);
    x.kill();
    
    for(int i = 0; i < aux_x.length(); i++)
        if(aux_x[i] > floor(q/2))
            NTL::sub(aux_x[i], aux_x[i], q);
    
    return aux_x;
    
}//end-SampleD()

int Samplers::CholeskyDecomposition(mat_RR& B, const mat_RR& A, int n) {

/*
 * Every symmetric, positive definite matrix A can be decomposed into a product 
 * of a unique lower triangular matrix L and its transpose.
 * 
 * If (A[i][j] - s) < 0, then A is not positive definite matrix or A is 
 * assymmetric, e.g. A is different of its transpose.
 * 
 */

/*
    mat_RR transposeA;
    NTL::transpose(transposeA, A);
    
    for(int i = 0; i < A.NumRows(); i++)
        for(int j = 0; j < A.NumCols(); j++)
            if(A[i][j] != transposeA[i][j]) {
                cout << "[!] The matrix A is not symmetric." << endl;
                transposeA.kill();
                return -1;
            }//end-if
    transposeA.kill();
*/
    
    int i, j, k;
    RR subt, mult, s;
    
    B.SetDims(n, n);
    
    for(i = 0; i < n; i++) {        
        for(j = 0; j <= i; j++) { 
            
            s = to_RR(0);            
            
            for(k = 0; k < j; k++) {
                NTL::mul(mult, B[i][k], B[j][k]);
                NTL::add(s, s, mult);
            }//end-for  
            
            NTL::sub(subt, A[i][j], s);
                        
            if(i == j) {
                if(subt < to_RR(0)) {
                    cout << "[!] The matrix A is not positive definite." << endl;
                    cout << subt << endl;
                    return -1;
                }//end-if
                B[i][j] = sqrt(subt);                
            } else
                NTL::div(B[i][j], subt, B[j][j]);
            
        }//end-for
    }//end-for
    
    return 0;
    
}//end-CholeskyDecomposition()

vec_RR Samplers::CompactGaussianSampler(const mat_RR& B, const vec_RR center, const vec_RR& BTildeN, const vec_RR& Vn, const vec_RR& I, const vec_RR& D, RR sigma, long precision) {

    cout << "\n[*] Compact-Gaussian-Sampler status: ";
    
    RR::SetPrecision(precision);
    
    vec_RR aux_b, b, c, mult, mult1, sample, sum, v;
    RR ci, innerp, H, sigmai;
    double z;
    int cols, i, rows;
    
    cols = B.NumCols();
    rows = B.NumRows();
    
    mult.SetLength(cols);
    mult1.SetLength(cols);
    
    c = center;
    b = BTildeN;
    v = Vn;
    
    for(i = (rows-1); i > 0; i--) {
        
        NTL::InnerProduct(innerp, c, b);
        NTL::div(ci, innerp, D[i]);
        NTL::div(sigmai, sigma, NTL::sqrt(D[i]));
        this->BuildProbabilityMatrix(precision, 13, sigmai, ci);
        z = (double)(this->KnuthYao(13, sigmai, ci));
        
        NTL::mul(mult, B[i], (double)z);
        NTL::sub(c, c, mult);
        
        NTL::div(H, D[i-1], D[i]);
        NTL::mul(mult, v, I[i-1]);
        NTL::mul(mult1, b, H);
        NTL::add(sum, mult1, mult); 
        aux_b = b;
        b = this->InvIsometry(sum);        

        NTL::mul(mult, v, H);
        NTL::mul(mult1, aux_b, I[i-1]);
        NTL::add(v, mult1, mult);                
        
    }//end-for
    
    aux_b.kill();
    mult1.kill();
    sum.kill();
    v.kill();

    NTL::InnerProduct(innerp, c, b);
    NTL::div(ci, innerp, D[i]);
    NTL::div(sigmai, sigma, NTL::sqrt(D[i]));
    this->BuildProbabilityMatrix(precision, 13, sigmai, ci);
    z = (double)(this->KnuthYao(13, sigmai, ci));

    NTL::mul(mult, B[i], (double)z);
    NTL::sub(c, c, mult);    
    
    b.kill();
    mult.kill();
    
    NTL::sub(sample, center, c);
    c.kill();
    
    cout << "Pass!" << endl;
    
    return sample;    
    
}//end-CompactGaussianSampler()

void Samplers::PrepareToSampleCGS(vec_RR& v, vec_RR& I, const mat_RR& BTilde, const vec_RR& D, const vec_RR& B1, long precision) {
    
    /**
     * 
     * @param v
     * @param I
     * @param BTilde - The orthogonal basis of B
     * @param D - It contains the squared norm of each vector in the orthogonal basis ~B
     * @param B1 - The first vector in the short basis
     * 
     */
    
    cout << "\n[*] PrepareToSampleCGS status: ";
    
    RR::SetPrecision(precision);    
        
    vec_RR isometry, aux_v, mult;
    RR C, di;
    int cols, i, rows;

    cols = BTilde.NumCols();
    rows = BTilde.NumRows();
    
    I.SetLength(rows-1);    
    isometry.SetLength(cols);
    aux_v.SetLength(cols);

    aux_v = B1;
    for(i = 0; i < (rows-1); i++) {
        isometry = this->Isometry(BTilde[i]);
        NTL::InnerProduct(C, aux_v, isometry);
        div(I[i], C, D[i+1]);
        div(di, C, D[i]);
        mul(mult, isometry, di);
        sub(aux_v, aux_v, mult);
    }//end-for    
    
    v = aux_v;
    aux_v.kill();
    isometry.kill();
    mult.kill();
    
    cout << "Pass!" << endl;
            
}//end-PrepareToSampleCGS()

RR Samplers::GramSchmidtProcess(mat_RR& T_ATilde, vec_RR& D, const mat_RR& T_A, long precision) {
    
    cout << "[*] Gram-Schmidt process status: ";
    
    RR::SetPrecision(precision);    
    
    vec_RR mult, sum;
    RR innerp, mu, norm;
    int cols, rows;
    
    cols = T_A.NumCols();
    rows = T_A.NumRows();
    
    T_ATilde.SetDims(rows, cols);
    D.SetLength(rows);
    mult.SetLength(cols);
    sum.SetLength(cols);
    
    T_ATilde[0] = T_A[0];

    for(int i = 1; i < rows; i++) {        
        
        clear(sum);
        
        for(int j = 0; j < (i-1); j++) {
            NTL::InnerProduct(innerp, T_A[i], T_ATilde[j]);
            div(mu, innerp, D[j]);
            mul(mult, T_ATilde[j], mu);
            add(sum, sum, mult);
        }//end-for

        NTL::InnerProduct(innerp, T_A[i], T_ATilde[i-1]);
        NTL::InnerProduct(D[i-1], T_ATilde[i-1], T_ATilde[i-1]);
        div(mu, innerp, D[i-1]);
        mul(mult, T_ATilde[i-1], mu);
        add(sum, sum, mult);
        
        sub(T_ATilde[i], T_A[i], sum);
        
    }//end-for
    
    NTL::InnerProduct(D[rows-1], T_ATilde[rows-1], T_ATilde[rows-1]);    
    norm = this->NormOfBasis(T_ATilde);
    
    mult.kill();
    sum.kill();
    
    cout << "Pass!" << endl;
    
    return norm;
    
}//end-GramSchmidtProcess() 

/* Expansion of matrix T_A, which was generated by IdealTrapGen, into the integer basis S */
void Samplers::RotBasis(mat_ZZ& T, const Vec< Vec<ZZX> >& S, int n) {
    
    /*
     * Input:
     * A (m x m)-dimension matrix of polynomials in R_0
     * Output:
     * A (nm x nm)-dimension integer matrix
     */
        
    int i, j, k, l;
    int m = S.length();
    
    T.SetDims(m*n, m*n);
    
    for(i = 0; i < m; i++) {
        Vec< Vec<ZZX> > outInnerRot;
        this->Rot(outInnerRot, S[i], m, n);        
        for(j = 0; j < n; j++) // For each row
            for(k = 0; k < m; k++) // For each column
                for(l = 0; l < n; l++)
                    T[i*n+j][k*n+l] = to_ZZ(outInnerRot[j][k][l]);                        
    }//end-for    
    
}//end-RotBasis()

/* Applies the rot operator component-wise in a row vector a */
void Samplers::Rot(Vec< Vec<ZZX> >& A, const Vec<ZZX>& a, int m, int n) {
    
    int i, j;
    
    A.SetLength(n);    
    for(i = 0; i < n; i++) {
        A[i].SetLength(m);
        for(j = 0; j < m; j++)
            A[i][j].SetLength(n);
    }//end-for
    
    for(j = 0; j < m; j++) {
        Vec<ZZX> out;
        this->rot(out, a[j], n);
        for(i = 0; i < n; i++)
            A[i][j] = out[i];
    }//end-for
    
}//end-Rot()

void Samplers::rot(Vec<ZZX>& out, const ZZX& b, int n) {
    
    ZZX isometry;   
    out.SetLength(n); 
    
    out[0].SetLength(n);
    isometry.SetLength(n);
    
    out[0] = b;
    isometry = b;
    
    for(int i = 1; i < n; i++) {
        out[i].SetLength(n);
        isometry = this->Isometry(isometry, n);
        out[i] = isometry;
    }//end-for
   
}//end-rot()

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

void Samplers::rot(mat_RR& out, const vec_RR& b, int n) {
        
    vec_RR isometry;
    isometry.SetLength(n);
    
    out.SetDims(n, n);    
    
    out[0] = b;    
    isometry = b;    
    
    for(int i = 1; i < n; i++) {
        isometry = this->Isometry(isometry);
        out[i] = isometry;
    }//end-for
    
}//end-rot()

vec_RR Samplers::Isometry(const vec_RR& b) {
    
    int n = b.length();    
    vec_RR out;    
    out.SetLength(n);
    
        
    out[0] = -b[n-1];    
    for(int i = 1; i < n; i++)
        out[i] = b[i-1];
    
    return out;
    
}//end-Isometry()

vec_RR Samplers::InvIsometry(const vec_RR& b) {
    
    int n = b.length();    
    vec_RR out;    
    out.SetLength(n);
    
        
    out[n-1] = -b[0];    
    for(int i = 1; i < n; i++)
        out[i-1] = b[i];
    
    return out;
    
}//end-InvIsometry()

void Samplers::SetCenter(vec_RR& c, const mat_ZZ& S) {
    
    int cols, rows;    
    cols = S.NumCols();
    rows = S.NumRows();
    
    c.SetLength(rows);
    RR acc = to_RR(0);
    
    for(int i = 0; i < rows; i++) {
        acc = to_RR(0);
        for(int j = 0; j < cols; j++)
            acc += to_RR(S[i][j]*S[i][j]);
        c[i] = sqrt(acc);
    }//end-for
    
}//end-SetCenter()

RR Samplers::NormOfBasis(const mat_RR& B) {
    
    RR innerp, norm, normB;    
    
    normB = to_RR(0);
    
    for(int i = 0; i < B.NumRows(); i++) {
        NTL::InnerProduct(innerp, B[i], B[i]);
        norm = sqrt(innerp);
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

void Samplers::PrintMatrix(const string& name, const Vec< Vec<int> >& matrix) {
    
    cout << "\n/** " << name << " **/" << endl;
    for(int i = 0; i < matrix.length(); i++) {
        for(int j = 0; j < matrix[0].length(); j++)
            cout << matrix[i][j] << " ";
        cout << endl;
    }//end-for
    
}//end-PrintVectorZZX()
