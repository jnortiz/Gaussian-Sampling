/* 
 * File:   Samplers.h
 * Author: jnortiz
 *
 * Created on April 24, 2015, 3:51 PM
 */

#include <NTL/mat_ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>
#include <NTL/RR.h>

#ifndef SAMPLERS_H
#define	SAMPLERS_H

using namespace std;
using namespace NTL;

class Samplers {
public:
       
    Samplers(int q, const ZZ_pX& f);
    virtual ~Samplers();
        
    vec_RR CompactGaussianSampler(const mat_RR& B, const vec_RR center, const vec_RR& BTildeN, const vec_RR& Vn, const vec_RR& D, const vec_RR& I, RR sigma, long precision);
    void PrepareToSampleCGS(vec_RR& v, vec_RR& I, const mat_RR& BTilde, const vec_RR& D, const vec_RR& B1, long precision);
    
    /* Generic algorithm for GSO basis construction */
    RR GramSchmidtProcess(mat_RR& T_ATilde, const mat_RR& T_A, long precision);
    RR GramSchmidtProcess(mat_RR& T_ATilde, vec_RR& D, const mat_RR& T_A, long precision);
    
    /* Sampling from the discrete Gaussian distribution D_{\Lambda, \sigma, c}*/
    vec_RR GaussianSamplerFromLattice(const mat_ZZ& B, const mat_RR& BTilde, RR sigma, int precision, int tailcut, const vec_RR center);
    
    /* Klein's method */
    vec_RR Klein(const mat_ZZ& B, const mat_RR& BTilde, RR sigma, long precision, int tailcut, const vec_RR center);
    
    RR ZCreatePartition(int m, RR sigma, long n, RR tail);
    
    /* Peikert's method */
    RR NormOfBasis(const mat_RR& B);
    int OfflinePeikert(mat_ZZ& Z, mat_RR& B2, const mat_ZZ& B, int q, RR r, const mat_RR& Sigma, int n, long precision);
    vec_ZZ RefreshPeikert(const mat_RR& B2, RR r, RR v, int n, long precision);
    vec_ZZ Peikert(const mat_ZZ_p& B, const mat_ZZ_p Z, const vec_ZZ_p& c, const vec_ZZ_p& x2, long q, RR r, long precision);
    
    void RotBasis(mat_ZZ& T, const Vec< Vec<ZZX> >& S, int n);
    
private:
        
    /* Attributes for sampling from lattice */
    ZZ_pX f;// Polynomial R = Z_p[X]/f    
    
    /* (Continuous) Ziggurat attributes */
    Vec<RR> X;
    
    /* Knuth-Yao attributes */
    Vec< Vec<int> > P;
    Vec<int> begin;
    
    /* Sampling from a discrete Gaussian distribution over the integers */
    int KnuthYao(int tailcut, RR sigma, RR c);
    void BuildProbabilityMatrix(int precision, int tailcut, RR sigma, RR c);
    RR Probability(RR x, RR sigma, RR c);
    
    /* Sampling from a continuous Gaussian distribution */
    RR Ziggurat(int m, RR sigma, long precision, RR v);
    RR ZRecursion(Vec<RR>& X, int m, RR r, RR sigma, RR& v);
    RR NewMarsagliaTailMethod(RR r);
        
    /* Gram-Schmidt orthogonalization */
    RR BlockGSO(mat_RR& BTilde, const Vec<ZZX>& B, int n, int precision);
    void FasterIsometricGSO(mat_RR& BTilde, const mat_RR& B);
    
    /* Sampling from a discrete Gaussian distribution over lattices */
    ZZX GaussianSamplerFromLattice(const Vec<ZZX>& B, const mat_RR& BTilde, RR sigma, int precision, int tailcut, ZZX center, int n);
    
    /* Procedure for T = BB^t used in Peikert algorithm - Lattice sampling */
    void SetCenter(vec_RR& c, const mat_ZZ& S);
    int CholeskyDecomposition(mat_RR& B, const mat_RR& A, int n);    
        
    void Rot(Vec< Vec<ZZX> >& A, const Vec<ZZX>& a, int m, int n); // RotBasis()
    void rot(Vec<ZZX>& out, const ZZX& b, int n); // Rot()
    ZZX Isometry(ZZX& b, int n);    
    
    void rot(mat_RR& out, const vec_RR& b, int n); // BlockGSO() and FasterIsometricGSO()
    vec_RR Isometry(const vec_RR& b);
    vec_RR InvIsometry(const vec_RR& b);
    
    RR Norm(const vec_RR& b, int n);
    
    RR NormOfBasis(const Vec<ZZX>& B, int m, int n);
    RR Norm(const ZZX& b, int n);

    void PrintMatrix(const string& name, const Vec< Vec<int> >& matrix);
    
};

#endif	/* SAMPLERS_H */

