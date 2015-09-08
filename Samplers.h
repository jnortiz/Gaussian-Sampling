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

#include <complex>

#ifndef SAMPLERS_H
#define	SAMPLERS_H

using namespace std;
using namespace NTL;

class Samplers {
public:
       
    Samplers(int q, const ZZ_pX& f);
    virtual ~Samplers();
    
    /* The center of distribution is the norm of public basis */
    void SetCenter(vec_RR& c, const mat_ZZ& S);
    
    /* Generic algorithm for GSO basis construction */
    RR GramSchmidtProcess(mat_RR& T_ATilde, const mat_RR& T_A, long precision);
    
    /* Sampling from the discrete Gaussian distribution D_{\Lambda, \sigma, c}*/
    vec_RR GaussianSamplerFromLattice(const mat_ZZ& B, const mat_RR& BTilde, RR sigma, int precision, int tailcut, const vec_RR center);
    
    /* Continuous sampling */
    RR Ziggurat(int m, RR sigma, int precision, RR tail);
    
    void RotBasis(mat_ZZ& T, const Vec< Vec<ZZX> >& S, int n);
    void RotBasis(Vec<ZZX>& T, const Vec< Vec<ZZX> >& S, int n);
    
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
    
    /* Auxiliary functions of Knuth-Yao algorithm */
    void BuildProbabilityMatrix(int precision, int tailcut, RR sigma, RR c);
    void BinaryExpansion(Vec< Vec<int> >& auxP, RR probability, int precision, int index);
    RR Probability(RR x, RR sigma, RR c);

    /* Auxiliary functions of Ziggurat algorithm */
    void ZCreatePartition(int m, RR sigma, int n, RR tail, RR& v);
    RR ZRecursion(Vec<RR>& X, int m, RR r, RR sigma, RR& v);
    RR NewMarsagliaTailMethod(RR r);
        
    /* Algorithms for Gram-Schmidt orthogonalization */
    double GramSchmidtProcess(Vec< Vec<double> >& T_ATilde, const Vec< Vec<int> >& T_A);
    /* Algorithm for generating the Gram-Schmidt reduced basis for block isometric basis */    
    RR BlockGSO(mat_RR& BTilde, const Vec<ZZX>& B, int n, int precision);
    void FasterIsometricGSO(mat_RR& BTilde, const mat_RR& B);
    
    /* Discrete Gaussian Sampling over lattices */
    ZZX GaussianSamplerFromLattice(const Vec<ZZX>& B, const mat_RR& BTilde, RR sigma, int precision, int tailcut, ZZX center, int n);
    vec_RR CompactGaussianSampler(const Vec< Vec<int> >& B, RR sigma, const vec_RR center, const vec_RR& BTilden, const vec_RR& Vn, const vec_RR& H, const vec_RR& I);    
    void PrepareToSampleCGS(const vec_RR& B1, const mat_RR& BTilde);
    
    /* Procedure for T = BB^t used in Peikert algorithm - Lattice sampling */
    void CholeskyDecomposition(Vec< Vec<double> >& B, const Vec< Vec<double> >& A, int n);    
        
    vec_RR Isometry(const vec_RR& b, int n);
    ZZX Isometry(ZZX& b, int n);    
    void Rot(Vec< Vec<ZZX> >& A, const Vec<ZZX>& a, int m, int n); // RotBasis()
    void rot(Vec<ZZX>& out, const ZZX& b, int n); // Rot()
    void rot(mat_RR& out, const vec_RR& b, int n); // BlockGSO()    
    RR InnerProduct(const vec_RR& a, const vec_RR& b);    
    double InnerProduct(const Vec<int>& a, const Vec<int>& b);
    double InnerProduct(const Vec<int>& a, const Vec<double>& b);
    double InnerProduct(const Vec<double>& a, const Vec<double>& b);    
    double NormOfBasis(const Vec< Vec<double> >& T_ATilde);
    double NormOfBasis(const Vec< Vec<int> >& T_A);
    RR NormOfBasis(const mat_RR& B);
    RR NormOfBasis(const Vec<ZZX>& B, int m, int n);
    RR Norm(const vec_RR& b, int n);
    RR Norm(const ZZX& b, int n);

    void PrintMatrix(const string& name, const Vec< Vec<int> >& matrix);
    void CopyIntToDoubleVec(Vec<double>& B, const Vec<int>& A);
    
};

#endif	/* SAMPLERS_H */

