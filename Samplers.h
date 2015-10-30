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
       
    Samplers();
    Samplers(int q, const ZZ_pX& f);
    virtual ~Samplers();
    
    /* Polynomial generation with coefficient sampled from a discrete Gaussian distribution */
    Vec<int> PolyGeneratorZiggurat(int dimension, int m, RR sigma, int omega, int n, int tail);
    Vec<int> PolyGeneratorKnuthYao(int dimension, int precision, float tailcut, RR sigma, RR c);
    
    void Rot(Vec<ZZX>& A, const Vec<ZZ_pX>& a, int m, int n);

    /* Algorithm for generating the Gram-Schmidt reduced basis */
    RR BlockGSO(mat_RR& BTilde, const Vec<ZZX>& B, int m, int n, int precision);
    
    /* Sampling from the discrete Gaussian distribution D_{\Lambda, \sigma, c}*/
    ZZX GaussianSamplerFromLattice(const Vec<ZZX>& B, const mat_RR& BTilde, RR sigma, int precision, float tailcut, ZZX center, int n);
        
private:
        
    /* Attributes for sampling from lattice */
    Vec< Vec< complex<double> > > V; //Vandermonde matrix
    ZZ_pX f;// Polynomial R = Z_p[X]/f    
        
    /* Knuth-Yao attributes */
    Vec< Vec<int> > P;
    Vec<int> begin;
    
    /* Ziggurat attributes */
    Vec<RR> X;
    Vec<RR> Y;
    Vec<int> X_ZZ;
    Vec<int> Y_ZZ;

    /* Subroutine of Block_GSO algorithm - it produces the reduced form of an isometric basis */
    void FasterIsometricGSO(mat_RR& BTilde, const mat_RR& B);    

    /* Sampling from a discrete Gaussian distribution over the integers */
    int KnuthYao(float tailcut, RR sigma, RR c);
    int Ziggurat(int m, int n, RR sigma, int omega);
    
    /* Auxiliary functions of Ziggurat algorithm */
    void DZCreatePartition(int m, RR sigma, int n, int tail);
    RR DZRecursion(Vec<RR>& X, Vec<RR>& Y, int m, int tail, RR c, RR sigma);
    RR Rho(RR sigma, RR x);
    ZZ sLine(ZZ x0, ZZ x1, ZZ y0, ZZ y1, ZZ x, long int i);
    
    /* Auxiliary functions of Knuth-Yao algorithm */
    void BinaryExpansion(Vec< Vec<int> >& auxP, RR probability, int precision, int index);
    void BuildProbabilityMatrix(int precision, float tailcut, RR sigma, RR c);
    RR Probability(RR x, RR sigma, RR c);

    /* Auxiliary functions of lattice sampler */    
    RR InnerProduct(const vec_RR& a, const vec_RR& b, int n);    
    vec_RR Isometry(vec_RR& b, int n);
    ZZX Isometry(ZZX& b, int n);
    RR Norm(const ZZX& b, int n);
    RR Norm(const vec_RR& b, int n);        
    RR NormOfBasis(const Vec<ZZX>& B, int m, int n);
    RR NormOfBasis(const mat_RR& B, int m, int n);    
    void rot(Vec<ZZX>& out, const ZZ_pX& b, int n);
    void rot(mat_RR& out, const vec_RR& b, int n);        

    void PrintMatrix(const string& name, const Vec< Vec<int> >& matrix);
    
};

#endif	/* SAMPLERS_H */

