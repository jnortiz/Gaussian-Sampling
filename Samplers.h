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
    
    /* Generic algorithm for GSO basis construction */
    void GramSchmidtProcess(Vec< Vec<double> >& T_ATilde, const Vec< Vec<int> >& T_A, int n);
    /* Sampling from the discrete Gaussian distribution D_{\Lambda, \sigma, c}*/
    ZZX GaussianSamplerFromLattice(const Vec<ZZX>& B, const mat_RR& BTilde, RR sigma, int precision, int tailcut, ZZX center, int n);

    double NormOfBasis(const Vec< Vec<double> >& T_ATilde);
    double NormOfBasis(const Vec< Vec<int> >& T_A);
    
private:
        
    /* Attributes for sampling from lattice */
    ZZ_pX f;// Polynomial R = Z_p[X]/f    
    
    /* Knuth-Yao attributes */
    Vec< Vec<int> > P;
    Vec<int> begin;
    
    /* Sampling from a discrete Gaussian distribution over the integers */
    int KnuthYao(int tailcut, RR sigma, RR c);
    
    /* Auxiliary functions of Knuth-Yao algorithm */
    void BuildProbabilityMatrix(int precision, int tailcut, RR sigma, RR c);
    void BinaryExpansion(Vec< Vec<int> >& auxP, RR probability, int precision, int index);
    RR Probability(RR x, RR sigma, RR c);

    /* Auxiliary functions of sampling from lattices */    
    RR InnerProduct(const vec_RR& a, const vec_RR& b, int n);    
    double InnerProduct(const Vec<int>& a, const Vec<int>& b);
    double InnerProduct(const Vec<int>& a, const Vec<double>& b);
    double InnerProduct(const Vec<double>& a, const Vec<double>& b);
    
    RR Norm(const vec_RR& b, int n);

    void PrintMatrix(const string& name, const Vec< Vec<int> >& matrix);
    void CopyIntToDoubleVec(Vec<double>& B, const Vec<int>& A);
    
};

#endif	/* SAMPLERS_H */

