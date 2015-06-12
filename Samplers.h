/* 
 * File:   Samplers.h
 * Author: jnortiz
 *
 * Created on April 24, 2015, 3:51 PM
 */

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/mat_ZZ.h>
#include <complex>
#include <NTL/ZZX.h>

#ifndef SAMPLERS_H
#define	SAMPLERS_H

using namespace std;
using namespace NTL;

class Samplers {
public:
    
    Samplers(int k);
    virtual ~Samplers();
    
    Vec<int> PolyGeneratorZiggurat(int dimension, int m, RR sigma, int omega, int n, int tail);
    Vec<int> PolyGeneratorKnuthYao(int dimension, int precision, int tailcut, RR sigma, RR c);
    void PolyGenerator(ZZX& b, int length, int q);
    ZZX GaussianSamplerFromLattice(const Vec<ZZ_pX>& B, const Vec<ZZX>& BTilde, RR sigma, int precision, int tailcut, ZZX c, int k);
    void FasterIsometricGSO(Vec<ZZX>& BTilde, Vec<ZZ>& C, Vec<double>& D, const Vec<ZZ_pX>& B, int k);
    void DZCreatePartition(int m, RR sigma, int n, int tail);
    RR CoverageAreaZiggurat();
        
private:
        
    /* Attributes for sampling from lattice */
    Vec< Vec< complex<double> > > V; //Vandermonde matrix
    ZZX phi;// Polynomial F = Q[X]/<phi>    
        
    /* Knuth-Yao attributes */
    Vec< Vec<int> > P;
    
    /* Ziggurat attributes */
    Vec<RR> X;
    Vec<RR> Y;
    Vec<int> X_ZZ;
    Vec<int> Y_ZZ;

    /* Sampling from a discrete Gaussian distribution over the integers */
    int Ziggurat(int m, RR sigma, int omega);
    int KnuthYao(int precision, int tailcut, RR sigma);
    
    /* Auxiliary functions of Ziggurat algorithm */
//    void DZCreatePartition(int m, RR sigma, int n, int tail);
    RR DZRecursion(Vec<RR>& X, Vec<RR>& Y, int m, int tail, RR c, RR sigma);
    RR Rho(RR sigma, RR x);
    ZZ sLine(ZZ x0, ZZ x1, ZZ y0, ZZ y1, ZZ x, long int i);
    
    /* Auxiliary functions of Knuth-Yao algorithm */
    RR Probability(RR x, RR sigma, RR c);
    void BinaryExpansion(Vec< Vec<int> >& auxP, RR probability, int precision, int index);
    void BuildProbabilityMatrix(int precision, int tailcut, RR sigma, RR c);

    /* Auxiliary functions of lattice sampler */
    void BuildVandermondeMatrix(int k);
    int EulerPhiPowerOfTwo(int k);
    void ConjugateOfMatrix(Vec< Vec< complex<double> > >& M);    
    void ComplexMatrixMult(Vec<ZZX>& c, const Vec< Vec< complex<double> > >& a, const Vec< Vec< complex<double> > >& b);
    ZZ InnerProduct(const ZZX& a, const ZZX& b);
    double Norm(const ZZX& b);
    ZZX Isometry(ZZX& b);
    
    void PrintMatrix(const string& name, const Vec< Vec<int> >& matrix);
    ZZX Mult(ZZX V, double c, int m);
//    RR CoverageAreaZiggurat();
    
};

#endif	/* SAMPLERS_H */

