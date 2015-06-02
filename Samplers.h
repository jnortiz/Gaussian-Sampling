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
    Vec<int> PolyGeneratorKnuthYao(int dimension, int precision, int tailcut, RR sigma);
    void InnerProduct(ZZ& out, const ZZX& a, const ZZX& b);
    void Norm(double& out, const ZZX& b);
    void Isometry(ZZX& out, ZZX& b);
    void PolyGenerator(ZZX& b, int length, int q);
        
private:
    
    ZZX phi;
    
    /* Attributes for sampling from lattice */
    Vec< Vec< complex<double> > > V; //Vandermonde matrix
        
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
    void DZCreatePartition(int m, RR sigma, int n, int tail);
    RR DZRecursion(int m, int c, RR sigma);
    RR Rho(RR sigma, RR x);
    ZZ sLine(ZZ x0, ZZ x1, ZZ y0, ZZ y1, ZZ x, long int i);
    
    /* Auxiliary functions of Knuth-Yao algorithm */
    RR Probability(RR x, RR sigma);
    void BinaryExpansion(RR probability, int precision, int index);
    void BuildProbabilityMatrix(int precision, int tailcut, RR sigma);

    /* Auxiliary functions of lattice sampler */
    void BuildVandermondeMatrix(int k);
    int EulerPhiPowerOfTwo(int k);
    void ConjugateOfMatrix(Vec< Vec< complex<double> > >& M);    
    void ComplexMatrixMult(Vec<ZZX>& c, const Vec< Vec< complex<double> > >& a, const Vec< Vec< complex<double> > >& b);
    
    void PrintMatrix(const string& name, const Vec< Vec<int> >& matrix);
    
};

#endif	/* SAMPLERS_H */

