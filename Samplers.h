/* 
 * File:   Samplers.h
 * Author: jnortiz
 *
 * Created on April 24, 2015, 3:51 PM
 */

#include <NTL/RR.h>

#ifndef SAMPLERS_H
#define	SAMPLERS_H

using namespace std;
using namespace NTL;

class Samplers {
public:
       
    Samplers();
    virtual ~Samplers();
    
    /* Polynomial generation with coefficients sampled from a discrete Gaussian distribution */
    Vec<int> PolyGeneratorZiggurat(int dimension, int m, RR sigma, int omega, int n, int tail);
    Vec<int> PolyGeneratorKnuthYao(int dimension, int precision, float tailcut, RR sigma, RR c);
    
private:
        
    /* Knuth-Yao attributes */
    Vec< Vec<int> > P;
    Vec<int> begin;
    
    /* Ziggurat attributes */
    Vec<RR> X;
    Vec<RR> Y;
    Vec<int> X_ZZ;
    Vec<int> Y_ZZ;

    /* Sampling from a discrete Gaussian distribution over the integers */
    int KnuthYao(float tailcut, RR sigma, RR c);
    int Ziggurat(int m, int n, RR sigma, int omega);
    
    /* Auxiliary functions of Ziggurat algorithm */
    void DZCreatePartition(int m, RR sigma, int n, int tail);
    RR DZRecursion(Vec<RR>& X, Vec<RR>& Y, int m, int tail, RR c, RR sigma);
    RR Rho(RR sigma, RR x);
    
    /* Auxiliary functions of Knuth-Yao algorithm */
    void BinaryExpansion(Vec< Vec<int> >& auxP, RR probability, int precision, int index);
    void BuildProbabilityMatrix(int precision, float tailcut, RR sigma, RR c);
    RR Probability(RR x, RR sigma, RR c);

    void PrintMatrix(const string& name, const Vec< Vec<int> >& matrix);
    
};

#endif	/* SAMPLERS_H */

