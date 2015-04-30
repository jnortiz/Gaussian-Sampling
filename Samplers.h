/* 
 * File:   Samplers.h
 * Author: jnortiz
 *
 * Created on April 24, 2015, 3:51 PM
 */

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/mat_ZZ.h>

#ifndef SAMPLERS_H
#define	SAMPLERS_H

using namespace std;
using namespace NTL;

class Samplers {
public:
    
    Samplers();
    Samplers(const Samplers& orig);
    virtual ~Samplers();
    
    Vec<int> PolyGeneratorZiggurat(int dimension, RR m, RR sigma, ZZ omega, RR n, RR tail);
    Vec<int> PolyGeneratorKnuthYao(int dimension, int precision, int tailcut, RR sigma);
    
private:
    /* Knuth-Yao attributes */
    mat_ZZ P; // Probability matrix such as each polynomial will be a binary expansion of a probability
    
    /* Ziggurat attributes */
    Vec<RR> X;
    Vec<RR> Y;
    Vec<ZZ> X_ZZ;
    Vec<ZZ> Y_ZZ;

    /* Sampling from a discrete Gaussian distribution over the integers */
    int ZiggutatO(RR m, RR sigma, ZZ omega);
    int KnuthYao(int precision, int tailcut, RR sigma);
    
    /* Auxiliary functions of Ziggurat algorithm */
    void DZCreatePartition(RR m, RR sigma, RR n, RR tail);
    RR DZRecursion(RR m, RR c, RR sigma);
    RR Rho(RR sigma, RR x);
    ZZ sLine(ZZ x0, ZZ x1, ZZ y0, ZZ y1, ZZ x, long int i);
    
    /* Auxiliary functions of Knuth-Yao algorithm */
    RR Probability(RR x, RR sigma);
    void BinaryExpansion(mat_ZZ& aux_P, RR probability, int precision, int index);
    void BuildProbabilityMatrix(int precision, int tailcut, RR sigma);

    void PrintMatrix(const string& name, const mat_ZZ& matrix);
    
};

#endif	/* SAMPLERS_H */

