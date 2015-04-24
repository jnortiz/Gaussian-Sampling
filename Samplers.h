/* 
 * File:   Samplers.h
 * Author: jnortiz
 *
 * Created on April 24, 2015, 3:51 PM
 */

#include <NTL/ZZ.h>
#include <NTL/RR.h>

#ifndef SAMPLERS_H
#define	SAMPLERS_H

using namespace NTL;

class Samplers {
public:
    
    Samplers();
    Samplers(const Samplers& orig);
    virtual ~Samplers();
    
    Vec<ZZ> PolyGeneratorZiggurat(int dimension, RR m, RR sigma, ZZ omega, RR n);
    
private:
    /* Ziggurat variables */
    Vec<RR> X;
    Vec<RR> Y;
    Vec<ZZ> X_ZZ;
    Vec<ZZ> Y_ZZ;

    /* Sampling from a discrete Gaussian distribution over the integers */
    ZZ ZiggutatO(RR m, RR sigma, ZZ omega);
    /* Auxiliary functions of Ziggurat algorithm */
    void DZCreatePartition(RR m, RR sigma, RR n);
    RR DZRecursion(RR m, RR c, RR sigma);
    RR Rho(RR sigma, RR x);
    ZZ sLine(ZZ x0, ZZ x1, ZZ y0, ZZ y1, ZZ x, long int i);
};

#endif	/* SAMPLERS_H */

