/* 
 * File:   HIBE.h
 * Author: jnortiz
 *
 * Created on March 10, 2015, 4:12 PM
 */

#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/RR.h>
#include "Samplers.h"

#ifndef HIBE_H
#define	HIBE_H

using namespace std;
using namespace NTL;

class HIBE {
public:
    
    void Setup(int h);
    HIBE(double q, int m1, int m2, int k);
    virtual ~HIBE();    
    
    long GetM() { return m; }    
    long GetN() { return n; }
    ZZ_pX GetF() const { return f; }
    long GetM1() const { return m1; }
    long GetM2() const { return m2; }    
    double GetR() const { return r; }    
    double GetQ() const { return q; }    
    double GetLambda() const { return lambda; }    
    int GetK() const { return k; }    
    Vec<ZZ_pX> GetA() const { return A; }    
    Vec<Vec<ZZ_pX> > GetA_prime() const { return A_prime; }
    Vec<ZZ_pX> GetB() const { return B; }
    ZZ_pX GetU() const { return u; }
    Vec<Vec<ZZX> > GetMsk() const { return msk; }
    Samplers* GetSampler() const { return sampler; }
    
private:
    /* Global parameters */
    double q;
    double r;
    double lambda;
    long m;
    long n;
    long m1;
    long m2;
    int k;
    ZZ_pX f; // R = Z_p/f and R_0 = Z/f
    
    /* Master public key */
    Vec<ZZ_pX> A;
    Vec< Vec<ZZ_pX> > A_prime;
    Vec<ZZ_pX> B;
    ZZ_pX u;
    
    /* Master secret key */
    Vec< Vec<ZZX> > msk;
    
    /* Sampling from Gaussian distributions */
    Samplers *sampler;
    
    /* Methods */
    int IdealTrapGen();   
    
    /* Auxiliary functions of IdealTrapGen algorithm */
    void Decomp(Vec< Vec<ZZX> >& _W, const Vec<ZZX>& _h);
    void Mult(ZZ_pX& c, const Vec<ZZX>& a, const Vec<ZZ_pX>& b, const ZZ_pX& f);
    void Mult(Vec<ZZX>& c, const ZZ_pX& a, const Vec<ZZX>& b);    
    void Mult(Vec<ZZX>& c, const ZZX& a, const Vec<ZZX>& b);
    void Mult(Vec< Vec<ZZX> >& c, const Vec< Vec<ZZX> >& a, const Vec< Vec<ZZX> >& b);
    void Mult(ZZX& c, const Vec<ZZX>& a, const Vec<ZZX>& b);     
    void Mult(Vec<ZZ_pX>& c, const Vec< Vec<ZZX> >& a, const Vec<ZZ_pX>& b);
    void Mult(ZZ_pX& c, const Vec<ZZX>& a, const Vec<ZZ_pX>& b);     
    void Add(Vec<ZZX>& c, const Vec<ZZX>& a, const Vec<ZZX>& b);
    void Concat(Vec< Vec<ZZX> >& S, const Vec< Vec<ZZX> >& V, const Vec< Vec<ZZX> >& P, const Vec< Vec<ZZX> >& D, const Vec< Vec<ZZX> >&B);
    void Concat(Vec<ZZ_pX>& A, const Vec<ZZ_pX>& A1, const Vec<ZZ_pX>& A2);
    int FinalVerification(const Vec<ZZ_pX>& A, const Vec< Vec<ZZX> >& S);
  
    /* Auxiliary functions */
    void PrintMatrixInt(const string& name, const Vec< Vec<int> >& M);
    void PrintVectorZZX(const string& name, const Vec<ZZX>& M);
    void PrintMatrixZZ_pX(const string& name, const Vec< Vec<ZZ_pX> >& M);
    void PrintVectorZZ_pX(const string& name, const Vec<ZZ_pX>& M);
    
};

#endif	/* HIBE_H */

