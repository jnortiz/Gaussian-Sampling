/* 
 * File:   HIBE.h
 * Author: jnortiz
 *
 * Created on March 10, 2015, 4:12 PM
 */

#include <NTL/ZZ.h>
#include <NTL/ZZ_pEX.h>

#ifndef HIBE_H
#define	HIBE_H

using namespace std;
using namespace NTL;

class HIBE {
public:
    HIBE(double q, int m1, int m2, int k);
    HIBE(const HIBE& orig);
    virtual ~HIBE();    
    void IdealTrapGen();
    
    int GetM() const {
        return m;
    }

    int GetN() const {
        return n;
    }
    
    long GetM() {
        return m;
    }

    void SetM(long m) {
        this->m = m;
    }

    long GetN() {
        return n;
    }

    void SetN(long n) {
        this->n = n;
    }
    
    ZZ_pX GetF() const {
        return f;
    }

    void SetF(ZZ_pX f) {
        this->f = f;
    }

    long GetM1() const {
        return m1;
    }

    void SetM1(long m1) {
        this->m1 = m1;
    }

    long GetM2() const {
        return m2;
    }

    void SetM2(long m2) {
        this->m2 = m2;
    }
    
    double GetR() const {
        return r;
    }
    
    void SetR(double r) {
        this->r = r;
    }
    
    double GetQ() const {
        return q;
    }

    void SetQ(double q) {
        this->q = q;
    }
    
    double GetLambda() const {
        return lambda;
    }

    void SetLambda(double lambda) {
        this->lambda = lambda;
    }
    
    int GetK() const {
        return k;
    }

    void SetK(int k) {
        this->k = k;
    }
    
    Vec<ZZ_pX> GetA() const {
        return A;
    }

    void SetA(Vec<ZZ_pX> A) {
        this->A = A;
    }
    
    Vec<Vec<ZZX> > GetMsk() const {
        return msk;
    }

    void SetMsk(Vec<Vec<ZZX> > msk) {
        this->msk = msk;
    }

    
private:
    double q;
    double r;
    double lambda;
    long m;
    long n;
    long m1;
    long m2;
    int k;
    ZZ_pX f; // R = Z_p/f and R_0 = Z/f
    
    Vec<ZZ_pX> A; // Part of mpk
    Vec< Vec<ZZX> > msk; // Master secret key of HIBE system
    
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
  
    void PrintMatrixZZX(const string& name, const Vec< Vec<ZZX> >& M);
    void PrintVectorZZX(const string& name, const Vec<ZZX>& M);
    void PrintMatrixZZ_pX(const string& name, const Vec< Vec<ZZ_pX> >& M);
    void PrintVectorZZ_pX(const string& name, const Vec<ZZ_pX>& M);
    
    int FinalVerification(const Vec<ZZ_pX>& A, const Vec< Vec<ZZX> >& S);
    
};

#endif	/* HIBE_H */

