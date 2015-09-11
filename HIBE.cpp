/* 
 * File:   HIBE.cpp
 * Author: jnortiz
 * 
 * Created on March 10, 2015, 4:12 PM
 * 
 */

#include "HIBE.h"
#include <math.h>
#include <cstdio>
#include <NTL/ZZ_p.h>
#include <NTL/ZZXFactoring.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/RR.h>
#include <cstring>

using namespace std;
using namespace NTL;

HIBE::HIBE(double q, int m1, int m2, int k) {
        
    this->lambda = ceil(1 + (log(q)/log(2)));    
    this->m = m1 + m2;
    this->n = pow(2, k);
    this->q = q;
    this->m1 = m1;
    this->m2 = m2;
    this->k = k;    
    this->r = ceil((double)(1 + (log(q)/log(3))));
    
    ZZ_p::init(conv<ZZ>(q)); // Coefficients modulo

    ZZ_pX f;
    f.SetLength((this->n)+1);
    SetCoeff(f, 0, 1);
    SetCoeff(f, this->n, 1);
    this->f = f;
    ZZ_pE::init(f); // Ring elements modulo
    this->sampler = new Samplers(q, f);
}

HIBE::~HIBE() {
    delete(this->sampler);
}

void HIBE::Setup(int h) {
    
    /* The IdealTrapGen algorithm selects uniformly a random vector A and a short basis S = msk */
    cout << "[*] IdealTrapGen status: ";
    if(this->IdealTrapGen())
        cout << "Pass!" << endl;
    else
        cout << "Error!" << endl;
    
    /* The A_prime matrix contains uniformly random vectors, one for each hierarchy level */
    int i, j;
    this->A_prime.SetLength(h);
    for(i = 0; i < A_prime.length(); i++) {
        this->A_prime[i].SetLength(this->m);
        for(j = 0; j < A_prime[0].length(); j++)
            this->A_prime[i][j] = random_ZZ_pX(this->n);
    }//end-for
    
    /* B is just a uniformly random vector of m elements */
    this->B.SetLength(this->m);
    for(i = 0; i < B.length(); i++)
        this->B[i] = random_ZZ_pX(this->n);
    
    /* u is a uniformly random polynomial in R = Z_p[x]/f */
    this->u = random_ZZ_pX(this->n);
    
    /* Uncomment the following lines if you want a preview of master keys */
//    cout << "/** Master public key (A, A_prime, B, u) **/" << endl;
//    this->PrintVectorZZ_pX("Vector A", this->A);
//    this->PrintMatrixZZ_pX("Matrix A_prime", this->A_prime);
//    this->PrintVectorZZ_pX("Vector B", this->B);
//    cout << "\n/** Polynomial u **/\n" << this->u << endl;    
    
    cout << "[*] Setup status: Pass!" << endl;
    
}//end-Setup()

int HIBE::IdealTrapGen() {
    
    Vec< Vec<ZZX> > T; // Matrix T_{lambda} -- it has elements in {-2, 1, 0}
    Vec< Vec<ZZX> > S;  // It will be the master secret key of HIBE system
    
    Vec<ZZ_pX> A; // A vector part of the master public key
    Vec<ZZ_pX> A1; // Vector of random polynomials in R
    Vec<ZZX> R; // Random vector with polynomials in {-1, 0, 1}
        
    /* Auxiliary vectors */
    ZZ_pX pX;
    ZZX p;
    
    int i, j;
    
    p.SetLength(this->n);
    pX.SetLength(this->n);
    
    /* Random vector R */
    R.SetLength(this->m1);    
    for(j = 0; j < this->r; j++) {
        for(i = 0; i < this->n; i++)
            p[i] = rand()%3 - 1;
        R[j] = p;
    }//end-for        
    
    for(j = this->r; j < this->m1; j++) {
        for(i = 0; i < this->n; i++)
            p[i] = 0;
        R[j] = p;
    }//end-for        
    
    /* A1 - vector of random polynomials */
    A1.SetLength(this->m1);
    for(i = 0; i < this->m1; i++) {
        random(pX, this->n);        
        A1[i] = pX;
    }//end-for    

    /* Matrix T_{lambda} */
    T.SetLength(this->lambda);    
    for(i = 0; i < this->lambda; i++)
        T[i].SetLength(this->lambda);
    
    for(i = 0; i < this->lambda; i++) {
        for(j = 0; j < this->lambda; j++) {            
            if(i == j)
                T[i][j] = 1;
            else if(i == j+1)
                    T[i][j] = -2;
                else
                    T[i][j] = 0;
        }//end-for
    }//end-for               
    
    /* A huge B matrix */    
    Vec< Vec<ZZX> > B;
    
    B.SetLength(this->GetM2());
    for(i = 0; i < this->GetM2(); i++)
        B[i].SetLength(this->GetM2());
    
    for(i = 0; i < this->GetM2(); i++) {      
        for(j = 0; j < this->GetM2(); j++) {
            B[i][j] = 0;
            if(i == j)
                B[i][j] = 1;
            if(i == j+1 && i < this->lambda*this->m1)
                B[i][j] = -2;
        }//end-for
    }//end-for
    
    /* H's construction -- it's required vectors e[i] and matrix A1 */    
    Vec< Vec<ZZX> > H;    
    Vec< Vec<ZZX> > e;    
    Vec<ZZ_pX> A1_copy;
    ZZ_pX hi, inv_a1;
    ZZX aux_swap;
    
    A1_copy.SetLength(this->m1);
    inv_a1.SetLength(this->n);
    hi.SetLength(this->n);    
    
    H.SetLength(this->m1);    
    for(i = 0; i < this->m1; i++) {
        H[i].SetLength(this->m1);    
        for(j = 0; j < this->m1; j++)
            H[i][j].SetLength(this->n);
    }//end-for
        
    /* Construction of all e_i at once */
    e.SetLength(this->m1);
    for(i = 0; i < this->m1; i++) {
        e[i].SetLength(this->m1);
        for(j = 0; j < this->m1; j++) {
            e[i][j].SetLength(this->n);
            e[i][j] = conv<ZZX>(0);
            if(i == j)
                e[i][j] = conv<ZZX>(1);
        }//end-for
    }//end-for    
    
    /* From now, we will manipulate A1_copy matrix in H construction */
    for(i = 0; i < A1.length(); i++)
        A1_copy[i] = A1[i];
             
    int iStar;
    for(iStar = 0; iStar < this->m1; iStar++) {
        /* If A1_aux[iStar] is invertible, InvModStatus returns 0 and calculates its inverse */
        if(InvModStatus(inv_a1, A1_copy[iStar], this->f) == 0) {
            
            if(iStar != 0) {
                /* When i_Star is not 1 (paper notation), swap rows 1 and iStar */
                A1_copy.put(0, A1.get(iStar));
                A1_copy.put(iStar, A1.get(0));
                
            }//end-if
            
            /* First row is qe1 (paper notation) -- that's is actually e[0] */
            H[0] = e[0];
            H[0][0] = conv<ZZX>(this->q);
            
            /* Computation of i-th row, with i = {2, ..., m1} */
            for(i = 1; i < this->m1; i++) {                
                MulMod(hi, -A1_copy[i], inv_a1, this->f);
                this->Mult(H[i], conv<ZZX>(hi), e[0]);
                this->Add(H[i], H[i], e[i]);
            }//end-for
            
            break;
            
        }//end-if        
    }//end-for
    
    /* Post-processing. If iStar is different of 1 (paper notation), we must swap columns 1 and iStar */    
    if(iStar != 0) {
        for(i = 0; i < H[0].length(); i++) {
            aux_swap = H[i][iStar];
            H[i][iStar] = H[i][0];
            H[i][0] = aux_swap;
        }//end-for
        
    }//end-if    
        
    Vec< Vec<ZZX> > H_prime;
    Vec< Vec<ZZX> > identity;
    
    /* H' = H - I_{m1} */
    H_prime.SetLength(this->m1);    
    identity.SetLength(this->m1);    
    for(i = 0; i < this->m1; i++) {
        identity[i].SetLength(this->m1);    
        H_prime[i].SetLength(this->m1);    
        for(j = 0; j < this->m1; j++) {
            identity[i][j].SetLength(this->n);
            H_prime[i][j].SetLength(this->n);
            if(i == j)
                identity[i][j] = conv<ZZX>(1);
            else 
                identity[i][j] = conv<ZZX>(0);
            sub(H_prime[i][j], H[i][j], identity[i][j]);
        }//end-for
    }//end-for
    
    /* Tk matrix inversion */
    mat_ZZ convTk, invTk;
    Vec< Vec<ZZX> > invT;
    
    convTk.SetDims(this->lambda, this->lambda);
    invTk.SetDims(this->lambda, this->lambda);
    invT.SetLength(this->lambda);
    
    for(i = 0; i < convTk.NumRows(); i++) {
        invT[i].SetLength(this->lambda);
        for(j = 0; j < convTk.NumCols(); j++)
            convTk[i][j] = ConstTerm(T[i][j]);
    }//end-for
    
    inv(invTk, convTk);
    
    /* Conversion of invTk from mat_ZZ to Vec< Vec<ZZX> >*/
    for(i = 0; i < this->lambda; i++)
        for(j = 0; j < this->lambda; j++)
            invT[i][j] = conv<ZZX>(invTk[i][j]); //Each ZZ becomes a ZZX
        
    Vec< Vec<ZZX> > W;    
    W.SetLength(this->m2); //Number of columns    
    for(i = 0; i < this->m2; i++) {
        W[i].SetLength(this->m1); //Number of lines
        for(j = 0; j < this->m1; j++)
            W[i][j].SetLength(this->n);
    }//end-for
    
    /* W's construction, such as G = B^{-1}*W */    
    Vec< Vec<ZZX> > inDecomp;
    Vec< Vec<ZZX> > Gj;
    Vec < Vec<ZZX> > G;
    int wIndex = 0;
    
    G.SetLength(this->m2);
    for(i = 0; i < this->m2; i++) {
        G[i].SetLength(this->m1); //Number of lines
        for(j = 0; j < this->m1; j++)
            G[i][j].SetLength(this->n);
    }//end-for    
    
    // Construction of W and G matrices    
    for(i = 0; i < this->m1; i++) { //For each row h'_j of H'
        // W's construction
        this->Decomp(inDecomp, H_prime[i]); //Run the binary decomposition algorithm for h'_j
        // G's construction -- Tk^{-1}*W_j, such as W_j = inDecomp
        this->Mult(Gj, invT, inDecomp); //Sub-matrix multiplication
        
        // Copying intermediate results to W and G
        for(j = 0; j < this->lambda; j++) { //Copy the block of (lambda, m1) polynomials to W
            for(int k = 0; k < this->m1; k++) {
                W[wIndex][k] = inDecomp[j][k];
                G[wIndex][k] = Gj[j][k];
            }
            wIndex++;
        }//end-for
        
    }//end-for
    
    // The remaining is filled with zero polynomials
    for(i = this->m1*this->lambda; i < this->m2; i++)
        for(j = 0; j < this->m1; j++)
            for(int k = 0; k < this->n; k++) {
                W[i][j][k] = 0; 
                G[i][j][k] = 0;
            }//end-for
            
    /* U's construction: U = G + R */
    Vec< Vec<ZZX> > U;
    Vec< Vec<ZZX> > minusU;
    
    minusU.SetLength(this->m2);
    U.SetLength(this->m2);
    for(i = 0; i < U.length(); i++) {
        U[i].SetLength(this->m1);
        minusU[i].SetLength(this->m1);
    }//end-for
    
    for(i = 0; i < this->m2; i++)
        Add(U[i], G[i], R);
    
    for(i = 0; i < minusU.length(); i++)
        for(j = 0; j < minusU[0].length(); j++)
            minusU[i][j] = -U[i][j];
    
    Vec<ZZ_pX> A2;
    A2.SetLength(this->m2);
    
    /* A2 = -U*A1 */
    this->Mult(A2, minusU, A1);
    
    /* P's construction: p_j = e_{lambda*j} \in R_0^{m2} */   
    Vec < Vec<ZZX> > P;
    
    P.SetLength(this->m1);    
    for(i = 0; i < this->m1; i++) {
        P[i].SetLength(this->m2);
        for(j = 0; j < this->m2; j++) {
            P[i][j].SetLength(this->n);
            P[i][j] = conv<ZZX>(0);
        }//end-for
        P[i][((i+1)*this->lambda) - 1] = conv<ZZX>(1);
    }//end-for
    
    Vec< Vec<ZZX> > D, V, Y, minusH;    
    minusH.SetLength(this->m1);
    V.SetLength(this->m1);    
    
    D.SetLength(this->m2);
    for(i = 0; i < D.length(); i++)
        D[i].SetLength(this->m1);
    
    // Y = P*U -- just an intermediate step...
    this->Mult(Y, P, U);
    
    for(i = 0; i < minusH.length(); i++) {
        minusH[i].SetLength(this->m1);
        V[i].SetLength(this->m1);
        for(int j = 0; j < minusH[0].length(); j++)
            minusH[i][j] = -H[i][j];
    }//end-for    

    /* V = -H + P*U */
    for(i = 0; i < minusH.length(); i++)
        this->Add(V[i], minusH[i], Y[i]);
    
    /* D = B*U */
    this->Mult(D, B, U);
    
    this->Concat(S, V, P, D, B);    
    this->Concat(A, A1, A2);
    
    if(this->FinalVerification(A, S)) {
        this->A = A;
        this->msk = S;
//        cout << "\n\n/** Short basis S **/\n" << S << endl;
//        cout << "\n\n/** Random vector A **/\n" << A << endl;
        return 1;
    } else
        return 0;
    
}//end-IdealTrapGen()

void HIBE::Decomp(Vec< Vec<ZZX> >& _W, const Vec<ZZX>& _h) {
    
    ZZ coeff;   
    
    _W.SetLength(this->m1);
    
    for(int i = 0; i < this->m1; i++) {
        _W[i].SetLength(this->lambda); 
        for(int j = 0; j < this->m1; j++)
            _W[i][j].SetLength(this->n);        
    }//end-for

    for(int i = 0; i < this->m1; i++) { //For each polynomial in h'_j
        for(int j = 0; j < this->n; j++) { //For each coefficient in each polynomial
            coeff = _h[i][j];
            for(int k = this->lambda-1; k >= 0; k--) { //For each power of 2
                if(coeff % 2 == 0)
                    _W[k][i][j] = 0;
                else 
                    _W[k][i][j] = 1;
                coeff = coeff / 2;
            }//end-for
        }//end-for
    }//end-for
    
}//end-Decomp()

void HIBE::Mult(ZZ_pX& c, const Vec<ZZX>& a, const Vec<ZZ_pX>& b, const ZZ_pX& f) {
    
    ZZ_pX acc;
    c.SetLength(this->n);
    acc.SetLength(this->n);
    c.zero();
    
    for(int i = 0; i < a.length(); i++) {
        MulMod(acc, conv<ZZ_pX>(a[i]), b[i], f);
        c += acc;
    }//end-for
                              
}

void HIBE::Mult(Vec<ZZX>& c, const ZZ_pX& a, const Vec<ZZX>& b) {
    
    for(int i = 0; i < b.length(); i++)
        mul(c[i], conv<ZZX>(a), b[i]);

}//end-Mult()

void HIBE::Mult(Vec<ZZX>& c, const ZZX& a, const Vec<ZZX>& b) {

    for(int i = 0; i < b.length(); i++)
        c[i] = a*b[i];    
    
}//end-Mult()

void HIBE::Mult(Vec< Vec<ZZX> >& c, const Vec< Vec<ZZX> >& a, const Vec< Vec<ZZX> >& b) {       
    
    Vec<ZZX> column;
    c.SetLength(a.length());
    column.SetLength(b.length());    

    for(int i = 0; i < c.length(); i++)
        c[i].SetLength(b[0].length());
    
    for(int i = 0; i < a.length(); i++) {
        for(int j = 0; j < b[0].length(); j++) { 
            for(int h = 0; h < b.length(); h++)
                column[h] = b[h][j];
            this->Mult(c[i][j], a[i], column);//[*]  
        }//end-for
    }//end-for
    
}//end-Mult()

void HIBE::Mult(ZZX& c, const Vec<ZZX>& a, const Vec<ZZX>& b) {
    
    ZZX out;
    c.SetLength(this->n);
    out.SetLength(this->n);
    c = conv<ZZX>(0);    
        
    for(int i = 0; i < a.length(); i++) {
        mul(out, a[i], b[i]);
        add(c, c, out);
    }//end-for
                              
}//end-Mult()

void HIBE::Mult(Vec<ZZ_pX>& c /*A2*/, const Vec< Vec<ZZX> >& a /*-U*/, const Vec<ZZ_pX>& b/*A1*/) {
    
    c.SetLength(a.length());

    for(int i = 0; i < a.length(); i++)
        this->Mult(c[i], a[i], b);//[1]
    
}//end-Mult()

void HIBE::Mult(ZZ_pX& c, const Vec<ZZX>& a, const Vec<ZZ_pX>& b) {
    
    ZZ_pX out;
    c.SetLength(this->n);
    out.SetLength(this->n);
    c = conv<ZZ_pX>(0);    
        
    for(int i = 0; i < a.length(); i++) {
        MulMod(out, conv<ZZ_pX>(a[i]), b[i], this->f);
        add(c, c, out);
    }//end-for
                              
}//end-Mult()

void HIBE::Add(Vec<ZZX>& c, const Vec<ZZX>& a, const Vec<ZZX>& b) {
    
    for(int i = 0; i < c.length(); i++)
        c[i] = a[i] + b[i];
    
}//end-Add()

void HIBE::Concat(Vec< Vec<ZZX> >& S, const Vec< Vec<ZZX> >& V, const Vec< Vec<ZZX> >& P, const Vec< Vec<ZZX> >& D, const Vec< Vec<ZZX> >&B) {
    
    int i, j, k;
    int iS, jS;
    
    S.SetLength(this->m);    
    for(i = 0; i < S.length(); i++) {
        S[i].SetLength(this->m);
        for(k = 0; k < this->m; k++) 
            S[i][k].SetLength(this->n);
    }//end-for
    
    /* Copying V matrix to top-left of S */
    for(i = 0; i < V.length(); i++)
        for(j = 0; j < V[0].length(); j++)
            S[i][j] = V[i][j];
    
    /* Copying P matrix to top-right of S */
    for(i = 0; i < P.length(); i++)
        for(j = 0, jS = V[0].length(); j < P[0].length() && jS < S[0].length(); j++, jS++)
            S[i][jS] = P[i][j];
    
    /* Copying D matrix to down-right of S */    
    for(i = 0, iS = V.length(); i < D.length() && iS < S.length(); i++, iS++)
        for(j = 0; j < D[0].length(); j++)
            S[iS][j] = D[i][j];
    
    /* Copying H matrix to down-left of S */    
    for(i = 0, iS = P.length(); i < B.length() && iS < S.length(); i++, iS++)
        for(j = 0, jS = D[0].length(); j < B[0].length(), jS < S[0].length(); j++, jS++)
            S[iS][jS] = B[i][j];
    
}//end-Concat()

void HIBE::Concat(Vec<ZZ_pX>& A, const Vec<ZZ_pX>& A1, const Vec<ZZ_pX>& A2) {
    
    A.SetLength(this->m);
    
    int i, iA;
    
    for(i = 0; i < A1.length(); i++)
        A[i] = A1[i];
    
    for(i = 0, iA = A1.length(); i < A2.length() && iA < A.length(); i++, iA++)
        A[iA] = A2[i];
    
}//end-Concat()

void HIBE::PrintMatrixInt(const string& name, const Vec< Vec<int> >& M) {
    
    cout << "\n/** " << name << " **/" << endl;
    for(int i = 0; i < M.length(); i++) {
        cout << M[i] << endl;
    }//end-for
    
}//end-PrintMatrixZZX()

void HIBE::PrintVectorZZX(const string& name, const Vec<ZZX>& V) {
    
    cout << "\n/** " << name << " **/" << endl;
    for(int i = 0; i < V.length(); i++)
        cout << V[i] << " ";
    cout << endl;
    
}//end-PrintVectorZZX()

void HIBE::PrintMatrixZZ_pX(const string& name, const Vec< Vec<ZZ_pX> >& M) {
    
    cout << "\n/** " << name << " **/" << endl;
    for(int i = 0; i < M.length(); i++) {
        for(int j = 0; j < M[0].length(); j++)
            cout << M[i][j] << " ";
        cout << endl;
    }//end-for
    
}//end-PrintMatrixZZ_pX()

void HIBE::PrintVectorZZ_pX(const string& name, const Vec<ZZ_pX>& V) {
    
    cout << "\n/** " << name << " **/" << endl;
    for(int i = 0; i < V.length(); i++)
        cout << V[i] << " ";
    cout << endl;
    
}//end-PrintVectorZZ_pX()

int HIBE::FinalVerification(const Vec<ZZ_pX>& A, const Vec< Vec<ZZX> >& S) {

    Vec<ZZ_pX> result;
    ZZ_pX acc;    
    result.SetLength(A.length());
    acc = conv<ZZ_pX>(0);
   
    for(int i = 0; i < S.length(); i++) {
        this->Mult(result[i], S[i], A);
        acc += result[i];
    }//end-for
    
    if(acc == acc.zero())
        return 1;
    
    return 0;
    
}//end-FinalVerification()