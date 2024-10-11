// Description:
// Computes the grid algorithm


//#include <iostream>
using namespace std;

#include <R.h>
#include "Rmath.h"
#include <R_ext/Rdynload.h>

extern "C" {
  
  SEXP simplegrid(SEXP Rr, SEXP RU, SEXP RY, SEXP RNg, SEXP RNc, SEXP RD, SEXP rho) {

    R_len_t p = Rf_length(Rr), n = Rf_nrows(RU);

    if(!Rf_isFunction(RD) & (RD != R_NilValue)) perror("RD must be a function");
    if(!Rf_isEnvironment(rho)) perror("rho must be an environment");

    SEXP R_fcall, ans, RUvect;

    int *Nc, *Ng;
    double *r, *U, *Y;

    PROTECT(RNc = AS_INTEGER(RNc));
    PROTECT(RNg = AS_INTEGER(RNg));

    PROTECT(Rr = Rf_coerceVector(Rr, REALSXP));
    PROTECT(RU = Rf_coerceVector(RU, REALSXP));
    PROTECT(RY = Rf_coerceVector(RY, REALSXP));

    Nc = INTEGER(RNc);
    Ng = INTEGER(RNg);

    r = REAL(Rr);
    U = REAL(RU);
    Y = REAL(RY);

    void deheuvels(double *r, double *s, int *ranking, int *n, double *resBn);

    int i, k, j, l, m;
    double *Dvect, *gridtheta, *vectr, *resBn, *Uvect;
    int *ranking;
    gridtheta = new double [Ng[0]];
    Dvect = new double [Ng[0]];
    double theta, theta0;
    vectr = new double[p];
    double norm;
    ranking = new int[1];
    ranking[0] = 1;
    resBn = new double[1];
    int indmaxDvect;
    double maxDvect;
    Uvect = new double [n];
    int *nptr;
    nptr = new int[1];
    nptr[0] = n;

    PROTECT(RUvect = Rf_allocVector(REALSXP, n));
    Uvect = REAL(RUvect);

    for (i=1;i<=Nc[0];i++) {
      //      for (k=0;k<Ng[0];k++) gridtheta[k] = - M_PI/(pow(2.0,(double)(i-1))) * (0.5 - ((double)k)/((double)(Ng[0]-1))); 
      for (k=0;k<Ng[0];k++) gridtheta[k] = - (M_PI/(pow(2.0,(double)(i-2)))) * (0.5 - ((double)k)/((double)(Ng[0]))); 
      for (j=1;j<=p;j++) {
	for (l=1;l<=Ng[0];l++) {
	  theta = gridtheta[l-1];
	  for (k=0;k<p;k++) {if ((k+1)==j) vectr[k] = cos(theta) * r[k] + sin(theta); else vectr[k] = cos(theta) * r[k];}
	  norm = 0.0;
	  for (k=0;k<p;k++) norm = norm + vectr[k]*vectr[k];
	  for (k=0;k<p;k++) vectr[k] = vectr[k]/sqrt(norm);	    
	  for (k=0;k<n;k++) {
	    Uvect[k] = 0.0;
	    for (m=0;m<p;m++) Uvect[k] = Uvect[k] + U[m*n+k] * vectr[m];
	  }
	  if (RD == R_NilValue) {
	    deheuvels(Uvect, Y, ranking, nptr, resBn);
	    Dvect[l-1] = resBn[0];
	  } else {
	    PROTECT( R_fcall = Rf_lang3(RD, RUvect, RY) );
	    ans = Rf_eval(R_fcall, rho);
	    UNPROTECT(1);
	    Dvect[l-1]= REAL(ans)[0];
	  }
	}
	indmaxDvect = 0;
	maxDvect = Dvect[0];
	for (k=1;k<Ng[0];k++) { if (Dvect[k]>maxDvect) {maxDvect = Dvect[k]; indmaxDvect = k; }}
	theta0 = gridtheta[indmaxDvect];
	for (k=0;k<p;k++) {if ((k+1)==j) vectr[k] = cos(theta0) * r[k] + sin(theta0); else vectr[k] = cos(theta0) * r[k];}
	norm = 0.0;
	for (k=0;k<p;k++) norm = norm + vectr[k]*vectr[k];
	for (k=0;k<p;k++) r[k] = vectr[k]/sqrt(norm);
      }
    }
      
    

    UNPROTECT(6);

    delete[] gridtheta;
    delete[] Dvect;
    delete[] vectr;
    delete[] ranking;
    delete[] resBn;
    //     delete[] Uvect;
    delete[] nptr;
    return(Rr);

  }

}

#include "deheuvels.cpp"
