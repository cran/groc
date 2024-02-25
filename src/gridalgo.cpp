// Description
// Performs the grid algorithm from Croux (2007)

// g++ -m64 -I/usr/include/R  -I/usr/local/include -fpic  -g -c gridalgo.cpp -o gridalgo.o
// g++ -m64 -shared -L/usr/local/lib64 -o gridalgo.so gridalgo.o -L/usr/lib64/R/lib -lR

#include <iostream>
using namespace std;

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

extern "C" {


  SEXP  grid(SEXP RNc, SEXP RNg, SEXP Rr, SEXP Rs, SEXP RU, SEXP RY, SEXP Rmaxiter, SEXP RD, SEXP rho) {

    if(!isFunction(RD) & (RD != R_NilValue)) error("RD must be a function");
    if(!isEnvironment(rho)) error("rho must be an environment");

    R_len_t p = length(Rr), q = length(Rs), n = nrows(RU);

    SEXP Rytmp, Rxtmp;

    int *maxiter;
    double *r, *s, *U, *Y;

    PROTECT(RNc = AS_INTEGER(RNc));
    PROTECT(RNg = AS_INTEGER(RNg));
    PROTECT(Rmaxiter = AS_INTEGER(Rmaxiter));

    PROTECT_INDEX ipr;
    PROTECT_WITH_INDEX(Rr = coerceVector(Rr, REALSXP), &ipr);

    PROTECT_INDEX ips;
    PROTECT_WITH_INDEX(Rs = coerceVector(Rs, REALSXP), &ips);

    //    PROTECT(Rr = coerceVector(Rr, REALSXP));
    //    PROTECT(Rs = coerceVector(Rs, REALSXP));
    PROTECT(RU = coerceVector(RU, REALSXP));
    PROTECT(RY = coerceVector(RY, REALSXP));

    maxiter = INTEGER(Rmaxiter);
    r = REAL(Rr);
    s = REAL(Rs);
    U = REAL(RU);
    Y = REAL(RY);

    SEXP simplegrid(SEXP Rr, SEXP RU, SEXP RY, SEXP RNg, SEXP RNc, SEXP RD, SEXP rho);
    void deheuvels(double *r, double *s, int *ranking, int *n, double *resBn);

    int i, k, j, l;
    double *ytmp, *xtmp;
    ytmp = new double[n];
    xtmp = new double[n];

    PROTECT(Rytmp = allocVector(REALSXP, n));
    PROTECT(Rxtmp = allocVector(REALSXP, n));
    ytmp = REAL(Rytmp);
    xtmp = REAL(Rxtmp);

    SEXP chariter;
    PROTECT(chariter = allocVector(STRSXP, 1));
    SET_STRING_ELT(chariter, 0, mkChar("failed"));

    double error;
    double *rold, *sold;
    rold = new double[p];
    sold = new double[q];
    double maxs, maxr, tmpr, tmps;

    int nbiter;

    if (q==1) {
      r[0] = 1.0;
      for (i=1;i<p;i++) r[i] = 0.0;
      REPROTECT(Rr = simplegrid(Rr,RU,RY,RNg,RNc,RD,rho), ipr);
      //      Rr = coerceVector(Rr, REALSXP), ipr);
      r = REAL(Rr);
    }

    if (q>1) {

      r[0] = 1.0;
      for (i=1;i<p;i++) r[i] = 0.0;
      s[0] = 1.0;
      for (i=1;i<q;i++) s[i] = 0.0;

      error = 1.0;
      nbiter = 1;
      while(error>0.0001) {

	if (nbiter > maxiter[0]) {
	  delete[] rold;
	  delete[] sold;
	  UNPROTECT(10);
	  return(chariter);
	}
					  
        for (k=0;k<p;k++) rold[k] = r[k];
	for (k=0;k<q;k++) sold[k] = s[k];
	
	for (j=0;j<n;j++) { 
	  ytmp[j] = 0.0;
	  for (l=0;l<q;l++) ytmp[j] = ytmp[j] + Y[l*n+j]*s[l];
	}

	//	REPROTECT(Rr = simplegrid(Rr,RU,Rytmp,RNg,RNc,RD,rho), ipr);
	Rr = simplegrid(Rr,RU,Rytmp,RNg,RNc,RD,rho);
	//	REPROTECT(Rr = coerceVector(Rr, REALSXP), ipr);
	r = REAL(Rr);
	
	for (j=0;j<n;j++) { 
	  xtmp[j] = 0.0;
	  for (l=0;l<p;l++) xtmp[j] = xtmp[j] + U[l*n+j]*r[l];
	}

	REPROTECT(Rs = simplegrid(Rs,RY,Rxtmp,RNg,RNc,RD,rho), ips);	  
	//	REPROTECT(Rs = coerceVector(Rs, REALSXP), ips);
	s = REAL(Rs);
	
	maxr = 0.0;
	maxs = 0.0;
	for (k=0;k<p;k++) {
	  tmpr = fabs(rold[k]-r[k]);
	  if (tmpr > maxr) maxr = tmpr;
	}
	for (k=0;k<q;k++) {
	  tmps = fabs(sold[k]-s[k]);
	  if (tmps > maxs) maxs = tmps;
	}
	error = maxs;
	if (maxr > maxs) error = maxr;

	nbiter = nbiter + 1;

      }
    }
    
    delete[] rold;
    delete[] sold;
    UNPROTECT(10);

    return(Rr);

  }

  
}

#include "Other/simplegrid.cpp"

