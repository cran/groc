// Description:
// Computes the dcov measure of dependence
// Distance correlation from Gabor et al., Annals of Stat, 2007, vol 35 (6), p.2769-2794


#include <iostream>
using namespace std;

#include <R.h>
#include "Rmath.h"
#include <R_ext/Rdynload.h>

extern "C" {
  
  void dcov(double *x, double *y, int *xlen, double *Vup) {

    int n=xlen[0], i, j;
    double meana=0.0, meanb=0.0;
    double *a, *b;
    a = new double[n*n];
    b = new double[n*n];
    double *meanak, *meanal;
    meanak = new double[n];
    meanal = new double[n];
    double *meanbk, *meanbl;
    meanbk = new double[n];
    meanbl = new double[n];
    
    for (i=1 ; i<=n ; i++) {
      for (j=1 ; j<=n ; j++) {
	a[(j-1)*n+(i-1)] = fabs(x[i-1]-x[j-1]);
	b[(j-1)*n+(i-1)] = fabs(y[i-1]-y[j-1]);
      }
    }

    for (i=1 ; i<=n ; i++) {
      meanak[i-1] = 0.0;
      meanbk[i-1] = 0.0;
    }
    for (i=1 ; i<=n ; i++) {
      for (j=1 ; j<=n ; j++) {
	meanak[i-1] = meanak[i-1] + a[(j-1)*n+(i-1)];
	meanbk[i-1] = meanbk[i-1] + b[(j-1)*n+(i-1)];
      }
      meanak[i-1] = meanak[i-1]/n;
      meanbk[i-1] = meanbk[i-1]/n;
    }

    for (i=1 ; i<=n ; i++) {
      meanal[i-1] = 0.0;
      meanbl[i-1] = 0.0;
    }
    for (i=1 ; i<=n ; i++) {
      for (j=1 ; j<=n ; j++) {
	meanal[i-1] = meanal[i-1] + a[(i-1)*n+(j-1)];
	meanbl[i-1] = meanbl[i-1] + b[(i-1)*n+(j-1)];
      }
      meanal[i-1] = meanal[i-1]/n;
      meanbl[i-1] = meanbl[i-1]/n;
    }

    for (i=1 ; i<=n ; i++) {
      for (j=1 ; j<=n ; j++) {
	meana = meana + a[(j-1)*n+(i-1)];
	meanb = meanb + b[(j-1)*n+(i-1)];
      }
    }
    meana = meana/n/n;
    meanb = meanb/n/n;

    for (i=1 ; i<=n ; i++) {
      for (j=1 ; j<=n ; j++) {
	a[(j-1)*n+(i-1)] = a[(j-1)*n+(i-1)] - meanak[i-1] - meanal[j-1] + meana;
	b[(j-1)*n+(i-1)] = b[(j-1)*n+(i-1)] - meanbk[i-1] - meanbl[j-1] + meanb;
      }
    }

    Vup[0] = 0.0;
    for (i=1 ; i<=n ; i++) {
      for (j=1 ; j<=n ; j++) {
	Vup[0] = Vup[0] + a[(j-1)*n+(i-1)]*b[(j-1)*n+(i-1)];
      }
    }

    delete[] a;
    delete[] meanak;
    delete[] meanal;
    delete[] b;
    delete[] meanbk;
    delete[] meanbl;

    return;

  }
  
}


