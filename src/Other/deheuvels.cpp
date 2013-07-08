// Description:
// Computes the Deheuvels measure of dependence


//#include <iostream>
using namespace std;

#include <R.h>
#include "Rmath.h"
#include <R_ext/Rdynload.h>

extern "C" {
  
  void deheuvels(double *r, double *s, int *ranking, int *n, double *resBn) {
    
    void Dn(double *ss, double *tt, int *n, double *resDn);
    void rank(int *n, double *data, int *irank);

    int i, j;
    double DnR, DnS;
    double *ss, *tt, *resDn, *rsave, *ssave;
    int *irank;
    irank = new int[n[0]];
    for (i = 0; i < n[0]; i++) irank[i] = 0;
    ss = new double[1];
    tt = new double[1];
    resDn = new double[1];
    resDn[0] = 0.0;
    rsave = new double[n[0]];
    ssave = new double[n[0]];

    for (i = 0; i < n[0]; i++) {
      rsave[i] = r[i];
      ssave[i] = s[i];
    }
    
    if (ranking[0] == 1) {
           
      rank(n,r,irank);
      for (i = 0; i < n[0]; i++) r[i] = (double)irank[i];
      
      rank(n,s,irank);
      for (i = 0; i < n[0]; i++) s[i] = (double)irank[i];
            
    }
    
    resBn[0] = 0.0;
    
    for (i = 0; i < n[0]; i++) {
      
      for (j = 0; j < n[0]; j++) {
	
	ss[0] = r[i];
	tt[0] = r[j];
	Dn(ss,tt,n,resDn);
	DnR = resDn[0];
	
	ss[0] = s[i];
	tt[0] = s[j];
	Dn(ss,tt,n,resDn);
       	DnS = resDn[0];
	
	resBn[0] = resBn[0] + DnR*DnS;
	
      }
      
    }
    
    resBn[0] = resBn[0] / (double)(n[0]);
    
    for (i = 0; i < n[0]; i++) {
      r[i] = rsave[i];
      s[i] = ssave[i];
    }


    delete[] irank;
    delete[] ss;
    delete[] tt;
    delete[] resDn;
    delete[] rsave;
    delete[] ssave;

    return;

  }

  void Dn(double *ss, double *tt, int *n, double *resDn) {
    
    double max=ss[0];
    if (tt[0]>max) max = tt[0];
    
    resDn[0] = (double)(2*n[0]+1)/(double)(6*n[0]) + ss[0]*(ss[0]-1.0)/(double)(2*n[0]*(n[0]+1)) + tt[0]*(tt[0]-1.0)/(double)(2*n[0]*(n[0]+1))-max/(double)(n[0]+1);

    return;
    
  }

  void rank(int *n, double *data, int *irank) {

    void order(double* data, int* index, int *size);
    
    int i, *index;
    index = new int[n[0]];
    for(i = 0; i < n[0]; i++) index[i] = i+1;  /* set the initial index values */
    
    order(data,index,n);
    
    int j;
    
    for (j=0;j<n[0];j++) irank[index[j]-1]=j+1;
    
    delete[] index;

    return;
    
  }
  
  void order(double* data, int* index, int *size) {

    void swap(int* a, int* b);
    
    int still_changing = 1;  /* flag to know when to stop the comparisons */
    int i;
    
    while(still_changing) {
      still_changing = 0;
      for (i = 0; i < (size[0] - 1); i++) {
	
	if ( data[index[i]-1] > data[index[i + 1]-1] ) {
	  swap(index + i, index + i + 1);
	  still_changing = 1;
	}
      }
    }

    return;

  }
  
  void swap(int* a, int* b) {
    
    int temp;
    temp = *a;
    *a = *b;
    *b = temp;

    return;

  }

    
  
}


