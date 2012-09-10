#ifndef QC_H
#define QC_H

/*
    High level functions for Gaussian integrals:

    getS
    getT
    getV
    getInts  // Two particle integrals
*/

void getS(int nbf, int *nprim, int *istart,
	  double *xcenter, double *ycenter, double *zcenter,
	  int *lpower,int *mpower, int *npower,
      int n2,
	  double *coef, double *alpha,
	  double *S);

void getT(int nbf, int *nprim, int *istart,
	  double *xcenter, double *ycenter, double *zcenter,
	  int *lpower,int *mpower, int *npower,
      int n2,
	  double *coef, double *alpha,
	  double *T);

void getV(int nbf, int *nprim, int *istart,
	  double *xcenter, double *ycenter, double *zcenter,
	  int *lpower,int *mpower, int *npower,
      int n2,
	  double *coef, double *alpha,
	  int nat, int *atno, double *x, double *y, double *z,
	  double *V);

void getInts(int nbf, int *nprim, int *istart,
	     double *xcenter, double *ycenter, double *zcenter,
	     int *lpower,int *mpower, int *npower,
         int n2,
	     double *coef, double *alpha, double *Ints);

#endif
