#ifndef CINTS_H
#define CINTS_H

/*
   Low level routines for calculating Gaussian integrals.
   Use qc.h for high level use.
*/

/*************************************************************************
 *
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
 **************************************************************************/

double Bfunc(int i, int r, double g);
double contr_coulomb(int ia, double *aexps, double *acoefs,
			    double xa, double ya, double za, int la, int ma, int na, 
			    int ib, double *bexps, double *bcoefs,
			    double xb, double yb, double zb, int lb, int mb, int nb, 
			    int ic, double *cexps, double *ccoefs,
			    double xc, double yc, double zc, int lc, int mc, int nc, 
			    int id, double *dexps, double *dcoefs,
			    double xd, double yd, double zd, int ld, int md, int nd);

double coulomb_repulsion(double xa, double ya, double za,
				int la, int ma, int na, double alphaa,
				double xb, double yb, double zb,
				int lb, int mb, int nb, double alphab,
				double xc, double yc, double zc,
				int lc, int mc, int nc, double alphac,
				double xd, double yd, double zd,
				int ld, int md, int nd, double alphad);

double *B_array(int l1, int l2, int l3, int l4, double p, double a,
		double b, double q, double c, double d,
		double g1, double g2, double delta);

double B_term(int i1, int i2, int r1, int r2, int u, int l1, int l2,
		     int l3, int l4, double Px, double Ax, double Bx,
		     double Qx, double Cx, double Dx, double gamma1,
		     double gamma2, double delta);
double kinetic(double alpha1, int l1, int m1, int n1,
		      double xa, double ya, double za,
		      double alpha2, int l2, int m2, int n2,
		      double xb, double yb, double zb);
double overlap(double alpha1, int l1, int m1, int n1,
		      double xa, double ya, double za,
		      double alpha2, int l2, int m2, int n2,
		      double xb, double yb, double zb);
double overlap_1D(int l1, int l2, double PAx,
			 double PBx, double gamma);
double nuclear_attraction(double x1, double y1, double z1,
				 int l1, int m1, int n1, double alpha1,
				 double x2, double y2, double z2,
				 int l2, int m2, int n2, double alpha2,
				 double x3, double y3, double z3);
double *A_array(int l1, int l2, double PA, double PB,
		       double CP, double g);

double dist2(double x1, double y1, double z1, 
		    double x2, double y2, double z2);
double binomial_prefactor(int s, int ia, int ib, double xpa, double xpb);

int ijkl2intindex(int i, int j, int k, int l);

double product_center_1D(double alphaa, double xa, 
			 double alphab, double xb);

#endif
