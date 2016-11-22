/*************************************************************************
*
*                             mathpack.h
*
*  Copyright (c) 1987-1998
*      John N. Sanders-Reed
*      SVS, Inc.
*      26 Meadow View Rd.
*      Cedar Crest, NM 87008
*      505-281-8563
*  All Rights Reserved.
*
*      This software is intended for single machine use.
*      Archival copies are permitted but other reproduction or distribution
*      of this software in any form is strictly prohibited. For complete
*      details, consult the License agreement.
*
*  Header file for mathpack library: function prototypes
*
***************************************************************************
*/
#ifndef MATHPACK_H
#define MATHPACK_H

/** normally don't use the following line **/
#define MATHPACK_USE_ERF

#define FFTDAT double
#define sqr(a) (a)*(a)

/****** mtrdprsv.c ******/
double * mtread(int * nr, int * nc);
void mtprint(double * array, int nr, int nc);
void mtsave(double * array, int nr, int nc);

/******** dtrdprsv.c *******/
double * dtread(char filename[], int * n, char * tfflag);
void dtprint(double * dataset, int n, char tfflag);
int dtsave(double * dataset, int n, char tfflag);

/***** probfunc.c *****/
double inverf(double y);
double gammaln(double x);
double gammaq(double a, double x);
double q_kolmogorov(double lambda);
double betapdf(double x, double m, double n);
double cum_beta(double x, double m, double n);

#ifdef MATHPACK_USE_ERF
double erf(double x);
double erfc(double x);
#endif

/****** fftift.c *******/
/* tfflag = T for time/space forward; F for freq reverse x-form */
/* dataset in fftift(), realfft(), pwrspec() are num_pts x 2 evenly spaced
   complex data pts */
double * fftift(double *dataset, int num_pts, char tfflag);
double * realfft(double dataset[], int num_pts, char tfflag);
double * pwrspec(double *dataset, int num_pts, double dtdx);
FFTDAT * ndimfft(FFTDAT * dataset, int num_dims, int dim_len[], char tfflag);
void optoshufl(FFTDAT * dataset, int nrows, int ncols);

/****** gaussref.c ******/
int gaussref
  (
   double * array1,  /* begins as original array, ends as identify matrix */
   int n1r,          /* number of rows */
   int n1c,          /* number of columns */
   double * idntinv, /* begins as identify matrix, ends as inverse array */
   double lambda[]   /* cofactors for determinant */
  );

/******* lstsqr.c ******/
void lstsqr
  (
   double * x,     /* f1(x) f2(x) ... fn(x) array. y = Af1(x) + Bf2(x) + ... */
   int n1r,          /* number of rows (data points) in x and y */
   int n1c,          /* number of columns in x (no. of parameters) */
   double y[],       /* vector of y values */
   double coeff[],   /* vector to hold coefficients A, B, ... */
   double errcoef[], /* vector to hold uncertainty of coefficients, dA,.. */
   double erest[],   /* erest[0] = correlation coeficient, erest[1] =
			standard estimate of error of fit */
   double * corr     /* nparam x nparam correlation matrix */
  );

/***** singValD.c ******/
int singValDecomp(double *aExt, int m, int n, double w[], double *vExt);
void svdInv(double * U, double * w, double * V, int nrows, int ncols, double * inv);

/****** svdLstSqr.c *******/
void svdLstSqr(double * M, int nrows, int ncols, double * y, double * coeff,
			   double * errcoef, double * erest, double * corr);
void svdLstSqrUwV(double * M, int nrows, int ncols, double * y, double * coeff,
			   double * errcoef, double * erest, double * corr, double * U,
			   double * w, double * V);

/******* mtmult.c ******/
int mtmult(double * array1, int n1r, int n1c, double * array2, int n2r,
		   int n2c, double * array3);   /** array1 * array2 = array3 **/

/******* eigen.c *******/
void eigen
  (
   double * a,   /* begins as array whose eigenvalues are to be
		    determined. Diagonalized so version is result
		    with eigenvalues as diagonal elements. */
   int n,        /** an n x n array **/
   int ni,       /* number of iterations to be performed */
   double sums[] /* sums[0] is sum of diagonal elements, sums[1] is sums
		    of off diagonal elements */
  );

/***** cmplxmth.c *****/
// double * cmplxroot(double x[], double y[]);   /** complex sqrt(complex x) **/
// double * cadd(double n1[], double n2[], double n3[]);  /** n3 = n1 + n2 **/
// double * csub(double n1[], double n2[], double n3[]);  /** n3 = n1 - n2 **/
// double * cmult(double n1[], double n2[], double n3[]);  /** n3 = n1 * n2 **/
// double * cdiv(double n1[], double n2[], double n3[]);  /** n3 = n1 / n2 **/
// double * smult(double sclr, double cmplx[], double rslt[]);  /** rslt = sclr * cmplx; a scalar times a complex **/
// double * cpow(double cmplx[], int pwr, double rslt[]);  /** rslt = cmplx^^pwr; cmplx raised to the scalar int pwr **/

/****** cbspline.c ******/
void cubic_spline
  (
   double x[],     /* array of x values */
   double *c,      /* y values passed in, in c[][0];
		      cubic coefficients for each interval returned here */
   int n,          /* number of data points */
   int end_cond    /* derivative conditions @ end */
  );
double interp_spline
(
   double xin,    /** x value to generate interpolated y value for **/
   int n,        /** number of points in tables **/
   double * x,   /** camera values **/
   double * c    /** cubic spline coefficients (from cubic_spline()) **/
   );  /* returns -1000 is xin is out of range */

/****** norm.c *****/
void norm
  (
   double mu,          /* mean of distribution */
   double sigma, /* width of gaussian distribution (4*sigma = wdth @ 1/e^2) */
   int n,              /* number of points returned */
   double normarray[], /* array of n points, all normally distributed */
   unsigned int seed   /* initialize random number sequence: 0 to 32767 */
  );

/****** nllstsqr.c *******/
int nllstsqr
    (
     double *data,     /* (x,y) data point array */
     int n,            /* number of data points */
     double (*pf[])(double x, double * dparam, int * iparam), /* pointers to function and partial derivatives */
     int nf,           /* number of function pointers passed */
     double (*wf)(double * dparam, int * iparam, int i, int n, double * data),   /* weight function */
     double dparam[],  /* contains initial guesses, nf-1 elements, returns best
			  estimate of parameters */
	 int iparam[],     /* integer parameters (pointers) to pass to functions */
     double er[],      /* returns errors on param[] */
     double tol,       /* tolerance for chi^2 */
     char lmflag       /* L or l to use Levenberg-Marquart Method */
    );

/****** gaussfuncs.c *****/
double gaussian_mpak(double x, double * dparam, int * iparam);
double dgaussianDamp(double x, double * dparam, int * iparam);
double dgaussianDmean(double x, double * dparam, int * iparam);
double dgaussianDsigma(double x, double * dparam, int * iparam);
double dgaussianDbase(double x, double * dparam, int * iparam);
double wf_gaussian(double * param, int * iparam, int i, int n, double * data);

/****** manyroot.c ******/
void manyroot
  (
   double (*eqnarray[])(double * dparam, int * iparam),   /* pointers to functions, nvar of them */
   double (*derivarray[])(double * dparam, int * iparam), /* pointers to partial derivative functions,
                                   nvar x nvar of them */
   double dparam[],          /* double precision vector of parameters */
   int nparam,               /* number of double precision parameters */
   int iparam[],             /* integer vector of parameters (for pointers) */
   int varvec[],             /* vector of indices of variables in param[] */
   int nvar,                 /* number of variables */
   double tol,               /* tolerance on dvar/var */
   int chkvar                /* variable to check tolerance */
  );

/******* simpson.c *****/
double simpson
  (
   double ll,      /* lower limit of integration */
   double ul,      /* upper limit of integration */
   int n,          /* number of intervals */
   double (*f)(double * dparam, int * iparam),  /* function to be integrated */
   double (*df)(double * dparam, int * iparam), /* derivative of function to be integrated */
   int ni,         /* index to parameter to be integrated */
   double dparam[],/* double precision parameter vector to pass to f() and df() */
   int iparam[]    /* integer parameter vector (for pointers) to pass to f() and df() */
  );

/****** simpsn2d.c ******/
double simpsn2d
  (
   double a,        /* lower limit of integration */
   double b,        /* upper limit of integration */
   int m,           /* # of x intervals  (must be even) */
   int n,           /* # of y intervals  (must be even) */
   double (*fn)(double x, double y, double * dparam, int * iparam),  /* function to be integrated */
   double (*f1)(double x, double * dparam, int * iparam),  /* derivative */
   double (*f2)(double x, double * dparam, int * iparam),  /* derivative */
   double dparam[], /* double precision parameter vector to pass to fn(), f1(), and f2() */
   int iparam[]     /* integer parameter vector (for pointers) to pass to fn(), f1(), and f2() */
  );

/***** gausquad.c *******/
double gaussquad
  (
   double (*fn)(double * dparam, int * iparam),  /* function to be integrated */
   double ll,       /* lower limit of integration */
   double ul,       /* upper limit of integration */
   double dparam[], /* double precision parameters to pass to fn() */
   int iparam[],    /* integer parameters (for pointers) to pass to fn() */
   int ni           /* index to variable to integrate */
  );

/***** romberg.c *******/
double romberg
  (
   double (*fn)(double * dparam, int * iparam),  /* function to be integrated */
   double ll,       /* lower limit of integration */
   double ul,       /* upper limit of integration */
   double dparam[], /* double precision parameters to pass to fn() */
   int iparam[],    /* integer parameters (for pointers) to pass to fn() */
   int ni           /* index to variable to integrate */
  );

double trap_intg
  (
   double (*fn)(double * dparam, int * iparam),  /* function to be integrated */
   double ll,       /* lower limit of integration */
   double ul,       /* upper limit of integration */
   double dparam[], /* double precision parameters to pass to fn() */
   int iparam[],    /* integer parameters (for pointers) to pass to fn() */
   int ni,          /* index to variable to integrate */
   int n            /* nth refinement of extended trapezoidal rule */
  );

void interp_neville
  (
   double xa[],     /** known x and y values **/
   double ya[],
   int n,          /** number of entries in xa and ya **/
   double x,       /** x value for which y is desired **/
   double *y,      /** estimated y(x) **/
   double *dy      /** uncertainty in y value **/
  );

/***** bbfuncs.c *****/
double PlanckBB(double * dparam, int * iparam);
double PlanckBB_phot(double * dparam, int * iparam);
double PlanckBB_atten(double * dparam, int * iparam);

/**************************************************************************/
#endif  /* MATHPACK_H */
