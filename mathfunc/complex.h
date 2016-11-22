/* cmplx.h - complex arithmetic declarations */

#ifndef CMPLX_H__
#define CMPLX_H__

#include <math.h>                         /* in MSC and TC/BC, it declarares: */
                                          /* \ttt{struct complex} and \ttt{cabs(z)} */

/* struct complex{double x, y;}; */       /* uncomment if neccessary */
/* double cabs(struct complex); */        /* uncomment if neccesary */
//struct complex {double x, double y};       /* uncomment if not MSC or TC/BC */

typedef struct complex_ {
	double x; 
	double y;
} complex;

complex cmplx(double, double);            /* define complex number */
complex conjg(complex);                   /* complex conjugate */

complex cadd(complex a, complex b);           /* complex addition */
complex csub(complex a, complex b);           /* complex subtraction */
complex cmul(complex a, complex b);           /* complex multiplication */
complex cdiv(complex a, complex b);           /* complex division */

complex rmul(double, complex);            /* multiplication by real */
complex rdiv(complex, double);            /* division by real */

double real(complex);                     /* real part */
double aimag(complex);                    /* imaginary part */

complex cexp_(complex);                    /* complex exponential */

#endif /* CMPLX_H__*/
