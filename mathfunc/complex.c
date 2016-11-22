
#include "complex.h"
/* complex.c - complex arithmetic functions */
                         /* for MSC and TC/BC, it declares: */
                                          /* \ttt{struct complex} and \ttt{cabs()} */


                                          /* uncomment if not MS or TC/BC */
//   double cabs(complex z){
//   return sqrt(z.x * z.x + z.y * z.y);
// }



//typedef struct complex complex;

complex cmplx(double x, double y)                              /* z = cmplx(x,y) = x+jy */
{
       complex z;

       z.x = x;  z.y = y;

       return z;
}

complex conjg(complex z)                                 /* complex conjugate of z=x+jy */
{
       return cmplx(z.x, -z.y);                  /* returns z* = x-jy */
}

complex cadd(complex a, complex b)                               /* complex addition */
{
       return cmplx(a.x + b.x, a.y + b.y);
}

complex csub(complex a, complex b)                               /* complex subtraction */
{
       return cmplx(a.x - b.x, a.y - b.y);
}

complex cmul(complex a, complex b)                                /* complex multiplication */
{
       return cmplx(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

complex rmul(double a, complex  z)                               /* multiplication by real */
{
       return cmplx(a * z.x, a * z.y);
}

complex cdiv(complex a, complex b)                                /* complex division */
{
   double D = b.x * b.x + b.y * b.y;

   return cmplx((a.x * b.x + a.y * b.y) / D, (a.y * b.x - a.x * b.y) / D);
}

complex rdiv(complex z, double a)                               /* division by real */
{
       return cmplx(z.x / a, z.y / a);
}

double real(complex z)                                   /* real part Re(z) */
{
       return z.x;
}

double aimag(complex z)                                  /* imaginary part Im(z) */
{
       return z.y;
}

complex cexp_(complex z)                                  /* complex exponential */
{
       double R_ = exp(z.x);

       return cmplx(R_ * cos(z.y), R_ * sin(z.y));
}
