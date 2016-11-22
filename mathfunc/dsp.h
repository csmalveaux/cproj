
#ifndef DSP_H__
#define DSP_H__


#include <math.h>
#include <stdlib.h>
#include "complex.h"


/* \def two(x) 
   (2\sp{x}\) by left-shifting
*/
#define two(x)       (1 << (x))                  /* \(2\sp{x}\) by left-shifting */
#define eps   (1.E-9)                     /* \(\ep=10\sp{-9}\) */
#define  a_    16807                              /* that is, \(a = 7\sp{5}\) */
#define  m_    2147483647                         /* that is, \(m = 2\sp{31}-1\) */
#define  q_    127773                             /* note, \(q = m/a\) = quotient */
#define  r_    2836                               /* note, \(r = m\%a\) = remainder */   


void adc(double x, int * b, double B, int R_);
double allpass(int D, double * w, double ** p, double a, double x);
int bitrev(int n, int B);
void blockcon(int M, double * h, int L, double * x, double * y, double * ytemp);
double can(int M, double * a, int L, double * b, double * w, double x);
double can2(int M, double *a, int L, double * b, double * w, double x);
double can3(int M, double * a, double * b, double * w, double x);
double cas(int K, double ** A, double ** B, double ** W, double x);
void cas2can(int K, double ** A, double * a);
double ccan(int M, double * a, double * b, double * w, double ** p, double x);
double ccan2(int M, double * a, double * b, double * w,int * q, double x);
double ccas(int K, double ** A, double ** B, double ** W, double ** P, double x);
double ccas2(int K, double ** A, double ** B, double ** W, int * Q, double x);
void cdelay(int D, double * w, double ** p);
void cdelay2(int D, int * q);
double cfir(int M, double * h, double * w, double ** p, double x);
double cfir1(int M, double * h, double * w, double ** p, double x);
double cfir2(int M, double * h, double * w, int * q, double x);
double cheby(int N, double x);
void conv(int M, double * h, int L, double * x, double * y);
void corr(int N, double * x, double * y, int M, double * R_);
double csos(double * a, double * b, double * w, double ** p, double x);
double csos2(double * a, double * b, double * w, int * q, double x);
double dac(int * b, int B, double R_);
void delay(int D, double * w);
double delta(int n);
void dft(int L, double * x, int N, complex * X);
void dftmerge(int N, complex * XF);
double dir(int M, double * a, int L, double * b, double * w, double * v, double x);
double dir2(int M, double * a, int L, double * b, double * w, double * v, double x);
double dot(int M, double * h, double * w);
complex dtft(int L, double * x, double w);
void dtftr(int L, double * x, int N, complex * X, double wa, double wb);
void fft(int N, complex * X);
double fir(int M, double * h, double * w, double x) ;
double fir2(int M, double * h, double * w, double x);
double fir3(int M, double * h, double * w, double x);
void gdelay2(int D, double c, double * q);
double gran(double m, double s, long * iseed);
double I0(double x);
void ifft(int N, complex * X);
double lowpass(int D, double * w, double ** p, int M, double * a, double * b, double * v, double x);
void modwrap(int L, double * x, int N, double * xtilde);
double plain(int D, double * w, double ** p, double a, double x);
double ran(long * iseed);
double ran1f(int B, double * u, int * q, long * iseed);
double ranh(int D, double * u, int * q, long * iseed);
double ranl(int D, double * u, int * q, long * iseed);
void shuffle(int N, complex * X);
double sine(int D, int i);
double sos(double * a, double * b, double * w, double x);
double square(int D1, int i);
void swap(complex * a, complex * b);
double tap(int D, double * w, double * p, int i);
double tap2(int D, double * w, int q, int i);
double tapi(int D, double * w, double * p, double d);
double tapi2(int D, double * w, int q, double d);
double trapez(int D, int  D1, int  D2, int i);
int u(double x);
double wavgen(int D, double * w, double A, double F, double * q);
double wavgeni(int D, double * w, double  A, double F, double * q);
double wavgenr(int D, double * w, double A, double F, double * q);
void wrap(int M, double * w, double ** p);
void wrap2(int M, int * q);

#endif /* DSP_H__*/
