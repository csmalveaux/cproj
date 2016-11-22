

#include "dsp.h"  

/* Prototypes */

void adc(double x, int * b, double B, int R_);
double allpass(int D, double * w, double ** p, double a, double x);
int bitrev(int n, int B);
void blockcon(int M, double * h, int L, double * x, double * y, double * ytemp);
double can(int M, double * a, int L, double * b, double * w, double x);
double can2(int M, double *a, int L, double * b, double * w, double x);

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

/* \fn void adc(double x, int * b, double B, int R)
   \brief successive approximation A/D converter
   \param double x
   \param int * b
   \param double B
   \param int R
*/
void adc(double x, int * b, double B, int R_){
       int i;
       double y, xQ, Q;

       Q = R_ / pow(2, B);                        /* quantization width \(Q=R/2\sp{B}\) */
       y = x + Q/2;                              /* rounding */

       for (i = 0; i < B; i++)                   /* initialize bit vector */
              b[i] = 0;

       b[0] = 1 - u(y);                          /* determine MSB */

       for (i = 1; i < B; i++) {                 /* loop starts with \(i=1\) */
              b[i] = 1;                          /* turn \(i\)th bit ON */
              xQ = dac(b, B, R_);                 /* compute DAC output */
              b[i] = u(y-xQ);                    /* test and correct bit */
              }
}

/* \fn double allpass(int D, double * w, double ** p, double a, double x
   \brief allpass reverberator with circular delay line
   \param 
   \return
*/
double allpass(int D, double * w, double ** p, double a, double x){                  /* usage: y=allpass(D,w,&p,a,x); */
       double y, s0, sD;

       sD = tap(D, w, *p, D);                   /* \(D\)th tap delay output */
       s0 = x + a * sD;
       y  = -a * s0 + sD;                       /* filter output */
       **p = s0;                                /* delay input */
       cdelay(D, w, p);                         /* update delay line */

       return y;
}

/* \fn int bitrev(int n, int B)
   \brief bit reverse of a B-bit integer n 
   \param
   \return
*/
int bitrev(int n, int B){
       int m, r;

       for (r=0, m=B-1; m>=0; m--)
          if ((n >> m) == 1) {                   /* if \(2\sp{m}\) term is present, then */
             r += two(B-1-m);                    /* add \(2\sp{B-1-m}\) to \(r\), and */
             n -= two(m);                        /* subtract \(2\sp{m}\) from \(n\) */
             }

       return(r);
}

/* \fn void blockcon(int M, double * h, int L, double * x, double * y, double * ytemp)
   \brief block convolution by overlap-add method
   \param
   \return
*/
void blockcon(int M, double * h, int L, double * x, double * y, double * ytemp)
//double *h, *x, *y, *ytemp;                    /* ytemp is tail of previous block */
//int M, L;                                     /* \(M\) = filter order, \(L\) = block size */
{
    int i;

    conv(M, h, L, x, y);                      /* compute output block y */

    for (i=0; i<M; i++) {
        y[i] += ytemp[i];                     /* add tail of previous block */
        ytemp[i] = y[i+L];                    /* update tail for next call */
        }
}

/* \fn double can(int M, double * a, int L, double * b, double * w, double x) 
   \brief IIR filtering in canonical form
   \param
   \return
*/
double can(int M, double * a, int L, double * b, double * w, double x)              /* usage: y = can(M, a, L, b, w, x); */
//double *a, *b, *w, x;                     /* \(w\) = internal state vector */
//int M, L;                                 /* denominator and numerator orders */
{
       int K, i;
       double y = 0;

       K = (L <= M) ? M : L;              /* \(K=\max(M,L)\) */

       w[0] = x;                          /* current input sample */

       for (i=1; i<=M; i++)               /* input adder */
              w[0] -= a[i] * w[i];

       for (i=0; i<=L; i++)               /* output adder */
              y += b[i] * w[i];

       for (i=K; i>=1; i--)               /* reverse updating of \(w\) */
              w[i] = w[i-1];

       return y;                          /* current output sample */
}

/* \fn double can(int M, double * a, int L, double * b, double * w, double x) 
   \brief IIR filtering in canonical form
   \param
   \return
*/
double can2(int M, double *a, int L, double * b, double * w, double x)             /* usage: y = can2(M, a, L, b, w, x); */
//double *a, *b, *w, x;
//int M, L;
{
       int K;
       double y;

       K = (L <= M) ? M : L;                     /* \(K=\max(M,L)\) */

       w[0] = 0;                                 /* needed for dot(M,a,w) */

       w[0] = x - dot(M, a, w);                  /* input adder */

       y = dot(L, b, w);                         /* output adder */

       delay(K, w);                              /* update delay line */

       return y;                                 /* current output sample */
}

/* \fn double can3(int M, double * a, double * b, double * w, double x)  
   \brief IIR filtering in canonical form, emulating a DSP chip
   \param double *a, *b, *w, x;                      \(w\) = internal state vector 
   \param int M;                                     \(a,b\) have common order \(M\)
   \return
*/
double can3(int M, double * a, double * b, double * w, double x)                /* usage: y = can3(M, a, b, w, x); */
{
       int i;
       double y;

       w[0] = x;                                 /* read input sample */

       for (i=1; i<=M; i++)                      /* forward order */
              w[0] -= a[i] * w[i];               /* MAC instruction */

       y = b[M] * w[M];

       for (i=M-1; i>=0; i--) {                  /* backward order */
              w[i+1] = w[i];                     /* data shift instruction */
              y += b[i] * w[i];                  /* MAC instruction */
              }

       return y;                                 /* output sample */
}

/* \fn double cas(int K, double ** A, double ** B, double ** W, double x)
   \brief IIR filtering in cascade of second-order sections
   \param int K
   \param int K
   \return double **A, **B, **W, x;                          \(A,B,W\) are \(Kx3\) matrices 
*/
double cas(int K, double ** A, double ** B, double ** W, double x){
       int i;
       double y;

       y = x;                                    /* initial input to first SOS */

       for (i=0; i<K; i++)
              y = sos(A[i], B[i], W[i], y);      /* output of \(i\)th section */

       return y;                                 /* final output from last SOS */
}


/* \fn void cas2can(int K, double ** A, double * a)
   \brief cascade to canonical
   \param double **A, *a;                                          \(A\) is \(Kx3\) matrix
   \param int K;                                                   \(K\) = no. of sections
   \return
*/
void cas2can(int K, double ** A, double * a)                                   /* \(a\) is \((2K+1)\)-dimensional */
{
       int i,j;
       double *d;

       d = (double *) calloc(2*K+1, sizeof(double));

       a[0] = 1;                                        /* initialize */

       for(i=0; i<K; i++) {
              conv(2, A[i], 2*i+1, a, d);               /* \(d = a[i] \ast a\) */
              for(j=0; j<2*i+3; j++)                    /* \(a = d\) */
                     a[j] = d[j];
              }
}


/* \fn double ccan(int M, double * a, double * b, double * w, double ** p, double x)
   \brief circular buffer implementation of canonical realization
   \param double *a, *b, *w, **p, x;                 \(p\) = circular pointer to buffer \(w\) 
   \param int M;                                     \(a,b\) have common order \(M\) 
   \return
*/
double ccan(int M, double * a, double * b, double * w, double ** p, double x)             /* usage: y=ccan(M, a, b, w, &p, x); */
{
       int i;
       double y = 0, s0;

       **p = x;                                 /* read input sample \(x\) */

       s0  = *(*p)++;                           /* \(s\sb{0}=x\) */
       wrap(M, w, p);                           /* \(p\) now points to \(s\sb{1}\) */

       for (a++, i=1; i<=M; i++) {              /* start with \(a\) incremented to \(a\sb{1}\) */
              s0 -= (*a++) * (*(*p)++);
              wrap(M, w, p);
              }

       **p = s0;                                /* \(p\) has wrapped around once */

       for (i=0; i<=M; i++) {                   /* numerator part */
              y += (*b++) * (*(*p)++);
              wrap(M, w, p);                    /* upon exit, \(p\) has wrapped */
              }                                 /* around once again */

       (*p)--;                                  /* update circular delay line */
       wrap(M, w, p);

       return y;                                /* output sample */
}

/* \fn double ccan2(int M, double * a, double * b, double * w,int * q, double x)
   \brief circular buffer implementation of canonical realization
   \param double *a, *b, *w, x;                 q = circular pointer offset index 
   \param int M, *q;                            a,b have common order M 
   \return
*/
double ccan2(int M, double * a, double * b, double * w,int * q, double x)
{
       int i;
       double y = 0;

       w[*q] = x;                               /* read input sample x */

       for (i=1; i<=M; i++)
              w[*q] -= a[i] * w[(*q+i)%(M+1)];

       for (i=0; i<=M; i++)
              y += b[i] * w[(*q+i)%(M+1)];

       (*q)--;                                  /* update circular delay line */
       wrap2(M, q);

       return y;                                /* output sample */
}

/* \fn double ccas(int K, double ** A, double ** B, double ** W, double ** P, double x)
   \brief circular buffer implementation of cascade realization
   \param double **A, **B, **W, **P, x;            \(P\) = array of circular pointers 
   \param int K
   \return
*/
double ccas(int K, double ** A, double ** B, double ** W, double ** P, double x){
       int i;
       double y;

       y = x;

       for (i=0; i<K; i++)
              y = csos(A[i], B[i], W[i], P+i, y);        /* note, \(P+i\) = &\(P[i]\) */

       return y;
}

/* \fn double ccas2(int K, double ** A, double ** B, double ** W, int * Q, double x)
   \brief circular buffer implementation of cascade realization
   \param int K, *Q;                            \(Q\) = array of circular pointer offsets 
   \param double **A, **B, **W, x; 
   \return
*/
double ccas2(int K, double ** A, double ** B, double ** W, int * Q, double x){
       int i;
       double y;

       y = x;

       for (i=0; i<K; i++)
              y = csos2(A[i], B[i], W[i], Q+i, y);       /* note, \(Q+i\) = &\(Q[i]\) */

       return y;
}

/* \fn void cdelay(int D, double * w, double ** p)
   \brief circular buffer implementation of D-fold delay
   \param int D;
   \param double *w, **p;
   \return
*/
void cdelay(int D, double * w, double ** p){
       (*p)--;                      /* decrement pointer and wrap modulo-\((D+1)\) */
       wrap(D, w, p);               /* when \(*p=w-1\), it wraps around to \(*p=w+D\) */
}


/* \fn void cdelay2(int D, int * q)
   \brief circular buffer implementation of D-fold delay
   \param int D, *q;
   \return
*/
void cdelay2(int D, int * q){
       (*q)--;                      /* decrement offset and wrap modulo-\((D+1)\) */
       wrap2(D, q);                 /* when \(*q=-1\), it wraps around to \(*q=D\) */
}

/* \fn double cfir(int M, double * h, double * w, double ** p, double x)
   \brief FIR filter implemented with circular delay-line buffer
   \param double *h, *w, **p, x;                          \(p\) = circular pointer to \(w\) 
   \param int M;                                          \(M\) = filter order 
   \return
*/
double cfir(int M, double * h, double * w, double ** p, double x){                        
       int i;
       double y;

       **p = x;                                /* read input sample \(x\) */

       for (y=0, i=0; i<=M; i++) {             /* compute output sample \(y\) */
              y += (*h++) * (*(*p)++);
              wrap(M, w, p);
              }

       (*p)--;                                 /* update circular delay line */
       wrap(M, w, p);

       return y;
}

/* \fn double cfir1(int M, double * h, double * w, double ** p, double x)
   \brief FIR filter implemented with circular delay-line buffer
   \param double *h, *w, **p, x;
   \param int M;
   \return
*/
double cfir1(int M, double * h, double * w, double ** p, double x){                        
       int i;
       double y;

       *(*p)-- = x;
       wrap(M, w, p);                          /* \(p\) now points to \(s\sb{M}\) */

       for (y=0, h+=M, i=M; i>=0; i--) {       /* \(h\) starts at \(h\sb{M}\) */
              y += (*h--) * (*(*p)--);
              wrap(M, w, p);
              }

       return y;
}

/* \fn double cfir2(int M, double * h, double * w, int * q, double x)
   \brief FIR filter implemented with circular delay-line buffer
   \param double *h, *w, x;                                 \(q\) = circular offset index 
   \param int M, *q;                                        \(M\) = filter order 
   \return
*/
double cfir2(int M, double * h, double * w, int * q, double x){                        
       int i;
       double y;

       w[*q] = x;                                /* read input sample \(x\) */

       for (y=0, i=0; i<=M; i++) {               /* compute output sample \(y\) */
              y += (*h++) * w[(*q)++];
              wrap2(M, q);
              }

       (*q)--;                                   /* update circular delay line */
       wrap2(M, q);
       
       return y;
}

/* \fn double cheby(int N, double x) 
   \brief Chebyshev polynomial \(C\sb{N}(x)\)
   \param int N;                                            \(N\) = polynomial order 
   \param double x;                                        \(x\) must be non-negative
   \return
*/
double cheby(int N, double x)                               /* usage: y = cheby(N, x); */
{
       if (x <= 1)
              return cos(N * acos(x));
       else
              return cosh(N * log(x + sqrt(x*x-1)));
}

/* \fn void conv(int M, double * h, int L, double * x, double * y)
   \brief convolution of x[n] with h[n], resulting in y[n]
   \param double *h, *x, *y;                         \(h,x,y\) = filter, input, output arrays 
   \param int M, L;                                  \(M\) = filter order, \(L\) = input length 
   \return
*/
void conv(int M, double * h, int L, double * x, double * y){
       int n, m;

       for (n = 0; n < L+M; n++)
              for (y[n] = 0, m = MAX(0, n-L+1); m <= MIN(n, M); m++)
                     y[n] += h[m] * x[n-m];
}


/* \fn void corr(int N, double * x, double * y, int M, double * R)
   \brief sample cross correlation of two length-N signals
   \param double *x, *y, *R;                         \(x,y\) are \(N\)-dimensional 
   \param int N, M;                                  \(R\) is \((M+1)\)-dimensional 
   \return
*/
void corr(int N, double * x, double * y, int M, double * R_)                  /* computes \(R[k]\), \(k = 0, 1,\dotsc, M\) */
{
       int k, n;

       for (k=0; k<=M; k++)
              for (R_[k]=0, n=0; n<N-k; n++)
                     R_[k] += x[n+k] * y[n] / N;
}

/* \fn double csos(double * a, double * b, double * w, double ** p, double x)  
   \brief circular buffer implementation of a single SOS
   \param double *a, *b, *w, **p, x;                       \(p\) is circular pointer to \(w\) 
   \return
*/
double csos(double * a, double * b, double * w, double ** p, double x)                      /* \(a,b,w\) are 3-dimensional */
{
       double y, s0;

       *(*p) = x;                               /* read input sample \(x\) */

       s0  = *(*p)++;           wrap(2, w, p);
       s0 -= a[1] * (*(*p)++);  wrap(2, w, p);
       s0 -= a[2] * (*(*p)++);  wrap(2, w, p);

       *(*p) = s0;                              /* \(p\) has wrapped around once */

       y  = b[0] * (*(*p)++);  wrap(2, w, p);
       y += b[1] * (*(*p)++);  wrap(2, w, p);
       y += b[2] * (*(*p));                     /* \(p\) now points to \(s\sb{2}\) */

       return y;
}

/* \fn double csos2(double * a, double * b, double * w, int * q, double x)
   \brief circular buffer implementation of a single SOS
   \param double *a, *b, *w, x;                    \(a,b,w\) are 3-dimensional arrays
   \param int *q;                                  \(q\) is circular offset relative to \(w\) 
   \return
*/
double csos2(double * a, double * b, double * w, int * q, double x){
       double y;

       w[*q] = x - a[1] * w[(*q+1)%3] - a[2] * w[(*q+2)%3];

       y = b[0] * w[*q] + b[1] * w[(*q+1)%3] + b[2] * w[(*q+2)%3];

       (*q)--;  
       wrap2(2, q);

       return y;
}

/* \fn double dac(int * b, int B, double R)
   \brief bipolar two's complement D/A converter
   \param int *b, B;                          bits are dimensioned as \(b[0], b[1], \dotsc, b[B-1]\) 
   \param double R;
   \return
*/
double dac(int * b, int B, double R_){
       int i;
       double dac = 0;

       b[0] = 1 - b[0];                          /* complement MSB */

       for (i = B-1; i >= 0; i--)                /* H\"orner's rule */
          dac = 0.5 * (dac + b[i]);

       dac = R_ * (dac - 0.5);                    /* shift and scale */

       b[0] = 1 - b[0];                          /* restore MSB */

       return dac;
}

/* \fn void delay(int D, double * w)
   \brief delay by D time samples
   \param int D;
   \param double *w;
   \return
*/
void delay(int D, double * w)                          /* \(w[0]\) = input, \(w[D]\) = output */
{
       int i;

       for (i=D; i>=1; i--)               /* reverse-order updating */
              w[i] = w[i-1];

}

/* \fn double delta(int n)
   \brief delta function
   \param int n;
   \return
*/
double delta(int n){
       if (n == 0)
              return 1;
       else
              return 0;
}

/* \fn void dft(int L, double * x, int N, complex * X)
   \brief N-point DFT of length-L real-valued signal
   \param double *x;                                        \(x\) is \(L\)-dimensional real 
   \param complex *X;                                       \(X\) is \(N\)-dimensional complex 
   \param int L, N;
   \return
*/
void dft(int L, double * x, int N, complex * X)                             /* usage: dft(L, x, N, X); */
{
       double pi = 4 * atan(1.0);

       dtftr(L, x, N, X, 0.0, 2*pi);             /* \(N\) frequencies over \([0,2\pi)\) */
}   

/* \fn void dftmerge(int N, complex * XF)
   \brief DFT merging for radix 2 decimation-in-time FFT
   \param complex *XF;
   \param int N;
   \return
*/
void dftmerge(int N, complex * XF){
       double pi = 4. * atan(1.0);
       int k, i, p, q, M;
       complex  A, B, V, W;

       M = 2;
       while (M <= N) {                          /* two \((M/2)\)-DFTs into one \(M\)-DFT  */
            W = cexp_(cmplx(0.0, -2 * pi / M));   /* order-\(M\) twiddle factor */
            V = cmplx(1.0, 0.0);                   /* successive powers of \(W\) */
            for (k = 0; k < M/2; k++) {          /* index for an \((M/2)\)-DFT */
                 for (i = 0; i < N; i += M) {    /* \(i\)th butterfly; increment by \(M\) */
                      p = k + i;                 /* absolute indices for */
                      q = p + M / 2;             /* \(i\)th butterfly */
                      A = XF[p];
                      B = cmul(XF[q], V);        /* \(V = W\sp{k}\) */
                      XF[p] = cadd(A, B);        /* butterfly operations */
                      XF[q] = csub(A, B);
                      }
                 V = cmul(V, W);                 /* \(V = VW = W\sp{k+1}\) */
                 }
            M = 2 * M;                           /* next stage */
            }
}

/* \fn double dir(int M, double * a, int L, double * b, double * w, double * v, double x)
   \brief IIR filtering in direct form 
   \param double *a, *b, *w, *v, x;                  \(v,w\) are internal states 
   \param int M, L;                                  denominator and numerator orders 
   \return
*/
double dir(int M, double * a, int L, double * b, double * w, double * v, double x)           /* usage: y = dir(M, a, L, b, w, v, x); */
{
       int i;

       v[0] = x;                          /* current input sample */
       w[0] = 0;                          /* current output to be computed */

       for (i=0; i<=L; i++)               /* numerator part */
              w[0] += b[i] * v[i];

       for (i=1; i<=M; i++)               /* denominator part */
              w[0] -= a[i] * w[i];

       for (i=L; i>=1; i--)               /* reverse-order updating of \(v\) */
              v[i] = v[i-1];

       for (i=M; i>=1; i--)               /* reverse-order updating of \(w\) */
              w[i] = w[i-1];

       return w[0];                       /* current output sample */
}

/* \fn double dir2(int M, double * a, int L, double * b, double * w, double * v, double x)
   \brief IIR filtering in direct form
   \param double *a, *b, *w, *v, x;
   \param int M, L;
   \return
*/
double dir2(int M, double * a, int L, double * b, double * w, double * v, double x)          /* usage: y = dir2(M, a, L, b, w, v, x); */
{
       v[0] = x;                                 /* current input sample */
       w[0] = 0;                                 /* needed for dot(M,a,w) */

       w[0] = dot(L, b, v) - dot(M, a, w);       /* current output */

       delay(L, v);                              /* update input delay line */

       delay(M, w);                              /* update output delay line */

       return w[0];
}


/* \fn double dot(int M, double * h, double * w)
   \brief dot product of two length-(M+1) vectors
   \param double *h, *w;                                  \(h\) = filter vector, \(w\) = state vector 
   \param int M;                                          \(M\) = filter order 
   \return
*/
double dot(int M, double * h, double * w)                            /* Usage: y = dot(M, h, w); */
{                        
       int i;
       double y;

       for (y=0, i=0; i<=M; i++)               /* compute dot product */
              y += h[i] * w[i];      

       return y;
}

/* \fn complex dtft(int L, double * x, double w)
   \brief DTFT of length-L signal at a single frequency w
   \param double *x, w;                                         \(x\) is \(L\)-dimensional 
   \param int L;
   \return
*/
complex dtft(int L, double * x, double w)                                /* usage: X=dtft(L, x, w); */
{
       complex z, X;
       int n;

       z = cexp_(cmplx(0, -w));                       /* set \(z=e\sp{-j\om}\) */

       X = cmplx(0,0);                               /* initialize \(X=0\) */

       for (n=L-1; n>=0; n--)
              X = cadd(cmplx(x[n], 0), cmul(z, X));

       return X;
}

/* \fn void dtftr(int L, double * x, int N, complex * X, double wa, double wb)
   \brief N DTFT values over frequency range [wa, wb)
   \param double *x, wa, wb;                                  \(x\) is \(L\)-dimensional real 
   \param complex *X;                                         \(X\) is \(N\)-dimensional complex 
   \param int L, N;
   \return
*/
void dtftr(int L, double * x, int N, complex * X, double wa, double wb)                     /* usage: dtftr(L, x, N, X, wa, wb); */
{
       int k;
       double dw = (wb-wa)/N;                      /* frequency bin width */

       for (k=0; k<N; k++)
              X[k] = dtft(L, x, wa + k*dw);        /* \(k\)th DTFT value \(X(\om\sb{k})\) */
}  

/* \fn void fft(int N, complex * X)
   \brief in-place decimation-in-time FFT
   \param complex *X;
   \param int N
   \return
*/
void fft(int N, complex * X){
       shuffle(N, X);
       dftmerge(N, X);
}   

/* \fn double fir(int M, double * h, double * w, double x) 
   \brief FIR filter in direct form
   \param double *h, *w, x;                             \(h\) = filter, \(w\) = state, \(x\) = input sample 
   \param int M;                                        \(M\) = filter order 
   \return
*/
double fir(int M, double * h, double * w, double x)                       /* Usage: y = fir(M, h, w, x); */
{                        
       int i;
       double y;                             /* output sample */

       w[0] = x;                             /* read current input sample \(x\) */

       for (y=0, i=0; i<=M; i++)
              y += h[i] * w[i];              /* compute current output sample \(y\) */

       for (i=M; i>=1; i--)                  /* update states for next call */
              w[i] = w[i-1];                 /* done in reverse order */

       return y;
}

/* \fn double fir2(int M, double * h, double * w, double x)
   \brief FIR filter in direct form
   \param double *h, *w, x;                              \(h\) = filter, \(w\) = state, \(x\) = input 
   \param int M;                                         \(M\) = filter order 
   \return
*/
double fir2(int M, double * h, double * w, double x)                       /* Usage: y = fir2(M, h, w, x); */
{                        
       double y;

       w[0] = x;                              /* read input */

       y = dot(M, h, w);                      /* compute output */

       delay(M, w);                           /* update states */

       return y;
}

/* \fn double fir3(int M, double * h, double * w, double x)
   \brief FIR filter emulating a DSP chip
   \param double *h, *w, x;
   \param int M;
   \return
*/
double fir3(int M, double * h, double * w, double x){                        
       int i;
       double y;

       w[0] = x;                                 /* read input */

       for (y=h[M]*w[M], i=M-1; i>=0; i--) {
              w[i+1] = w[i];                     /* data shift instruction */
              y += h[i] * w[i];                  /* MAC instruction */
              }

       return y;
}

/* \fn void gdelay2(int D, double c, double * q)
   \brief generalized circular delay with real-valued shift
   \param int D;
   \param double c, *q;                              \(c\)=shift, \(q\)=offset index 
   \return
*/
void gdelay2(int D, double c, double * q){
       *q -= c;                           /* decrement by \(c\) */

       if (*q < 0)  
              *q += D+1;

       if (*q > D)
              *q -= D+1;
}

/* \fn double gran(double m, double s, long * iseed)
   \brief gaussian random number generator
   \param double m, s;                              \(m\) = mean, \(s\sp{2}\) = variance 
   \param long *iseed;                              iseed passed by address 
   \return
*/
double gran(double m, double s, long * iseed)                 /* usage: x = gran(m, s, &iseed); */
{
       double v = 0;
       int i;

       for (i = 0; i < 12; i++)          /* sum 12 uniform random numbers */
              v += ran(iseed);

       return s * (v - 6) + m;           /* adjust mean and variance */
}   

/* \fn double I0(double x)
   \brief Modified Bessel Function I0(x)
   \param double x;
   \return
*/
double I0(double x)                              /* usage: y = I0(x) */
{
       int n = 1;
       double S = 1, D = 1, T;

       while (D > eps * S) {
              T = x / (2 * n++);
              D *= T * T;
              S += D;
              }

       return S;
}

/* \fn void ifft(int N, complex * X)
   \brief inverse FFT
   \param complex *X;
   \param int N;
   \return
*/
void ifft(int N, complex * X){
    int k;

    for (k=0; k<N; k++)
         X[k] = conjg(X[k]);                     /* conjugate input */

    fft(N, X);                                   /* compute FFT of conjugate */

    for (k=0; k<N; k++)
         X[k] = rdiv(conjg(X[k]), (double)N);    /* conjugate and divide by \(N\) */
}

/* \fn double lowpass(int D, double * w, double ** p, int M, double * a, double * b, double * v, double x)
   \brief lowpass reverberator with feedback filter G(z)
   \param double *w, **p, *a, *b, *v, x;                    \(v\) = state vector for \(G(z)\) 
   \param int D;                                            \(a,b,v\) are \((M+1)\)-dimensional 
   \return
*/
double lowpass(int D, double * w, double ** p, int M, double * a, double * b, double * v, double x){
       double y, sD;

       sD = tap(D, w, *p, D);                    /* delay output is \(G(z)\) input */
       y = x + can(M, a, M, b, v, sD);           /* reverb output */
       **p = y;                                  /* delay input */
       cdelay(D, w, p);                          /* update delay line */

       return y;
}

/* \fn void modwrap(int L, double * x, int N, double * xtilde)
   \brief modulo-N wrapping of length-L signal
   \param int L, N;                                     \(x\) is \(L\)-dimensional 
   \param double *x, *xtilde;                           xtilde is \(N\)-dimensional 
   \return
*/
void modwrap(int L, double * x, int N, double * xtilde)                /* usage: modwrap(L, x, N, xtilde); */
{
    int n, r, m, M;

    r = L % N;                               /* remainder \(r=0,1,\dotsc,N-1\) */
    M = (L-r) / N;                           /* quotient of division \(L/N\) */

    for (n=0; n<N; n++) {
         if (n < r)                          /* non-zero part of last block */
             xtilde[n] = x[M*N+n];           /* if \(L<N\), this is the only block */
         else
             xtilde[n] = 0;                  /* if \(L<N\), pad \(N-L\) zeros at end */

         for (m=M-1; m>=0; m--)              /* remaining blocks */
              xtilde[n] += x[m*N+n];         /* if \(L<N\), this loop is skipped */
         }
}

/* \fn double plain(int D, double * w, double ** p, double a, double x)
   \brief plain reverberator with circular delay line
   \param double *w, **p, a, x;                            \(p\) is passed by address 
   \param int D;
   \return
*/
double plain(int D, double * w, double ** p, double a, double x)                     /* usage: y=plain(D,w,&p,a,x); */
{
       double y, sD;

       sD = tap(D, w, *p, D);                   /* \(D\)th tap delay output */
       y = x + a * sD;                          /* filter output */
       **p = y;                                 /* delay input */
       cdelay(D, w, p);                         /* update delay line */

       return y;
}

/* \fn double ran(long * iseed)  
   \brief uniform random number generator in [0, 1)
   \param long *iseed;                                     iseed passed by address 
   \return
*/
double ran(long * iseed)                                /* usage: u = ran(&iseed); */
{
    *iseed = a_ * (*iseed % q_) - r_ * (*iseed / q_);          /* update seed */

    if (*iseed < 0)                              /* wrap to positive values */
           *iseed += m_;

    return (double) *iseed / (double) m_;
}  

/* \fn double ran1f(int B, double * u, int * q, long * iseed)
   \brief 1/f random number generator
   \param int B, *q;                                   \(q, u\) are \(B\)-dimensional 
   \param double *u;
   \param long *iseed;                                 passed by address 
   \return
*/
double ran1f(int B, double * u, int * q, long * iseed)                /* usage: y = ran1f(B, u, q, &iseed); */
{
    double y;
    int b;

    for(y=0, b=0; b<B; b++)
          y += ranh(1<<b, u+b, q+b, iseed);         /* period = (1<<b) = 2\(\sp{b}\) */

    return y / B;
}

/* \fn double ranh(int D, double * u, int * q, long * iseed)
   \brief hold random number generator of period D
   \param int D, *q;                                 \(q\) is cycled modulo-\(D\) 
   \param double *u;                                 \(u\) = 1-dimensional array 
   \param long *iseed;                               \(q\), iseed are passed by address 
   \return
*/
double ranh(int D, double * u, int * q, long * iseed)               /* usage: y = ranh(D, u, &q, &iseed); */
{
       double y;

       y = u[0];                              /* hold sample for \(D\) calls */

       cdelay2(D-1, q);                       /* decrement \(q\) and wrap mod-\(D\) */

       if (*q == 0)                           /* every \(D\) calls, */
              u[0] = ran(iseed) - 0.5;        /* get new \(u[0]\) (zero mean) */

       return y;
}

/* \fn double ranl(int D, double * u, int * q, long * iseed)
   \brief linearly interpolated random generator of period D
   \param int D, *q;                                 \(q\) is cycled modulo-\(D\) 
   \param double *u;                                 \(u\) = 2-dimensional array
   \param long *iseed;                               \(q\), iseed are passed by address 
   \return
*/
double ranl(int D, double * u, int * q, long * iseed)               /* usage: y = ranl(D, u, &q, &iseed); */
{
       double y;
       int i;

       i = (D - *q) % D;                      /* interpolation index */

       y = u[0] + (u[1] - u[0]) * i / D;      /* linear interpolation */

       cdelay2(D-1, q);                       /* decrement \(q\) and wrap mod-\(D\) */

       if (*q == 0) {                         /* every \(D\) calls, */
              u[0] = u[1];                    /* set new \(u[0]\) and */
              u[1] = ran(iseed) - 0.5;        /* get new \(u[1]\) (zero mean) */
              }

       return y;
}

/* \fn void shuffle(int N, complex * X)
   \brief in-place shuffling (bit-reversal) of a complex array
   \param complex *X;
   \param int N;                                     \(N\) must be a power of 2 
   \return
*/
void shuffle(int N, complex * X){
       int n, r, B=1;

       while ( (N >> B) > 0 )             /* \(B\) = number of bits */
              B++;

       B--;                               /* \(N = 2\sp{B}\) */

       for (n = 0; n < N; n++) {
           r = bitrev(n, B);              /* bit-reversed version of \(n\) */
           if (r < n) continue;           /* swap only half of the \(n\)s */
           swap(X+n, X+r);                /* swap by addresses */
           }
}

/* \fn double sine(int D, int i)
   \brief sine wavetable of length D 
   \param int D, i;
   \return
*/
double sine(int D, int i){
       double pi = 4 * atan(1.0);

       return sin(2 * pi * i / D);
}

/* \fn double sos(double * a, double * b, double * w, double x)
   \brief IIR filtering by single second order section
   \param double *a, *b, *w, x;                      \(a[0]=1\) always 
   \return
*/
double sos(double * a, double * b, double * w, double x)                    /* \(a, b, w\) are 3-dimensional */
{
       double y;

       w[0] = x - a[1] * w[1] - a[2] * w[2];
       y = b[0] * w[0] + b[1] * w[1] + b[2] * w[2];

       w[2] = w[1];
       w[1] = w[0];

       return y;
}

/* \fn double square(int D1, int i)
   \brief square wavetable of length D, with D1 ones
   \param int D1, i;
   \return
*/
double square(int D1, int i){
       if (i < D1)
              return 1;
       else
              return 0;
}

/* \fn void swap(complex * a, complex * b)
   \brief swap two complex numbers (by their addresses)
   \param complex *a, *b;
   \return
*/
void swap(complex * a, complex * b){
       complex t;

        t = *a;
       *a = *b;
       *b =  t;
}

/* \fn double tap(int D, double * w, double * p, int i)
   \brief i-th tap of circular delay-line buffer
   \param w, *p;                             \(p\) passed by value 
   \param int D, i;                           \(i=0,1,\dotsc, D\) 
   \return
*/
double tap(int D, double * w, double * p, int i)                    /* usage: si = tap(D, w, p, i); */
{
       return w[(p - w + i) % (D + 1)];
}

/* \fn double tap2(int D, double * w, int q, int i)
   \brief i-th tap of circular delay-line buffer
   \param double *w;
   \param int D, q, i;                               \(i=0,1,\dotsc,D\) 
   \return
*/
double tap2(int D, double * w, int q, int i)                   /* usage: si = tap2(D, w, q, i); */
{
       return w[(q + i) % (D + 1)];
}

/* \fn double tapi(int D, double * w, double * p, double d)
   \brief interpolated tap output of a delay line
   \param double *w, *p, d;                          \(d\) = desired non-integer delay
   \param int D;                                     \(p\) = circular pointer to \(w\) 
   \return
*/
double tapi(int D, double * w, double * p, double d)                   /* usage: sd = tapi(D, w, p, d); */
{
       int i, j;
       double si, sj;

       i = (int) d;                       /* interpolate between \(s\sb{i}\) and \(s\sb{j}\) */
       j = (i+1) % (D+1);                 /* if \(i=D\), then \(j=0\); otherwise, \(j=i+1\) */

       si = tap(D, w, p, i);              /* note, \(s\sb{i}(n) = x(n-i)\) */
       sj = tap(D, w, p, j);              /* note, \(s\sb{j}(n) = x(n-j)\) */

       return si + (d - i) * (sj - si);
}

/* \fn double tapi2(int D, double * w, int q, double d) 
   \brief interpolated tap output of a delay line
   \param double *w, d;                              \(d\) = desired non-integer delay
   \param int D, q;                                  \(q\) = circular offset index 
   \return
*/
double tapi2(int D, double * w, int q, double d)                  /* usage: sd = tapi2(D, w, q, d); */
{
       int i, j;
       double si, sj;

       i = (int) d;                       /* interpolate between \(s\sb{i}\) and \(s\sb{j}\) */
       j = (i+1) % (D+1);                 /* if \(i=D\), then \(j=0\); otherwise, \(j=i+1\) */

       si = tap2(D, w, q, i);             /* note, \(s\sb{i}(n) = x(n-i)\) */
       sj = tap2(D, w, q, j);             /* note, \(s\sb{j}(n) = x(n-j)\) */

       return si + (d - i) * (sj - si);
}

/* \fn double trapez(int D, int  D1, int  D2, int i)
   \brief trapezoidal wavetable: D1 rising, D2 steady
   \param int D, D1, D2, i;
   \return
*/
double trapez(int D, int  D1, int  D2, int i){
       if (i < D1)
              return i/(double) D1;
       else
            if (i < D1+D2)
                 return 1;
            else
                 return (D - i)/(double) (D - D1 - D2);
}

/* \fn int u(x)
   \brief unit step function
   \param double x;
   \return
*/
int u(double x){
       if (x >= 0)
              return 1;
       else
              return 0;
}

/* \fn double wavgen(int D, double * w, double A, double F, double * q)
   \brief wavetable generator (truncation method)
   \param int D;                               \(D\) = wavetable length 
   \param double *w, A, F, *q;                 \(A\) = amplitude, \(F\) = frequency, \(q\) = offset index 
   \return
*/
double wavgen(int D, double * w, double A, double F, double * q)        /* usage: y = wavgen(D, w, A, F, &q); */
{
       double y;
       int i;

       i = (int) (*q);                             /* truncate down */

       y = A * w[i];

       gdelay2(D-1, D*F, q);                       /* shift  \(c = DF\) */

       return y;
}


/* \fn double wavgeni(int D, double * w, double  A, double F, double * q)
   \brief wavetable generator (interpolation method)
   \param int D;                              \(D\) = wavetable length 
   \param double *w, A, F, *q;                \(A\) = amplitude, \(F\) = frequency, \(q\) = offset index 
   \return
*/
double wavgeni(int D, double * w, double  A, double F, double * q)      /* usage: y = wavgeni(D, w, A, F, &q); */
{
       double y;
       int i, j;

       i = (int) *q;                        /* interpolate between \(w[i], w[j]\) */
       j = (i + 1) % D;                     

       y = A * (w[i] + (*q - i) * (w[j] - w[i]));

       gdelay2(D-1, D*F, q);                     /* shift  \(c = DF\) */

       return y;
}

/* \fn double wavgenr(int D, double * w, double A, double F, double * q)
   \brief wavetable generator (rounding method)
   \param int D;                               \(D\) = wavetable length 
   \param double *w, A, F, *q;                 \(A\) = amplitude, \(F\) = frequency, \(q\) = offset index 
   \return
*/
double wavgenr(int D, double * w, double A, double F, double * q)       /* usage: y = wavgenr(D, w, A, F, &q); */
{
       double y;
       int k;

       k = (int) (*q + 0.5);                     /* round */

       y = A * w[k];

       gdelay2(D-1, D*F, q);                     /* shift  \(c = DF\) */

       return y;
}

/* \fn void wrap(int M, double * w, double ** p)
   \brief circular wrap of pointer p, relative to array w 
   \param double *w, **p;
   \param int M;
   \return
*/
void wrap(int M, double * w, double ** p){
       if (*p > w + M)  
              *p -= M + 1;          /* when \(*p=w+M+1\), it wraps around to \(*p=w\) */

       if (*p < w)  
              *p += M + 1;          /* when \(*p=w-1\), it wraps around to \(*p=w+M\) */
}

/* \fn void wrap2(int M, int * q)
   \brief circular wrap of pointer offset q, relative to array w
   \param int M, *q;
   \return
*/
void wrap2(int M, int * q){
       if (*q > M)  
              *q -= M + 1;          /* when \(*q=M+1\), it wraps around to \(*q=0\) */

       if (*q < 0)  
              *q += M + 1;          /* when \(*q=-1\), it wraps around to \(*q=M\) */
}

