/////////////////////////////////////////////////////
//
//                              singValD.c
//
//  Copyright (c) 2003
//      John N. Sanders-Reed
//      26 Meadow View Rd.
//      Cedar Crest, NM 87008
//      505-281-8563
//  All Rights Reserved.
//
//      This software is intended for single machine use.
//      Archival copies are permitted but other reproduction or distribution
//      of this software in any form is strictly prohibited. For complete
//      details, consult the License agreement.
//
//  Singular Value Decomposition: A = U w Vt
//  Matrix A has m rows and n columns.
//  On output, input A is replaced by U of same dimesnions
//  On output, V is an n x n matrix (Note that V, not Vt is returned)
//  w is an n x n strictly diagonal matrix, so only the n diagonal elements are returned
//
//  Other properites:
//     U * Ut = V * Vt = 1 (identity matrix
//     inv(A) = V 1/wi Ut
//       For m != n, this is strictly an inv(A) * A = 1 inverse. The other order does NOT work
//
//  Least Squares:
//   Solve A x = b for unknown x, given A and b
//   Since inv(A) = V 1/wi Ut, we have
//    x = V 1/wi Ut b
//   For ill-conditioned or singular matrices, set small wi to zero.
//     Find wimax, compute wimin = wimax * 1e-6, set any smaller wi to zero
//       This sets any wi which is more than 1e6 smaller than the largest wi to zero
//
//  Reference: Numerical Recipes in C
//
/////////////////////////////////////////////////////
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "mathpack.h"

///////////////////////////////////////////////////////////
#define SIGN(a,b) ((b) >= 0. ? fabs(a) : -fabs(a))

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

double pythag(double a, double b);

/////////////////////////////////////////////////////
//void svdcmp(double **a, int m, int n, double w[], double **v)
int singValDecomp(double *aExt, int m, int n, double * w, double *vExt)
{
	// Given a matrix a[1 .. m][1 .. n], this routine computes its singular
	// value decomposition, A = U W Vt. The matrix U replaces a on output.
	//m=rows, n=columns
	// The diagonal matrix of singular values W is output as a vector w[1 .. n].
	// The matrix V (not the transpose Vt) is output as v[1 .. n][1 .. n].

	int flag, i, its, j, jj, k, L, nm;
	double anorm, c, f, g, h, s, scale, x, y, z, *rv1;

	////////////////////////////////////
	double **a, **v;

	a = (double **) calloc(m, sizeof(double *));
	for (i=0; i<m; ++i)
		a[i] = aExt + i*n - 1;
	--a;

	v = (double **) calloc(n, sizeof(double *));
	for (i=0; i<n; ++i)
		v[i] = vExt + i*n - 1;
	--v;

	w--;

	////////////////////////////

	rv1 = (double *) malloc(n * sizeof(double));
	--rv1;
	g = scale = anorm = 0.;  // Householder reduction to bidiagonal form
	for (i=1; i<=n; ++i)
	{
		L = i + 1;
		rv1[i] = scale * g;
		g = s = scale = 0.;
		if (i <= m)
		{
			for (k=i; k<=m; ++k)
				scale += fabs(a[k][i]);
			if (scale)
			{
				for (k=i; k<=m; ++k)
				{
					a[k][i] /= scale;
					s += a[k][i] * a[k][i];
				}
				f = a[i][i];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][i] = f - g;
				for (j=L; j<=n; ++j)
				{
					for (s=0., k=i; k<=m; ++k)
						s += a[k][i] * a[k][j];
					f = s / h;
					for (k=i; k<=m; ++k)
						a[k][j] += f * a[k][i];
				}  // end of for loop over j
				for (k=i; k<=m; ++k)
					a[k][i] *= scale;
			}  // end of if (scale)
		}  // end of if (i <= m)

		w[i] = scale * g;
		g = s = scale = 0.;
		if ((i <= m) && (i != n))
		{
			for (k=L; k<=n; ++k)
				scale += fabs(a[i][k]);
			if (scale)
			{
				for (k=L; k<=n; ++k)
				{
					a[i][k] /= scale;
					s += a[i][k] * a[i][k];
				}
				f = a[i][L];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][L] = f - g;
				for (k=L; k<=n; ++k)
					rv1[k] = a[i][k] / h;
				for (j=L; j<=m; ++j)
				{
					for (s=0.,k=L; k<=n; ++k)
						s += a[j][k] * a[i][k];
					for (k=L; k<=n; ++k)
						a[j][k] += s * rv1[k];
				}
				for (k=L; k<=n; ++k)
					a[i][k] *= scale;

			}  // end of if (scale)
		}  // end of if ((i <= m) && (i != n))
		anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
	}  // end of i-loop over n, incrementing

	// Accumulation of right-hand transformations
	for (i=n; i>=1; --i)
	{
		if (i < n)
		{
			if (g)
			{
				// Double division to avoid possible underflow
				for (j=L; j<=n; ++j)
					v[j][i] = (a[i][j] / a[i][L]) / g;
				for (j=L; j<=n; ++j)
				{
					for (s=0.,k=L; k<=n; ++k)
						s += a[i][k] * v[k][j];
					for (k=L; k<=n; ++k)
						v[k][j] += s * v[k][i];
				}
			}
			for (j=L; j<=n; ++j)
				v[i][j] = v[j][i] = 0.;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		L = i;
	}  // end of i-loop over n, decrementing

	// Accumulation of left-hand transformations
	for (i=MIN(m,n); i>=1; --i)
	{
		L = i + 1;
		g = w[i];
		for (j=L; j<=n; ++j)
			a[i][j] = 0.;
		if (g)
		{
			g = 1.0 / g;
			for (j=L; j<=n; ++j)
			{
				for (s=0.,k=L; k<=m; ++k)
					s += a[k][i] * a[k][j];
				f = (s / a[i][i]) * g;
				for (k=i; k<=m; ++k)
					a[k][j] += f * a[k][i];
			}
			for (j=i; j<=m; ++j)
				a[j][i] *= g;
		}

		else
		{
			for (j=i; j<=m; ++j)
				a[j][i] = 0.;
		}
		++a[i][i];
	}  // end of i-loop from MIN

	// Diagonalization of the bidiagonal form: Loop over
	//   singular values, and over allowed iterations
	for (k=n; k>=1; --k)
	{
		for (its=1; its<=30; ++its)
		{
			flag = 1;
			// Test for splitting. Note that rv1[1] is always zero
			for (L=k; L>=1; --L)
			{
				nm = L - 1;
				if ((fabs(rv1[L]) + anorm) == anorm)
				{
					flag = 0;
					break;
				}
				if ((fabs(w[nm]) + anorm) == anorm)
				{
					break;
				}
			}  // end of l-loop over k

			// Cancellation of rv1[L], if L > 1
			if (flag)
			{
				c = 0.;
				s = 1.0;
				for (i=L; i<=k; ++i)
				{
					f = s * rv1[i];
					rv1[i] = c * rv1[i];
					if ((fabs(f) + anorm) == anorm)
					{
						break;
					}
					g = w[i];
					h = pythag(f, g);
					w[i] = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;
					for (j=1; j<=m; ++j)
					{
						y = a[j][nm];
						z = a[j][i];
						a[j][nm] = y * c + z * s;
						a[j][i] = z * c - y * s;
					}  // end of j-loop over m
				}  // end of i-loop
			}  // end of if (flag)

			z = w[k];
			if (L == k)  // Convergence
			{
				if (z < 0.)  // Singular value is made nonnegative
				{
					w[k] = -z;
					for (j=1; j<=n; ++j)
						v[j][k] = -v[j][k];
				}
				break;
			}

			if (its == 30)
			{
				++rv1;
				free(rv1);
				++a;
				++v;
				free(a);
				free(v);
				return(1);  // failed
			}

			// Shift from bottom 2-by-2 minor.
			x = w[L];
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y-z) * (y+z) + (g-h) * (g+h)) / (2. * h * y);
			g = pythag(f, 1.0);
			f = ((x-z) * (x+z) + h * ((y / (f + SIGN(g,f))) - h)) / x;
			c = s = 1.0;

			// Next QR tranformation
			for (j=L; j<=nm; ++j)
			{
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x*c + g*s;
				g = g*c - x*s;
				h = y * s;
				y *= c;

				for (jj=1; jj<=n; ++jj)
				{
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x*c + z*s;
					v[jj][i] = z*c - x*s;
				}

				z = pythag(f, h);
				w[j] = z;

				// Rotation can be arbitrary if z = 0
				if (z)
				{
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}

				f = c*g + s*y;
				x = c*y - s*g;

				for (jj=1; jj<=m; ++jj)
				{
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y*c + z*s;
					a[jj][i] = z*c - y*s;
				} // end of jj-loop
			}  // end of j-loop

			rv1[L] = 0.;
			rv1[k] = f;
			w[k] = x;

		}  // end of its loop
	}  // end of k-loop

	++rv1;
	free(rv1);

	// JSR added
	++a;
	++v;
	free(a);
	free(v);

	return(0);  // success
}  // end of svdcmp()

/////////////////////////////////////////////////////
double pythag(double a, double b)
{
	// computes sqrt(a*a + b*b) without destructive overflow or underflow
	double absa, absb;
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb)
		return(absa * sqrt(1. + sqr(absb / absa)));
	else
		return(absb == 0. ? 0. : absb * sqrt(1. + sqr(absa / absb)));
}  // end of pythag()

/////////////////////////////////////////////////////
void svdInv(double * U, double * w, double * V, int nrows, int ncols, double * inv)
{
	// Call this after svdcmp() to compute the inverse
	int i, j, k;
	double * Ut, * wInv;
	double wmin, wmax;
	
	// generate Ut
	Ut = (double *) malloc(ncols * nrows * sizeof(double));

	for (i=0; i<ncols; ++i)
	{
		for (j=0; j<nrows; ++j)
		{
			*(Ut + i*nrows + j) = *(U + j*ncols + i);
		}
	}

	// compute 1/w[i], but set small w[i] values to zero instead
	wmax = 0.;
	for (i=0; i<ncols; ++i)
	{
		if (fabs(w[i]) > wmax)
			wmax = fabs(w[i]);
	}

	wmin = wmax * 1e-6;
	wInv = (double *) malloc(ncols * sizeof(double));
	for (i=0; i<ncols; ++i)
	{
		if (fabs(w[i]) < wmin)
			wInv[i] = 0.;
		else
			wInv[i] = 1. / w[i];
	}

	// compute invM
	for (i=0; i<ncols; ++i)
	{
		for (j=0; j<nrows; ++j)
		{
			*(inv + i*nrows + j) = 0.;
			for (k=0; k<ncols; ++k)
			{
				*(inv + i*nrows + j) += *(V + i*ncols + k) * (wInv[k] * *(Ut + k*nrows + j));
			}
		}
	}

	free(wInv);
	free(Ut);

}  // end of svdInv()

/////////////////////////////////////////////////////
