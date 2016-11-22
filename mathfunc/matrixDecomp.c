#include "matrix.h"
#include "dsp.h"
#include "array.h"
#include <string.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

int LUPdecompose(double ** A, int size, int * P);
int LUPdecomposefloat(float ** A, int size, int * P);

int LUPinverse(double ** LU, double ** B, int size, int * P, double * X, double * Y);
int LUPinversefloat(float ** LU, float ** B, int size, int * P, float * X, float * Y);

void LUextraction(double ** LU, int n, double ** L, double ** U);
void LUextractionfloat(float ** LU, int n, float ** L, float ** U);

void LUinversion(int n, double ** m, double ** inv);
void LUinversionfloat(int n, float ** m, float ** inv);

void U_inverse(double ** U, int n);
void U_inversefloat(float ** U, int n);

void L_inverse(double ** L, int n);
void L_inversefloat(float ** L, int n);

void backwardsSub(double ** U, int n, double * x, double * b);
void backwardsSubfloat(float ** U, int n, float * x, float * b);

int dsvd(double **a, int m, int n, double *w, double **v);
double PYTHAG(double a, double b);
void svdInvf(double ** U, double * w, double ** V, int nrows, int ncols, double ** inv);

int LUPdecompose(double ** A, int size, int * P){
  int i, j, k, kd = 0, T;  
  double p, t;  
   
  /* Finding the pivot of the LUP decomposition. */  
  for(i = 0; i < size; i++) P[i] = i; //Initializing.  
   
    for(k = 0; k < size-1; k++){  
      p = 0;  
      for(i = k; i < size; i++){  
        t = A[i][k];  
        if(t < 0) t *= -1; //Abosolute value of 't'.  
        if(t > p){  
          p = t;  
          kd = i;  
        }  
      }  
   
      if(p == 0){  
        printf("\nLUPdecompose(): ERROR: A singular matrix is supplied.\n\tRefusing to proceed any further.\n");  
        return -1;  
      }  
   
      /* Exchanging the rows according to the pivot determined above. */  
      T = P[kd];  
      P[kd] = P[k];  
      P[k] = T;  
      for(i = 0; i < size; i++){  
        t = A[kd][i];  
        A[kd][i] = A[k][i];  
        A[k][i] = t;  
      }  
      
      for(i = k + 1; i < size; i++){ //Performing substraction to decompose A as LU.  
        A[i][k] = A[i][k]/(A[k][k] + EPS) + EPS; 
        for(j = k + 1; j < size; j++) A[i][j] -= A[i][k]*A[k][j];  
       }  
     } //Now, 'A' contains the L (without the diagonal elements, which are all 1)  
      //and the U.  
   
   return 0;  
}

int LUPdecomposefloat(float ** A, int size, int * P){
  int i, j, k, kd = 0;//, T;  
  float p, t;  
   
  /* Finding the pivot of the LUP decomposition. */  
  for(i = 0; i < size; i++) P[i] = i; //Initializing.  
   
    for(k = 0; k < size-1; k++){  
      p = 0;  
      for(i = k; i < size; i++){  
        t = A[i][k];  
        if(t < 0) t *= -1; //Abosolute value of 't'.  
        if(t > p){  
          p = t;  
          kd = i;  
        }  
      }  
   
      if(p == 0){  
        printf("\nLUPdecompose(): ERROR: A singular matrix is supplied.\n\tRefusing to proceed any further.\n");  
        return -1;  
      }  
   
      // /* Exchanging the rows according to the pivot determined above. */  
      // T = P[kd];  
      // P[kd] = P[k];  
      // P[k] = T;  
      // for(i = 0; i < size; i++){  
      //   t = A[kd][i];  
      //   A[kd][i] = A[k][i];  
      //   A[k][i] = t;  
      // }  
      
      for(i = k + 1; i < size; i++){ //Performing substraction to decompose A as LU.  
        A[i][k] = A[i][k]/(A[k][k] + EPS) + EPS; 
        for(j = k + 1; j < size; j++) A[i][j] -= A[i][k]*A[k][j];  
       }  
     } //Now, 'A' contains the L (without the diagonal elements, which are all 1)  
      //and the U.  
   
   return 0;  
}


int LUPinverse(double ** LU, double ** B, int size, int * P, double * X, double * Y){

  int i, j, n, m;  
  double t;  
   
  //Initializing X and Y.  
  for(n = 0; n < size; n++) X[n] = Y[n] = EPS;  
   
  /* Solving LUX = Pe, in order to calculate the inverse of 'A'. Here, 'e' is a column  
  * vector of the identity matrix of size 'size-1'. Solving for all 'e'. */  
  for(i = 0; i < size; i++) {  
    //Storing elements of the i-th column of the identity matrix in i-th row of 'B'.  
    for(j = 0; j < size; j++) B[i][j] = 0;  
      B[i][i] = 1;  
   
    //Solving Ly = Pb.  
    for(n = 0; n < size; n++){  
      t = 0;  
      for(m = 0; m < n-1; m++) t += LU[n][m]*Y[m];  
      Y[n] = B[i][P[n]]-t;  
    }  
   
   //Solving Ux = y.  
    for(n = size - 1; n >= 0; n--){  
      t = 0;  
      for(m = n+1; m < size; m++) t += LU[n][m]*X[m];  
      X[n] = (Y[n]-t)/LU[n][n];  
     }//Now, X contains the solution.
     //printArrayDouble("X", size, X);
   
    for(j = 0; j < size; j++) B[i][j] = X[j]; //Copying 'X' into the same row of 'B'.
  } //Now, 'B' the transpose of the inverse of 'A'.  
   
  /* Copying transpose of 'B' into 'LU', which would the inverse of 'A'. */  
  for(i = 0; i < size; i++) for(j = 0; j < size; j++) LU[i][j] = B[j][i];  
   
  return 0;  

}

int LUPinversefloat(float ** LU, float ** B, int size, int * P, float * X, float * Y){

  int i, j, n, m;  
  float t;  
   
  //Initializing X and Y.  
  for(n = 0; n < size; n++) X[n] = Y[n] = EPS;  
   
  /* Solving LUX = Pe, in order to calculate the inverse of 'A'. Here, 'e' is a column  
  * vector of the identity matrix of size 'size-1'. Solving for all 'e'. */  
  for(i = 0; i < size; i++) {  
    //Storing elements of the i-th column of the identity matrix in i-th row of 'B'.  
    for(j = 0; j < size; j++) B[i][j] = 0;  
      B[i][i] = 1;  
   
    //Solving Ly = Pb.  
    for(n = 0; n < size; n++){  
      t = 0;  
      for(m = 0; m < n-1; m++) t += LU[n][m]*Y[m];  
      Y[n] = B[i][P[n]]-t;  
    }  
   
   //Solving Ux = y.  
    for(n = size - 1; n >= 0; n--){  
      t = 0;  
      for(m = n+1; m < size; m++) t += LU[n][m]*X[m];  
      X[n] = (Y[n]-t)/LU[n][n];  
     }//Now, X contains the solution.
     //printArrayDouble("X", size, X);
   
    for(j = 0; j < size; j++) B[i][j] = X[j]; //Copying 'X' into the same row of 'B'.
  } //Now, 'B' the transpose of the inverse of 'A'.  
   
  /* Copying transpose of 'B' into 'LU', which would the inverse of 'A'. */  
  for(i = 0; i < size; i++) for(j = 0; j < size; j++) LU[i][j] = B[j][i];  
   
  return 0;  

}

void LUinversion(int n, double ** m, double ** inv){
	int i, j;
  int * P = calloc(n, sizeof(int));

	double ** temp_m = calloc(n, sizeof(double *));
	for(i = 0; i < n; i++){
  	temp_m[i] = calloc(n, sizeof(double));
  	for(j = 0; j < n; j++)
    		temp_m[i][j] = m[i][j];
	}

	LUPdecompose(m, n, P);

	double ** L = calloc(n, sizeof(double *));
	for(i = 0; i < n; i++)
  	L[i] = calloc(n, sizeof(double));

	double ** U = calloc(n, sizeof(double *));
	for(i = 0; i < n; i++)
  	U[i] = calloc(n, sizeof(double));

  LUextraction(m, n, L, U);

	U_inverse(U, n);
	L_inverse(L, n);

  multiplyMatrixDouble(n, n, U, n, n, L, inv);

	for(i = 0; i < n; i++)
  	for(j = 0; j < n; j++)
    		m[i][j] = temp_m[i][j];

  free(temp_m);
  free(L);
  free(U);
  free(P);
}

void LUinversionfloat(int n, float ** m, float ** inv){
  int i, j;
  int * P = calloc(n, sizeof(int));

  float ** temp_m = calloc(n, sizeof(float *));
  for(i = 0; i < n; i++){
    temp_m[i] = calloc(n, sizeof(float));
    for(j = 0; j < n; j++)
        temp_m[i][j] = m[i][j];
  }

  LUPdecomposefloat(m, n, P);

  float ** L = calloc(n, sizeof(float *));
  for(i = 0; i < n; i++)
    L[i] = calloc(n, sizeof(float));

  float ** U = calloc(n, sizeof(float *));
  for(i = 0; i < n; i++)
    U[i] = calloc(n, sizeof(float));

  LUextractionfloat(m, n, L, U);

  U_inversefloat(U, n);
  L_inversefloat(L, n);

  multiplyMatrixFloat(n, n, U, n, n, L, inv);

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
        m[i][j] = temp_m[i][j];

  free(temp_m);
  free(L);
  free(U);
  free(P);
}

void LUextraction(double ** LU, int n, double ** L, double ** U){
	int i, j;

	for(i = 0; i < n; i++){
		L[i][i] = 1;
		for(j = i; j < n; j++)
			U[i][j] = LU[i][j];
	}

	for(i = n - 1; i >= 0; i--)
   		for(j = i - 1; j >= 0; j--)
   			L[i][j] = LU[i][j]; 

}

void LUextractionfloat(float ** LU, int n, float ** L, float ** U){
  int i, j;

  for(i = 0; i < n; i++){
    L[i][i] = 1;
    for(j = i; j < n; j++)
      U[i][j] = LU[i][j];
  }

  for(i = n - 1; i >= 0; i--)
      for(j = i - 1; j >= 0; j--)
        L[i][j] = LU[i][j]; 

}

void L_inverse(double ** L, int n){
   int i, j, k;
//         Invert the subdiagonal part of the matrix L row by row where
//         the diagonal elements are assumed to be 1.0.

   for(i = 1; i < n; i++){
   		for(j = 0; j < i; j++){
   			L[i][j] = -L[i][j];	
   			for(k = j + 1; k < i; k++)
   				L[i][j] -= L[i][k] * L[k][j];
   		}
   }
}

void L_inversefloat(float ** L, int n){
   int i, j, k;
//         Invert the subdiagonal part of the matrix L row by row where
//         the diagonal elements are assumed to be 1.0.

   for(i = 1; i < n; i++){
      for(j = 0; j < i; j++){
        L[i][j] = -L[i][j]; 
        for(k = j + 1; k < i; k++)
          L[i][j] -= L[i][k] * L[k][j];
      }
   }
}

void U_inverse(double ** U, int n){
  int i, j;  

  for(i = n - 1; i >= 0; i--){

  		double * x = calloc(i + 1, sizeof(double));
  		double * b = calloc(i + 1, sizeof(double));
  		
  		b[i] = 1.0;
  		backwardsSub(U, i + 1, x, b);

  		for(j = 0; j < i + 1; j++){

  			U[j][i] = x[j];
  		}

  		free(x);
  		free(b);

  }
   
}

void U_inversefloat(float ** U, int n){
  int i, j;  

  for(i = n - 1; i >= 0; i--){

      float * x = calloc(i + 1, sizeof(float));
      float * b = calloc(i + 1, sizeof(float));
      
      b[i] = 1.0;
      backwardsSubfloat(U, i + 1, x, b);

      for(j = 0; j < i + 1; j++){

        U[j][i] = x[j];
      }

      free(x);
      free(b);

  }
   
}

void backwardsSub(double ** U, int n, double * x, double * b){
	int i, j;

	for(i = n - 1; i >= 0; i--){
		
		x[i] = b[i];
		for(j = i + 1; j < n; j++)
			x[i] -= U[i][j]*x[j];

		x[i] /= U[i][i];

	}

}

void backwardsSubfloat(float ** U, int n, float * x, float * b){
  int i, j;

  for(i = n - 1; i >= 0; i--){
    
    x[i] = b[i];
    for(j = i + 1; j < n; j++)
      x[i] -= U[i][j]*x[j];

    x[i] /= U[i][i];

  }

}

double PYTHAG(double a, double b)
{
    double at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}


int dsvd(double **a, int m, int n, double *w, double **v)
{
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;
  
    if (m < n){
        fprintf(stderr, "#rows must be > #cols \n");
        return(0);
    }
  
    rv1 = (double *)malloc((unsigned int) n*sizeof(double));

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++){
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m){
            for (k = i; k < m; k++) 
                scale += fabs(a[k][i]);
            if (scale){
                for (k = i; k < m; k++){
                    a[k][i] /= scale;
                    s += a[k][i] * a[k][i];
                }
                f = a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = f - g;
                // if (i != n - 1){
                    for (j = l; j < n; j++){
                        for (s = 0.0, k = i; k < m; k++) 
                            s += (a[k][i] * a[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++) 
                            a[k][j] += (f * a[k][i]);
                    }
                // }
                for (k = i; k < m; k++) 
                    a[k][i] = (a[k][i]*scale);
            }
        }
        w[i] = (scale * g);
    
        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1){
            for (k = l; k < n; k++) 
                scale += fabs(a[i][k]);
            if (scale){
                for (k = l; k < n; k++){
                    a[i][k] = (a[i][k]/scale);
                    s += (a[i][k] * a[i][k]);
                }
                f = a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = (f - g);
                for (k = l; k < n; k++) 
                    rv1[k] = a[i][k] / h;
                if (i != m - 1){
                    for (j = l; j < m; j++){
                        for (s = 0.0, k = l; k < n; k++) 
                            s += (a[j][k] * a[i][k]);
                        for (k = l; k < n; k++) 
                            a[j][k] += (s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++) 
                    a[i][k] = (a[i][k]*scale);
            }
        }
        anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }
  
    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--){
        if (i < n - 1){
            if (g){
                for (j = l; j < n; j++)
                    v[j][i] = ((a[i][j] / a[i][l]) / g);
                    /* double division to avoid underflow */
                for (j = l; j < n; j++){
                    for (s = 0.0, k = l; k < n; k++) 
                        s += (a[i][k] * v[k][j]);
                    for (k = l; k < n; k++) 
                        v[k][j] += (s * v[k][i]);
                }
            }
            for (j = l; j < n; j++) 
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
  
    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--){
        l = i + 1;
        g = w[i];
        if (i < n - 1) 
            for (j = l; j < n; j++) 
                a[i][j] = 0.0;
        if (g){
            g = 1.0 / g;
            if (i != n - 1){
                for (j = l; j < n; j++){
                    for (s = 0.0, k = l; k < m; k++) 
                        s += (a[k][i] * a[k][j]);
                    f = (s / a[i][i]) * g;
                    for (k = i; k < m; k++) 
                        a[k][j] += (f * a[k][i]);
                }
            }
            for (j = i; j < m; j++) 
                a[j][i] = (a[j][i]*g);
        }
        else 
        {
            for (j = i; j < m; j++) 
                a[j][i] = 0.0;
        }
        ++a[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--){ /* loop over singular values */
        for (its = 0; its < 30; its++){ /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--){ /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm){
                    flag = 0;
                    break;
                }
                if (fabs(w[nm]) + anorm == anorm) 
                    break;
            }
            if (flag) {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm) {
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h; 
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++) {
                            y = a[j][nm];
                            z = a[j][i];
                            a[j][nm] = (y * c + z * s);
                            a[j][i] = (z * c - y * s);
                        }
                    }
                }
            }
            z = w[k];
            if (l == k) { /* convergence */
                if (z < 0.0) { /* make singular value nonnegative */
                    w[k] = (-z);
                    for (j = 0; j < n; j++) 
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                free((void*) rv1);
                fprintf(stderr, "No convergence after 30,000! iterations \n");
                return(0);
            }
    
            /* shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
          
            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++) {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++) {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = (x * c + z * s);
                    v[jj][i] = (z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = z;
                if (z) {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++) {
                    y = a[jj][j];
                    z = a[jj][i];
                    a[jj][j] = (y * c + z * s);
                    a[jj][i] = (z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }
    free((void*) rv1);
    return(1);
}

void svdInvf(double ** U, double * w, double ** V, int nrows, int ncols, double ** inv)
{
  // Call this after svdcmp() to compute the inverse
  int i, j, k;
  double * wInv;
  double wmin, wmax;
  
  // generate Ut
  //Ut = (double *) malloc(ncols * nrows * sizeof(double));
  double ** Ut = malloc(ncols*sizeof(double *));
  for(i = 0; i < ncols; i++)
    Ut[i] = malloc(nrows*sizeof(double));

  for (i=0; i<ncols; ++i){
    for (j=0; j<nrows; ++j){
      //*(Ut + i*nrows + j) = *(U + j*ncols + i);
      Ut[i][j] = U[j][i];
    }
  }

  // compute 1/w[i], but set small w[i] values to zero instead
  wmax = 0.;
  for (i=0; i<ncols; ++i){
    if (fabs(w[i]) > wmax)
      wmax = fabs(w[i]);
  }

  // wmin = wmax * 1e-6;
  wInv = (double *) malloc(ncols * sizeof(double));
  for (i=0; i<ncols; ++i){
    if (fabs(w[i]) < wmin)
      wInv[i] = 0.;
    else
      wInv[i] = 1. / w[i];
  }

  // compute invM
  for (i=0; i<ncols; ++i){
    for (j=0; j<nrows; ++j){
  //     //*(inv + i*nrows + j) = 0.;
      inv[i][j] = 0.;
      for (k=0; k<ncols; ++k){
  //       //*(inv + i*nrows + j) += *(V + i*ncols + k) * (wInv[k] * *(Ut + k*nrows + j));
        inv[i][j] = inv[i][j] + V[i][k] * wInv[k] * Ut[k][j];

      }
    }
  }

  free(wInv);
  free(Ut);

}  // end of svdInv()
