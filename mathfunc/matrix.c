
#include "matrix.h"
#include "dsp.h"
#include "array.h"
#include <string.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

void printMatrixShort(char * name, int * dim, short ** arr);
void printMatrixUShort(char * name, int * dim, unsigned short ** arr);
void printMatrixInt(char * name, int * dim, int ** arr);
void printMatrixUInt(char * name, int * dim, unsigned int ** arr);
void printMatrixFloat(char * name, int * dim, float ** arr);
void printMatrixDouble(char * name, int * dim, double ** arr);

void printMatrix3DShort(char * name, int * dim, int plane, short *** arr);
void printMatrix3DUShort(char * name, int * dim, int plane, unsigned short *** arr);
void printMatrix3DInt(char * name, int * dim, int plane, int *** arr);
void printMatrix3DUInt(char * name, int * dim, int plane, unsigned int *** arr);
void printMatrix3DFloat(char * name, int * dim, int plane, float *** arr);
void printMatrix3DDouble(char * name, int * dim, int plane, double *** arr);

void setMatrixFloat(int r, int c, float val, float ** m);
void setMatrixDouble(int r, int c, double val, double ** m);
void setMatrixInt(int r, int c, int val, int ** m);
void setMatrixUInt(int r, int c, unsigned int val, unsigned int ** m);

void setMatrix3DFloat(float *** m, int r, int c, int z, float val);
void setMatrix3DDouble(double *** m, int r, int c, int z, double val);
void setMatrix3DInt(int *** m, int r, int c, int z, int val);
void setMatrix3DUInt(unsigned int *** m, int r, int c,  int z, unsigned int val);

void getMatrixColumnUInt(unsigned int ** mat, int r, int c, int col, unsigned int * column);
void getMatrixColumnInt(int ** mat, int r, int c, int col, int * column);
void getMatrixColumnFloat(float ** mat, int r, int c, int col, float * column);
void getMatrixColumnDouble(double ** mat, int r, int c, int col, double * column);

void getMatrixRowUInt(unsigned int ** mat, int r, int c, int rw, unsigned int * row);
void getMatrixRowInt(int ** mat, int r, int c, int rw, int * row);
void getMatrixRowFloat(float ** mat, int r, int c, int rw, float * row);
void getMatrixRowDouble(double ** mat, int r, int c, int rw, double * row);

void multiplyMatrixFloat(int r1, int c1, float ** m1, int r2, int c2, float ** m2, float ** prd);
void multiplyMatrixDouble(int r1, int c1, double ** m1, int r2, int c2, double ** m2, double ** prd);
void multiplyMatrixInt(int r1, int c1, int ** m1, int r2, int c2, int ** m2, int ** prd);
void multiplyMatrixUInt(int r1, int c1, unsigned int ** m1, int r2, int c2, unsigned int ** m2, unsigned int ** prd);

void multiplyMatrixVectorFloat(int r, int c, float ** m, float * v, int l,  float * prd);
void multiplyMatrixVectorDouble(int r, int c, double ** m, double * v, int l,  double * prd);
void multiplyMatrixVectorInt(int r, int c, int ** m, int * v, int l, int * prd);
void multiplyMatrixVectorUInt(int r, int c, unsigned int ** m, unsigned int * v, int l, unsigned int * prd);

void transposeFloat(float ** m, int r, int c, float ** trans);
void transposeDouble(double ** m, int r, int c, double ** trans);
void transposeInt(int ** m, int r, int c, int ** trans);
void transposeUInt(unsigned int ** m, int r, int c, unsigned int ** trans);

float determinantf(float ** m, int n);
double determinantd(double ** m, int n);

void pinvf(int nr, int nc, float ** M, float ** Minv);
void pinvd(int nr, int nc, double ** M, double ** Minv);

void inversef(int n, float ** m, float ** inv);

void mean_2Df(float ** m, int h, int w, float * r);
void mean_2Dd(double ** m, int h, int w, double * r);

void interp2(unsigned int ** input, int rows, int cols, double * x, double * y, double sample_size, double * output);
double bilininterp(unsigned int ** input, int rows, int cols, double x, double y);

void convert2D21Dfloat(int r, int c, float ** mat2D, float * mat1D);
void convert1D22Dfloat(int r, int c, float * mat1D, float ** mat2D);

void convert2D21Ddouble(int r, int c, double ** mat2D, double * mat1D);
void convert1D22Ddouble(int r, int c, double * mat1D, double ** mat2D);

void convert2D21Dint(int r, int c, int (* mat2D)[c], int * mat1D);
void convert1D22Dint(int r, int c, int * mat1D, int (* mat2D)[c]);

void convert2D21Duint(int r, int c, unsigned int ** mat2D, unsigned int * mat1D);
void convert1D22Duint(int r, int c, unsigned int * mat1D, unsigned int ** mat2D);

void convert2D21Dushort(int r, int c, unsigned short ** mat2D, unsigned short * mat1D);
void convert1D22Dushort(int r, int c, unsigned short * mat1D, unsigned short ** mat2D);

void printmatrixfloat(int r, int c, float (* mat)[c]);
void printmatrixdouble(int r, int c, double (* mat)[c]);

double min2Dd(int r, int c, double ** mat);
double max2Dd(int r, int c, double ** mat);

double sum2Dd(int r, int c, double ** mat);
double std2Dd(int r, int c, double ** mat);

float sum2Df(int r, int c, float ** mat);
float std2Df(int r, int c, float ** mat);
float nonzerostd2Df(int r, int c, float ** mat);

double nonzerostd2Dd(int r, int c, double ** mat);



void getPlane(double *** mat3d, int r, int c, int plane, double ** matout);
void subtractMatrixDouble(int r, int c, double ** mat1, double ** mat2, double ** matout);
void subtractMatrixFloat(int r, int c, float ** mat1, float ** mat2, float ** matout);

void getRowUInt(int irow, unsigned int ** mat, int rows, int cols, unsigned int * row);
void getRowInt(int irow, int ** mat, int rows, int cols, int * row);
void getRowFloat(int irow, float ** mat, int rows, int cols, float * row);
void getRowDouble(int irow, double ** mat, int rows, int cols, double * row);


void nonzero_elements_matrix_uint(int r, int c, unsigned int ** mat, unsigned int * elements, int * size);
void nonzero_elements_matrix_int(int r, int c, int ** mat, int * elements, int * size);
void nonzero_elements_matrix_float(int r, int c, float ** mat, float * elements, int * size);
void nonzero_elements_matrix_double(int r, int c, double ** mat, double * elements, int * size);

void fliplruint(int r, int c, unsigned int ** mat);
void fliplrint(int r, int c, int ** mat);
void fliplrfloat(int r, int c, float ** mat);
void fliplrdouble(int r, int c, double ** mat);

void printMatrixShort(char * name, int * dim, short ** arr){
  int i, j;
  printf("%s\n", name);
  for(i = 0; i < dim[0]; i++){
    printf("[ ");
    for(j = 0; j < dim[1]; j++)
      printf("%d ", arr[i][j]);
    printf("]\n");
  }
}

void printMatrixUShort(char * name, int * dim, unsigned short ** arr){
  int i, j;
  printf("%s\n", name);
  for(i = 0; i < dim[0]; i++){
    printf("[ ");
    for(j = 0; j < dim[1]; j++)
      printf("%d ", arr[i][j]);
    printf("]\n");
  }
}

void printMatrixInt(char * name, int * dim, int ** arr){
  int i, j;
  printf("%s\n", name);
  for(i = 0; i < dim[0]; i++){
    printf("[ ");
    for(j = 0; j < dim[1]; j++)
      printf("%d ", arr[i][j]);
    printf("]\n");
  }
}

void printMatrixUInt(char * name, int * dim, unsigned int ** arr){
  int i, j;
  printf("%s\n", name);
  for(i = 0; i < dim[0]; i++){
    printf("[ ");
    for(j = 0; j < dim[1]; j++)
      printf("%d ", arr[i][j]);
    printf("]\n");
  }
}

void printMatrixFloat(char * name, int * dim, float ** arr){
  int i, j;
  int r, c;
  r = dim[0];
  c = dim[1];
  printf("%s[%d][%d]: \n", name, r, c);
  for(i = 0; i < dim[0]; i++){
    printf("%d: [ ", i);
    for(j = 0; j < dim[1]; j++)
      printf("%f ", arr[i][j]);
    printf("]\n");
  }
}

void printMatrixDouble(char * name, int * dim, double ** arr){
  int i, j;
  int r, c;
  r = dim[0];
  c = dim[1];
  printf("%s[%d][%d]: \n", name, r, c);
  for(i = 0; i < dim[0]; i++){
    printf("%d: [ ", i);
    for(j = 0; j < dim[1]; j++)
      printf("%e ", arr[i][j]);
    printf("]\n");
  }
}

void printMatrix3DShort(char * name, int * dim, int plane, short *** arr){
  int i, j;
  int r, c;
  r = dim[0];
  c = dim[1];
  printf("%s[%d][%d][%d]: \n", name, r, c, plane);
  for(i = 0; i < dim[0]; i++){
    printf("[ ");
    for(j = 0; j < dim[1]; j++)
      printf("%d ", arr[i][j][plane]);
    printf("]\n");
  }
}

void printMatrix3DUShort(char * name, int * dim, int plane, unsigned short *** arr){
  int i, j;
  int r, c;
  r = dim[0];
  c = dim[1];
  printf("%s[%d][%d][%d]: \n", name, r, c, plane);
  for(i = 0; i < dim[0]; i++){
    printf("[ ");
    for(j = 0; j < dim[1]; j++)
      printf("%d ", arr[i][j][plane]);
    printf("]\n");
  }
}

void printMatrix3DInt(char * name, int * dim, int plane, int *** arr){
  int i, j;
  int r, c;
  r = dim[0];
  c = dim[1];
  printf("%s[%d][%d][%d]: \n", name, r, c, plane);
  for(i = 0; i < dim[0]; i++){
    printf("[ ");
    for(j = 0; j < dim[1]; j++)
      printf("%d ", arr[i][j][plane]);
    printf("]\n");
  }
}

void printMatrix3DUInt(char * name, int * dim, int plane, unsigned int *** arr){
  int i, j;
  int r, c;
  r = dim[0];
  c = dim[1];
  printf("%s[%d][%d][%d]: \n", name, r, c, plane);
  for(i = 0; i < dim[0]; i++){
    printf("[ ");
    for(j = 0; j < dim[1]; j++)
      printf("%d ", arr[i][j][plane]);
    printf("]\n"); 
  }
}

void printMatrix3DFloat(char * name, int * dim, int plane, float *** arr){
  int i, j;
  int r, c;
  r = dim[0];
  c = dim[1];
  printf("%s[%d][%d]: \n", name, r, c);
  for(i = 0; i < dim[0]; i++){
    printf("%d: [ ", i);
    for(j = 0; j < dim[1]; j++)
      printf("%e ", arr[i][j][plane]);
    printf("]\n");
  }
}

void printMatrix3DDouble(char * name, int * dim, int plane, double *** arr){
  int i, j;
  int r, c;
  r = dim[0];
  c = dim[1];
  printf("%s[%d][%d]: \n", name, r, c);
  for(i = 0; i < dim[0]; i++){
    printf("%d: [ ", i);
    for(j = 0; j < dim[1]; j++)
      printf("%e ", arr[i][j][plane]);
    printf("]\n");
  }
}

void setMatrixFloat(int r, int c, float val, float ** m){
  int i, j;
  for(i = 0; i < r; i++)
    for(j = 0; j < c; j++)
      m[i][j] = val;
}
void setMatrixDouble(int r, int c, double val, double ** m){
	int i, j;
	for(i = 0; i < r; i++)
		for(j = 0; j < c; j++)
			m[i][j] = val;
}

void setMatrixInt(int r, int c, int val, int ** m){
	int i, j;
	for(i = 0; i < r; i++)
		for(j = 0; j < c; j++)
			m[i][j] = val;
}

void setMatrixUInt(int r, int c, unsigned int val, unsigned int ** m){
	int i, j;
	for(i = 0; i < r; i++)
		for(j = 0; j < c; j++)
			m[i][j] = val;
}

void setMatrix3DFloat(float *** m, int r, int c, int z, float val){
  int i, j, k;
  for(i = 0; i < r; i++)
    for(j = 0; j < c; j++)
      for(k = 0; k < z; k++)
        m[i][j][k] = val;
}

void setMatrix3DDouble(double *** m, int r, int c, int z, double val){
	int i, j, k;
	for(i = 0; i < r; i++)
		for(j = 0; j < c; j++)
			for(k = 0; k < z; k++)
				m[i][j][k] = val;
}

void setMatrix3DInt(int *** m, int r, int c, int z, int val){
		int i, j, k;
	for(i = 0; i < r; i++)
		for(j = 0; j < c; j++)
			for(k = 0; k < z; k++)
				m[i][j][k] = val;
}

void setMatrix3DUInt(unsigned int *** m, int r, int c,  int z, unsigned int val){
		int i, j, k;
	for(i = 0; i < r; i++)
		for(j = 0; j < c; j++)
			for(k = 0; k < z; k++)
				m[i][j][k] = val;
}

void getMatrixColumnUInt(unsigned int ** mat, int r, int c, int col, unsigned int * column){
  int i;

  if(col > c || col < 0) {fputs ("getMatrixColumn: column selected out of bounds.",stderr); exit (2);}
  
  for(i = 0; i < r; i++)
    column[i] = mat[i][col];

}

void getMatrixColumnInt(int ** mat, int r, int c, int col, int * column){
  int i;

  if(col > c || col < 0) {fputs ("getMatrixColumn: column selected out of bounds.",stderr); exit (2);}
  
  for(i = 0; i < r; i++)
    column[i] = mat[i][col];

}

void getMatrixColumnFloat(float ** mat, int r, int c, int col, float * column){
  int i;

  if(col > c || col < 0) {fputs ("getMatrixColumn: column selected out of bounds.",stderr); exit (2);}
  
  for(i = 0; i < r; i++)
    column[i] = mat[i][col];

}

void getMatrixColumnDouble(double ** mat, int r, int c, int col, double * column){
  int i;

  if(col > c || col < 0) {fputs ("getMatrixColumn: column selected out of bounds.",stderr); exit (2);}
  
  for(i = 0; i < r; i++)
    column[i] = mat[i][col];

}


void getMatrixRowUInt(unsigned int ** mat, int r, int c, int rw, unsigned int * row){
  int i;

  if(rw > r || rw < 0) {fputs ("getMatrixRow: row selected out of bounds.",stderr); exit (2);}
  
  for(i = 0; i < c; i++)
    row[i] = mat[rw][i];

}

void getMatrixRowInt(int ** mat, int r, int c, int rw, int * row){
  int i;

  if(rw > r || rw < 0) {fputs ("getMatrixRow: row selected out of bounds.",stderr); exit (2);}
  
  for(i = 0; i < c; i++)
    row[i] = mat[rw][i];

}

void getMatrixRowFloat(float ** mat, int r, int c, int rw, float * row){
  int i;

  if(rw > r || rw < 0) {fputs ("getMatrixRow: row selected out of bounds.",stderr); exit (2);}
  
  for(i = 0; i < c; i++)
    row[i] = mat[rw][i];

}

void getMatrixRowDouble(double ** mat, int r, int c, int rw, double * row){
  int i;

  if(rw > r || rw < 0) {fputs ("getMatrixRow: row selected out of bounds.",stderr); exit (2);}
  
  for(i = 0; i < c; i++)
    row[i] = mat[rw][i];

}

void multiplyMatrixFloat(int r1, int c1, float ** m1, int r2, int c2, float ** m2, float ** prd){
  int c, d, k;
  float sum = 0;

  if(c1 != r2){
    printf("Matrices incompatible for multiplication col_1 %d != %d row_2\n", c1, r2);
    return;
  }

  for (c = 0; c < r1; c++) {
      for (d = 0; d < c2; d++) {
        for (k = 0; k < r2; k++) {
          //printf("m1[%d][%d] m2[%d][%d] ", c, k, k, d);
          sum = sum + m1[c][k]*m2[k][d];
        }
        //printf("\n prd[%d][%d] = %f\n", c, d, sum);
        prd[c][d] = sum;
        sum = 0;
      }
    }
}

void multiplyMatrixDouble(int r1, int c1, double ** m1, int r2, int c2, double ** m2, double ** prd){
  int c, d, k;
  double sum = 0;

  if(c1 != r2){
    printf("Matrices incompatible for multiplication col_1 %d != %d row_2\n", c1, r2);
    return;
  }

  for (c = 0; c < r1; c++) {
      for (d = 0; d < c2; d++) {
        for (k = 0; k < r2; k++) {
          sum = sum + m1[c][k]*m2[k][d];
        }
 
        prd[c][d] = sum;
        sum = 0;
      }
    }
}

void multiplyMatrixInt(int r1, int c1, int ** m1, int r2, int c2, int ** m2, int ** prd){
  int c, d, k;
  int sum = 0;

  if(c1 != r2){
    printf("Matrices incompatible for multiplication col_1 %d != %d row_2\n", c1, r2);
    return;
  }

  for (c = 0; c < r1; c++) {
      for (d = 0; d < c2; d++) {
        for (k = 0; k < r2; k++) {
          sum = sum + m1[c][k]*m2[k][d];
        }
 
        prd[c][d] = sum;
        sum = 0;
      }
    }
}

void multiplyMatrixUInt(int r1, int c1, unsigned int ** m1, int r2, int c2, unsigned int ** m2, unsigned int ** prd){
  int c, d, k;
  unsigned int sum = 0;

  if(c1 != r2){
    printf("Matrices incompatible for multiplication col_1 %d != %d row_2\n", c1, r2);
    return;
  }

  for (c = 0; c < r1; c++) {
      for (d = 0; d < c2; d++) {
        for (k = 0; k < r2; k++) {
          sum = sum + m1[c][k]*m2[k][d];
        }
 
        prd[c][d] = sum;
        sum = 0;
      }
    }
}

void multiplyMatrixVectorFloat(int r, int c, float ** m, float * v, int l,  float * prd){
  int i, j;
  float sum = 0;

  if(c != l){
    printf("Matrices incompatible for multiplication col_1 %d != %d row_2\n", c, r);
    return;
  }

  for (i = 0; i < r; i++) {
    for (j = 0; j < l; j++) 
      sum = sum + m[i][j]*v[j];
        
    prd[i] = sum;
    sum = 0;
  }
}

void multiplyMatrixVectorDouble(int r, int c, double ** m, double * v, int l,  double * prd){
  int i, j;
  double sum;

  if(c != l){
    printf("Matrices incompatible for multiplication col_1 %d != %d row_2\n", c, r);
    return;
  }

  for (i = 0; i < r; i++) {
    sum = 0;
    
  	for (j = 0; j < l; j++) 
  		sum += m[i][j]*v[j];
        
    prd[i] = sum;
  }
}

void multiplyMatrixVectorInt(int r, int c, int ** m, int * v, int l, int * prd){
  int i, j;
  int sum = 0;

  if(c != l){
    printf("Matrices incompatible for multiplication col_1 %d != %d row_2\n", c, r);
    return;
  }

  for (i = 0; i < r; i++) {
  	for (j = 0; j < l; j++) 
  		sum = sum + m[i][j]*v[j];
        
    prd[i] = sum;
    sum = 0;
  }
}

void multiplyMatrixVectorUInt(int r, int c, unsigned int ** m, unsigned int * v, int l, unsigned int * prd){
  int i, j;
  unsigned int sum = 0;

  if(c != l){
    printf("Matrices incompatible for multiplication col_1 %d != %d row_2\n", c, r);
    return;
  }

  for (i = 0; i < r; i++) {
  	for (j = 0; j < l; j++) 
  		sum = sum + m[i][j]*v[j];
        
    prd[i] = sum;
    sum = 0;
  }
}

void getPlane(double *** mat3d, int r, int c, int plane, double ** matout){
  int i, j;
  for(i = 0; i < r; i++)
    for(j = 0; j < c; j++)
      matout[i][j] = mat3d[i][j][plane];
}

void subtractMatrixDouble(int r, int c, double ** mat1, double ** mat2, double ** matout){
  int i, j;
  for(i = 0; i < r; i++)
    for(j = 0; j < c; j++)
      matout[i][j] = mat1[i][j] - mat2[i][j];
}

void subtractMatrixFloat(int r, int c, float ** mat1, float ** mat2, float ** matout){
  int i, j;
  for(i = 0; i < r; i++)
    for(j = 0; j < c; j++)
      matout[i][j] = mat1[i][j] - mat2[i][j];
}

float determinantf(float ** m, int n){
  int i, j, j1, j2;
    float det = 0;
    float ** temp;

   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      det = m[0][0];
   } else if (n == 2) {
      det = m[0][0] * m[1][1] - m[1][0] * m[0][1];
   } else {
      det = 0;
      for (j1 = 0; j1 < n; j1++) {
         temp = malloc((n - 1)*sizeof(float *));
         for (i = 0; i < n - 1; i++)
            temp[i] = malloc((n - 1)*sizeof(float));
         for (i = 1; i < n; i++) {
            j2 = 0;
            for (j=  0; j < n; j++) {
               if (j == j1)
                  continue;
               temp[i-1][j2] = m[i][j];
               j2++;
            }
         }
         det += pow(-1.0, j1 + 2.0) * m[0][j1] * determinantf(temp, n - 1);
         for (i = 0; i < n - 1; i++)
            free(temp[i]);
         free(temp);
      }
   }
   return det;
}

double determinantd(double ** m, int n){
	int i, j, j1, j2;
   	double det = 0;
   	double ** temp;
    //printf("determinantd: n = %d\n", n);
   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      det = m[0][0];
   } else if (n == 2) {
      det = m[0][0] * m[1][1] - m[1][0] * m[0][1];
   } else {
      det = 0;
      for (j1 = 0; j1 < n; j1++) {
         temp = malloc((n - 1)*sizeof(double *));
         for (i = 0; i < n - 1; i++)
            temp[i] = malloc((n - 1)*sizeof(double));
         for (i = 1; i < n; i++) {
            j2 = 0;
            for (j=  0; j < n; j++) {
               if (j == j1)
                  continue;
               temp[i-1][j2] = m[i][j];
               j2++;
            }
         }
         det += pow(-1.0, j1 + 2.0) * m[0][j1] * determinantd(temp, n - 1);
         for (i = 0; i < n - 1; i++)
            free(temp[i]);
         free(temp);
      }
   }
   return det;
}

void transposeFloat(float ** m, int r, int c, float ** trans){
  int i,j;

    for(i = 0; i < r; i++)
        for(j = 0; j < c; j++) 
          trans[j][i] = m[i][j];
}

void transposeDouble(double ** m, int r, int c, double ** trans){
	int i,j;

   	for(i = 0; i < r; i++)
      	for(j = 0; j < c; j++) 
         	trans[j][i] = m[i][j];
}

void transposeInt(int ** m, int r, int c, int ** trans){
	int i,j;

   	for(i = 0; i < r; i++)
      	for(j = 0; j < c; j++) 
         	trans[j][i] = m[i][j];
}

void transposeUInt(unsigned int ** m, int r, int c, unsigned int ** trans){
	int i,j;

   	for(i = 0; i < r; i++)
      	for(j = 0; j < c; j++) 
         	trans[j][i] = m[i][j];
}


void inversef(int n, float ** m, float ** inv){
  int i,j,ii,jj,i1,j1;
    float det;
    float **c;

    c = malloc((n-1)*sizeof(float *));
    for (i = 0; i < n - 1; i++)
      c[i] = malloc((n-1)*sizeof(float));

    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) {

          /* Form the adjoint a_ij */
          i1 = 0;
          for (ii = 0; ii < n; ii++) {
              if (ii == i)
                  continue;
              j1 = 0;
              for (jj = 0; jj < n; jj++) {
                  if (jj == j)
                      continue;
                  c[i1][j1] = m[ii][jj];
                  j1++;
              }
              i1++;
          }

          /* Calculate the determinate */
          det = determinantf(c, n-1);

          /* Fill in the elements of the cofactor */
          inv[i][j] = pow(-1.0, i + j + 2.0) * det;
      }
   }
   
   for (i=0;i<n-1;i++)
      free(c[i]);
   free(c);
}

void pinvf(int nr, int nc, float ** M, float ** Minv){
    // int nElements = nc*nr;
    // int i;
    
    // double ** V = calloc(nc, sizeof(double *));
    // for(i = 0; i < nc; i++)
    //   V[i] = calloc(nr, sizeof(double));

    // double W[nElements];

    // dsvd(M2D, nr, nc, W, V);
    // svdInvf(M2D, W, V, nr, nc, Minv2D);

    // free(V);

    int i;

    float ** Mt = calloc(nc, sizeof(float *));
    for(i = 0; i < nc; i++)
      Mt[i] = calloc(nr, sizeof(float));

    float ** MtM = calloc(nc, sizeof(float *));
    for(i = 0; i < nc; i++)
      MtM[i] = calloc(nc, sizeof(float));

    transposeFloat(M, nr, nc, Mt);
    multiplyMatrixFloat(nc, nr, Mt, nr, nc, M, MtM);

    float ** MtMinv = calloc(nc, sizeof(float *));
    for(i = 0; i < nc; i++)
      MtMinv[i] = calloc(nc, sizeof(double));
    
    LUinversionfloat(nc, MtM, MtMinv);
    multiplyMatrixFloat(nc, nc, MtMinv, nc, nr, Mt, Minv);

    free(Mt);
    free(MtM);
    free(MtMinv);
}

void pinvd(int nr, int nc, double ** M, double ** Minv){
    int i;

    double ** Mt = calloc(nc, sizeof(double *));
    for(i = 0; i < nc; i++)
      Mt[i] = calloc(nr, sizeof(double));

    double ** MtM = calloc(nc, sizeof(double *));
    for(i = 0; i < nc; i++)
      MtM[i] = calloc(nc, sizeof(double));

    transposeDouble(M, nr, nc, Mt);
    multiplyMatrixDouble(nc, nr, Mt, nr, nc, M, MtM);

    double ** MtMinv = calloc(nc, sizeof(double *));
    for(i = 0; i < nc; i++)
      MtMinv[i] = calloc(nc, sizeof(double));
    
    LUinversion(nc, MtM, MtMinv);
    multiplyMatrixDouble(nc, nc, MtMinv, nc, nr, Mt, Minv);

    free(Mt);
    free(MtM);
    free(MtMinv);
}

void mean_2Df(float ** m, int h, int w, float * r){

    float sum = 0;
    int i,j;
    for(j = 0; j < w; j++){
        for(i = 0; i < h; i++)
            sum = sum + m[i][j];

        r[j] = sum/w;
        sum = 0;
    }
}

void mean_2Dd(double ** m, int h, int w, double * r){

    double sum = 0;
    int i,j;
    for(j = 0; j < w; j++){
        for(i = 0; i < h; i++)
            sum = sum + m[i][j];

        r[j] = sum/w;
        sum = 0;
    }
}

 void interp2(unsigned int ** input, int rows, int cols, double * x, double * y, double sample_size, double * output){

  int i;
  for(i = 0; i < sample_size; i++){
    output[i] = bilininterp(input, rows, cols, x[i], y[i]);
  }

}

double bilininterp(unsigned int ** input, int rows, int cols, double x, double y){

  int x1, x2, y1, y2;
  double Q11, Q12, Q21, Q22;
  double value;

  //printf("(x,y) = %f, %f\n", x, y);

  x1 = floor(x); //printf("x1 = %d\n", x1);
  x2 = x1 + 1;   //printf("x2 = %d\n", x2);

  if( x1 < 0) x1 = 0;
  if( x2 >= cols) x2 = cols - 1;

  y1 = floor(y); //printf("y1 = %d\n", y1);
  y2 = y1 + 1;   //printf("y2 = %d\n", y2);

  if( y1 < 0) y1 = 0;
  if( y2 >= rows) y2 = rows - 1;

  if(x1 <= 0 || x2 <= 0 || y1 <= 0 || y2 <= 0)
    return 0;

 if(x1 > cols || x2 > cols || y1 > rows || y2 > rows)
    return 0;
  
  Q11 = (double) input[x1 - 1][y1 - 1]; //printf("Q11 = %f\n", Q11);
  Q12 = (double) input[x1 - 1][y2 - 1]; //printf("Q12 = %f\n", Q11);
  Q21 = (double) input[x2 - 1][y1 - 1]; //printf("Q21 = %f\n", Q11);
  Q22 = (double) input[x2 - 1][y2 - 1]; //printf("Q23 = %f\n", Q11);

  value = Q11*((double)x2 - x)*((double)y2 - y) + Q21*(x - (double)x1)*((double)y2 - y) + Q12*((double)x2 - x)*(y - (double)y1) + Q22*(x - (double)x1)*(y - (double)y1);

  return value;

}

void convert2D21Dfloat(int r, int c, float ** mat2D, float * mat1D){

    int i,j;
    for(i = 0; i < r; i++)
        for(j = 0; j < c; j++)
            mat1D[i*c + j] = mat2D[i][j];

    return;

}

void convert1D22Dfloat(int r, int c, float * mat1D, float ** mat2D){

    int i,j;
    for(i = 0; i < r; i++)
        for(j = 0; j < c; j++)
           mat2D[i][j] = mat1D[i*c + j];

    return;
    
}

void convert2D21Ddouble(int r, int c, double ** mat2D, double * mat1D){

    int i,j;
    for(i = 0; i < r; i++)
        for(j = 0; j < c; j++)
            mat1D[i*c + j] = mat2D[i][j];

    return;

}

void convert1D22Ddouble(int r, int c, double * mat1D, double ** mat2D){

    int i,j;
    for(i = 0; i < r; i++)
        for(j = 0; j < c; j++)
           mat2D[i][j] = mat1D[i*c + j];

    return;
    
}

void convert2D21Dint(int r, int c, int (* mat2D)[c], int * mat1D){

    int i,j;
    for(i = 0; i < r; i++)
        for(j = 0; j < c; j++)
            mat1D[i*c + j] = mat2D[i][j];

    return;

}

void convert1D22Dint(int r, int c, int * mat1D, int (* mat2D)[c]){

    int i,j;
    for(i = 0; i < r; i++)
        for(j = 0; j < c; j++)
           mat2D[i][j] = mat1D[i*c + j];

    return;
    
}

void convert2D21Duint(int r, int c, unsigned int ** mat2D, unsigned int * mat1D){

    int i,j;
    for(i = 0; i < r; i++)
        for(j = 0; j < c; j++)
            mat1D[i*c + j] = mat2D[i][j];

    return;

}

void convert1D22Duint(int r, int c, unsigned int * mat1D, unsigned int ** mat2D){

    int i,j;
    for(i = 0; i < r; i++)
        for(j = 0; j < c; j++)
           mat2D[i][j] = mat1D[i*c + j];

    return;
    
}

void convert2D21Dushort(int r, int c, unsigned short ** mat2D, unsigned short * mat1D){

    int i,j;
    for(i = 0; i < r; i++)
        for(j = 0; j < c; j++)
            mat1D[i*c + j] = mat2D[i][j];

    return;

}

void convert1D22Dushort(int r, int c, unsigned short * mat1D, unsigned short ** mat2D){

    int i,j;
    for(i = 0; i < r; i++)
        for(j = 0; j < c; j++)
           mat2D[i][j] = mat1D[i*c + j];

    return;
    
}

double min2Dd(int r, int c, double ** mat){
    int i, j;
    double min_arr[r];
    double rowarr[c];

    for(i = 0; i < r; i++){
      for(j = 0; j < c; j++)
        rowarr[j] = mat[i][j];
      min_arr[i] = minarrayd(rowarr, c);
    }

    return minarrayd(min_arr, r);
}

double max2Dd(int r, int c, double ** mat){
    int i, j;
    double max_arr[r];
    double rowarr[c];

    for(i = 0; i < r; i++){
      for(j = 0; j < c; j++)
        rowarr[j] = mat[i][j];
      max_arr[i] = maxarrayd(rowarr, c);
    }

    return maxarrayd(max_arr, r);
}

double sum2Dd(int r, int c, double ** mat){
  int i, j;
  double sum = 0;
  for(i = 0; i < r; i++)
      for(j = 0; j < c; j++)
        sum = sum + mat[i][j];

  return sum;
}

double std2Dd(int r, int c, double ** mat){
    int i, j;
    double mat_sum = sum2Dd(r, c, mat);
    double mat_avg = mat_sum/(r*c);
    double diff = 0;

    for(i = 0; i < r; i++)
      for(j = 0; j < c; j++)
        diff = fabs(mat_avg - mat[i][j]);

    return diff/(r*c);
}

float sum2Df(int r, int c, float ** mat){
  int i, j;
  float sum = 0;
  for(i = 0; i < r; i++)
      for(j = 0; j < c; j++)
        sum = sum + mat[i][j];

  return sum;
}

float std2Df(int r, int c, float ** mat){
    int i, j;
    float mat_sum = sum2Df(r, c, mat);
    float mat_avg = mat_sum/(r*c);
    float diff = 0;

    for(i = 0; i < r; i++)
      for(j = 0; j < c; j++)
        diff = fabsf(mat_avg - mat[i][j]);

    return diff/(r*c);
}

float nonzerostd2Df(int r, int c, float ** mat){
    int i, j, num;
    float mat_sum = sum2Df(r, c, mat);
    float mat_avg;
    float diff = 0;

    num = 0;

    for(i = 0; i < r; i++)
      for(j = 0; j < c; j++)
        if(mat[i][j] != 0)
          num++;

    mat_avg = mat_sum/num;

    for(i = 0; i < r; i++)
      for(j = 0; j < c; j++)
        if(mat[i][j] != 0)
          diff += fabsf(mat_avg - mat[i][j]);      
      

    return diff/(num);
}

double nonzerostd2Dd(int r, int c, double ** mat){
    int i, j, num;
    double mat_sum = sum2Dd(r, c, mat);
    double mat_avg;
    double diff = 0;

    num = 0;

    for(i = 0; i < r; i++)
      for(j = 0; j < c; j++)
        if(mat[i][j] != 0)
          num++;

    mat_avg = mat_sum/num;

    for(i = 0; i < r; i++)
      for(j = 0; j < c; j++)
        if(mat[i][j] != 0)
          diff += fabs(mat_avg - mat[i][j]);      
      

    return diff/(num);
}

void fliplr2D(int r, int c, double ** mat){
  int i, j, k;
  double temp;
  for(k = 0; k < r; k++){
    j = c - 1;
    for(i = 0; i < j; i++){
        temp = mat[k][i];
        mat[k][i] = mat[k][j];
        mat[k][j] = temp;
        j--;
    }
  }
}



void getRowUInt(int irow, unsigned int ** mat, int rows, int cols, unsigned int * row){
  int i;
  for(i = 0; i < cols; i++)
    row[i] = mat[irow][i];
}

void getRowInt(int irow, int ** mat, int rows, int cols, int * row){
  int i;
  for(i = 0; i < cols; i++)
    row[i] = mat[irow][i];
}

void getRowFloat(int irow, float ** mat, int rows, int cols, float * row){
  int i;
  for(i = 0; i < cols; i++)
    row[i] = mat[irow][i];
}

void getRowDouble(int irow, double ** mat, int rows, int cols, double * row){
  int i;
  for(i = 0; i < cols; i++)
    row[i] = mat[irow][i];
}


void nonzero_elements_matrix_uint(int r, int c, unsigned int ** mat, unsigned int * elements, int * size){
  int i, j;
  int ind = 0;

  for(i = 0; i < r; i++)
    for(j = 0; j < c; j++)
      if(mat[i][j] != 0){
        elements[ind] = mat[i][j];
        ind++;
      }

  * size = ind;
}

void nonzero_elements_matrix_double(int r, int c, double ** mat, double * elements, int * size){
  int i, j;
  int ind = 0;

  for(i = 0; i < r; i++)
    for(j = 0; j < c; j++)
      if(mat[i][j] != 0){
        elements[ind] = mat[i][j];
        ind++;
      }

  * size = ind;
}

void nonzero_elements_matrix_float(int r, int c, float ** mat, float * elements, int * size){
  int i, j;
  int ind = 0;

  for(i = 0; i < r; i++)
    for(j = 0; j < c; j++)
      if(mat[i][j] != 0){
        elements[ind] = mat[i][j];
        ind++;
      }

  * size = ind;
}

void nonzero_elements_matrix_int(int r, int c, int ** mat, int * elements, int * size){
  int i, j;
  int ind = 0;

  for(i = 0; i < r; i++)
    for(j = 0; j < c; j++)
      if(mat[i][j] != 0){
        elements[ind] = mat[i][j];
        ind++;
      }

  * size = ind;
}

void fliplruint(int r, int c, unsigned int ** mat){
  int i, j;
  unsigned int temp[r];
  for(i = 0; i < c/2; i++){
    if(i == c - i) return;
    getMatrixColumnUInt(mat, r, c, i, temp);
    for(j = 0; j < r; j++){
      mat[j][i] = mat[j][c - i]; 
      mat[j][c - i] = temp[j];
    }
  }
}

void fliplrint(int r, int c, int ** mat){
  int i, j;
  int temp[r];
  for(i = 0; i < c/2; i++){
    if(i == c - i) return;
    getMatrixColumnInt(mat, r, c, i, temp);
    for(j = 0; j < r; j++){
      mat[j][i] = mat[j][c - i]; 
      mat[j][c - i] = temp[j];
    }
  }
}

void fliplrfloat(int r, int c, float ** mat){
  int i, j;
  float temp[r];
  for(i = 0; i < c/2; i++){
    if(i == c - i) return;
    getMatrixColumnFloat(mat, r, c, i, temp);
    for(j = 0; j < r; j++){
      mat[j][i] = mat[j][c - i]; 
      mat[j][c - i] = temp[j];
    }
  }
}

void fliplrdouble(int r, int c, double ** mat){
  int i, j;
  double temp[r];
  for(i = 0; i < c/2; i++){
    if(i == c - i) return;
    getMatrixColumnDouble(mat, r, c, i, temp);
    for(j = 0; j < r; j++){
      mat[j][i] = mat[j][c - i]; 
      mat[j][c - i] = temp[j];
    }
  }
}
