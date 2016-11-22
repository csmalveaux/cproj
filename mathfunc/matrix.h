
#ifndef MATRIX_H__
#define MATRIX_H__

#include <stdio.h>

#define EPS 			   2.2204e-16

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

void nonzero_elements_matrix_uint(int r, int c, unsigned int ** mat, unsigned int * elements, int * size);
void nonzero_elements_matrix_int(int r, int c, int ** mat, int * elements, int * size);
void nonzero_elements_matrix_float(int r, int c, float ** mat, float * elements, int * size);
void nonzero_elements_matrix_double(int r, int c, double ** mat, double * elements, int * size);

void fliplruint(int r, int c, unsigned int ** mat);
void fliplrint(int r, int c, int ** mat);
void fliplrfloat(int r, int c, float ** mat);
void fliplrdouble(int r, int c, double ** mat);

#endif /* MATRIX_H__ */
