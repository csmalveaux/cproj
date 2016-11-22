
#ifndef ARRAY_H__
#define ARRAY_H__

#include "matrix.h"

#define SIGN(a,b) ((b) >= 0. ? fabs(a) : -fabs(a))

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

void setarrayfloat(float * arr, int l, float val);
void setarraydouble(double * arr, int l, double val);
void setarrayint(int * arr, int l, int val);
void setarrayuint(unsigned int * arr, int l, unsigned int val);

void everyotheruint(unsigned int * arr, int l, int offset, int gap, unsigned int * set);
void everyotherint(int * arr, int l, int offset, int gap, int * set);
void everyotherfloat(float * arr, int l, int offset, int gap, float * set);
void everyotherdouble(double * arr, int l, int offset, int gap, double * set);

void seteveryotheruint(unsigned int * set, int l, int offset, int gap, unsigned int * arr);
void seteveryotherint(int * set, int l, int offset, int gap, int * arr);
void seteveryotherfloat(float * set, int l, int offset, int gap, float * arr);
void seteveryotherdouble(double * set, int l, int offset, int gap, double * arr);

unsigned int maxarrayuint(unsigned int * a, int size);
int   maxarrayint(int * a, int size);
double maxarrayd(double * a, int size);
float maxarrayf(float * a, int size);

unsigned int minarrayuint(unsigned int * a, int size);
int   minarrayint(int * a, int size);
double minarrayd(double * a, int size);
float minarrayf(float * a, int size);

void sortf(float * arr, int l);
void sortd(double * arr, int l);

void sortindexf(float * arr, int l, int * index);
void sortindexd(double * arr, int l, int * index);

float medianf(float * arr, int l);
double mediand(double * arr, int l);

void median2Ddouble(double ** m, int row, int w, double * medi);

void medfilt1(double * s, int n, int w, double * mr);

void arraydiffd(double * arr, int size, double * diff);
void arraydifff(float * arr, int size, float * diff);
void arraydiffint(int * arr, int size, int * diff);
void arraydiffuint(unsigned int * arr, int size, unsigned int * diff);

void flipuduint(int n, unsigned int * arr);
void flipudint(int n, int * arr);
void flipuddouble(int n, double * arr);
void flipudfloat(int n, float * arr);

void crossf(float * a, float * b, float * c);
void crossd(double * a, double * b, double * c);

float absf(float * a);
double absd(double * a);

void normf(float * a, float * b);
void normd(double * a, double * b);

double trapz(double * x, double * y, int size);
void cumsum(double * a, double * b, int size);
double sum(double * a, int size);
unsigned int sumuint(unsigned int * a, int size);
void find(int *a, int size, int * indArray, int * num);
void findf(float * a, int size, int * indArray, int * num);
void findd(double * a, int size, int * indArray, int * num);
void diff(double * x, int size, double * y);
int any(int * a, int size, int * num);
int sign(double input);

double std(double * arr, int size);
float stdf(float * arr, int size);

void multiply(double * x, double * y, int size, double * z);

void printArrayShort(char * name, int length, short * arr);
void printArrayUShort(char * name, int length, unsigned short * arr);
void printArrayInt(char * name, int length, int * arr);
void printArrayUInt(char * name, int length, unsigned int * arr);
void printArrayDouble(char * name, int length, double * arr);
void printArrayFloat(char * name, int length, float * arr);

void abs_int(int * arr, int length);
void abs_float(float * arr, int length);
void abs_double(double * arr, int length);

void move_elements_end_double(double * arr, int * index, int l, int * beg, int ibeg, int * end, int iend);
void find_all_inf(double * arr,  int l, int * indexinf, int * i_inf, int * index, int * ind);

void move_elements_end_float(float * arr, int * index, int l, int * beg, int ibeg, int * end, int iend);
void find_all_inf_float(float * arr,  int l, int * indexinf, int * i_inf, int * index, int * ind);

#endif /* ARRAY_H__ */
