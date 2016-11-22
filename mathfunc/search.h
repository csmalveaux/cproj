#ifndef SEARCH_H__
#define SEARCH_H__

int linearsearchuint(unsigned int * arr, int l, unsigned int val, int * index);
int linearsearchint(int * arr, int l, int val, int * index);
int linearsearchfloat(float * arr, int l, float val, float tol, int * index);
int linearsearchdouble(double * arr, int l, double val, double tol, int * index);

int binarysearchuint(unsigned int * arr, int l, unsigned int val, int * index);
int binarysearchint(int * arr, int l, int val, int * index);
int binarysearchfloat(float * arr, int l, float val, float tol, int * index);
int binarysearchdouble(double * arr, int l, double val, double tol, int * index);

int interpolationsearchuint(unsigned int * arr, int l, unsigned int val, int * index);
int interpolationsearchint(int * arr, int l, int val, int * index);
int interpolationsearchfloat(float * arr, int l, float val, float tol, int * index);
int interpolationsearchdouble(double * arr, int l, double val, double tol, int * index);

int findminuint(unsigned int * arr, int l, unsigned int * val, int * index);
int findminint(int * arr, int l, int * val, int  * index);
int findminfloat(float * arr, int l, float * val, int * index);
int findmindouble(double * arr, int l, double * val, int * index);

int findmaxuint(unsigned int * arr, int l, unsigned int * val, int * index);
int findmaxint(int * arr, int l, int * val, int * index);
int findmaxfloat(float * arr, int l, float * val, int * index);
int findmaxdouble(double * arr, int l, double * val, int * index);

int findminmaxuint(unsigned int * arr, int l, int offset, int * minindex, int * imin, int * maxindex, int * imax, unsigned int * min_val, unsigned int * max_val);
int findminmaxint(int * arr, int l, int offset, int * minindex, int * imin, int * maxindex, int * imax, int * min_val, int * max_val);
int findminmaxfloat(float * arr, int l, int offset, int * minindex, int * imin, int * maxindex, int * imax, float * min_val, float * max_val);
int findminmaxdouble(double * arr, int l, int offset, int * minindex, int * imin, int * maxindex, int * imax, double * min_val, double * max_val);

#endif /* SEARCH_H__ */
