
#include "array.h"
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#include "sort.h"
#include "matrix.h"
#include "binarytree.h"
#include "search.h"

// #define MIN(a,b) (((a)<(b))?(a):(b))
// #define MAX(a,b) (((a)>(b))?(a):(b))

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
void multiply(double * x, double * y, int size, double * z);

double std(double * arr, int size);
float stdf(float * arr, int size);

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


void setarrayfloat(float * arr, int l, float val){
    int i;
    for(i = 0; i < l; i++)
        arr[i] = val;
}

void setarraydouble(double * arr, int l, double val){
	int i;
	for(i = 0; i < l; i++)
		arr[i] = val;
}

void setarrayint(int * arr, int l, int val){
	int i;
	for(i = 0; i < l; i++)
		arr[i] = val;
}

void setarrayuint(unsigned int * arr, int l, unsigned int val){
	int i;
	for(i = 0; i < l; i++)
		arr[i] = val;
}

void everyotheruint(unsigned int * arr, int l, int offset, int gap, unsigned int * set){
    int i, j, next;

    for(i = offset, j = 0, next = offset; i < l; i++){
        if(i == next) {set[j] = arr[i]; next += gap; j++;}
    }
}

void everyotherint(int * arr, int l, int offset, int gap, int * set){
    int i, j, next;

    for(i = offset, j = 0, next = offset; i < l; i++){
        if(i == next) {set[j] = arr[i]; next += gap; j++;}
    }
}

void everyotherfloat(float * arr, int l, int offset, int gap, float * set){
    int i, j, next;

    for(i = offset, j = 0, next = offset; i < l; i++){
        if(i == next) {set[j] = arr[i]; next += gap; j++;}
    }
}

void everyotherdouble(double * arr, int l, int offset, int gap, double * set){
    int i, j, next;

    for(i = offset, j = 0, next = offset; i < l; i++){
        if(i == next) {set[j] = arr[i]; next += gap; j++;}

    }
}

void seteveryotheruint(unsigned int * set, int l, int offset, int gap, unsigned int * arr){
    int i;
    for(i = 0; i < l; i++){
        arr[offset + i*gap] = set[i];
    }
}

void seteveryotherint(int * set, int l, int offset, int gap, int * arr){
    int i;
    for(i = 0; i < l; i++){
        arr[offset + i*gap] = set[i];
    }
}

void seteveryotherfloat(float * set, int l, int offset, int gap, float * arr){
    int i;
    for(i = 0; i < l; i++){
        arr[offset + i*gap] = set[i];
    }
}

void seteveryotherdouble(double * set, int l, int offset, int gap, double * arr){
    int i;
    for(i = 0; i < l; i++){
        arr[offset + i*gap] = set[i];
    }
}

unsigned int maxarrayuint(unsigned int * a, int size){
    unsigned int m = a[0];
    int i;
    for(i = 1; i < size; i++)
        if(a[i] > m) 
            m = a[i];
    
    return m;
}

int   maxarrayint(int * a, int size){
    int m = a[0];
    int i;
    for(i = 1; i < size; i++)
        if(a[i] > m) 
            m = a[i];
    
    return m;
}

double maxarrayd(double * a, int size){
    double m = a[0];
    int i;
    for(i = 1; i < size; i++)
        if(a[i] > m) 
            m = a[i];
    
    return m;
}

float maxarrayf(float * a, int size){
    float m = a[0];
    int i;
    for(i = 1; i < size; i++)
        if(a[i] > m) 
            m = a[i];
    
    return m;
}

unsigned int minarrayuint(unsigned int * a, int size){
    unsigned int m = a[0];
    int i;
    for(i = 1; i < size; i++)
        if(a[i] < m) 
            m = a[i];
    
    return m;
}

int   minarrayint(int * a, int size){
    int m = a[0];
    int i;
    for(i = 1; i < size; i++)
        if(a[i] < m) 
            m = a[i];
    
    return m;
}

double minarrayd(double * a, int size){
    double m = a[0];
    int i;
    for(i = 1; i < size; i++)
        if(a[i] < m) 
            m = a[i];
    
    return m;
}

float minarrayf(float * a, int size){
    float m = a[0];
    int i;
    for(i = 1; i < size; i++)
        if(a[i] < m) 
            m = a[i];
    
    return m;
}

void sortf(float * arr, int l){
    int i;
    float temp;
    int order = 1;

    while(order != 0){
        order = 0;
        for(i = 0; i < l - 1; i++){
            if(MIN(arr[i], arr[i + 1]) == arr[i + 1]){
                order = 1;
                temp = arr[i];
                arr[i] = arr[i + 1];
                arr[i + 1] = temp;
            }
        }

    }
}


void sortd(double * arr, int l){
	int i, j, imin;

    for(i = 0; i < l; i++){

        double temparr[l - i];
        for(j = 0; j < l - i; j++)
            temparr[j] = arr[j + i];

        //imin = firstminindexdouble(temparr, l - i) + i;
        findmindouble(temparr, l + i, NULL, &imin);
        imin += i;

        switchdouble(arr, l, i, imin);
    }
}

void sortindexf(float * arr, int l, int * index){

    float temp[l];

    memcpy(temp, arr, l*sizeof(float));

    int temp_ind[l], temp_inf[l];
    int itemp_ind, itemp_inf;

    find_all_inf_float(arr, l, temp_inf, &itemp_inf, temp_ind, &itemp_ind);
    move_elements_end_float(arr, index, l, temp_ind, itemp_ind, temp_inf, itemp_inf);

    selectionsortfloat(arr, itemp_ind);
    sortedarraymapfloat(temp, arr, l, index);

    memcpy(arr, temp, l);
}


void sortindexd(double * arr, int l, int * index){

    double temp[l];

    memcpy(temp, arr, l*sizeof(double));

    int temp_ind[l], temp_inf[l];
    int itemp_ind, itemp_inf;

    find_all_inf(arr, l, temp_inf, &itemp_inf, temp_ind, &itemp_ind);
    move_elements_end_double(arr, index, l, temp_ind, itemp_ind, temp_inf, itemp_inf);

    selectionsortdouble(arr, itemp_ind);
    sortedarraymapdouble(temp, arr, l, index);

    memcpy(arr, temp, l);

}

float medianf(float * arr, int l){
    float temp[l];
    memcpy(temp, arr, sizeof(float)*l);
    selectionsortfloat(temp, l);
    return temp[l/2];
}

double mediand(double * arr, int l){
    double temp[l];
    memcpy(temp, arr, sizeof(double)*l);
    selectionsortdouble(temp, l);
    return temp[l/2];
}

void median2Ddouble(double ** m, int row, int w, double * medi){
    int i, j;
    double * temp = calloc(w, sizeof(double));

    for(i = 0; i < row; i++){
        for(j = 0; j < w; j++)
            temp[j] = m[i][j];
        medi[i] = mediand(temp, w);
    }

    free(temp);
}

void medfilt1(double * s, int n, int w, double * mr){

    int i,j,k = 0;
    int w2 = w/2;
    w = 2*w2 + 1;

    double ** m = calloc((n + w - 1), sizeof(double *));
    for(i = 0; i < (n + w - 1); i++)
        m[i] = calloc(w, sizeof(double));
    double * mdfilt = calloc(n + w - 1, sizeof(double));

    double s0 = s[0];
    double s1 = s[n - 1];

    for(i = 0; i < w; i++){
        for(j = 0; j < (n + w - 1); j++){
            if(i > 0 && j < i)
                m[j][i] = s0;
            else if(k < n){
                m[j][i] = s[k];
                k++;
            }
            else
                m[j][i] = s1;
        }
    }

    median2Ddouble(m, (n + w - 1), w, mdfilt);

    for(i = 0; i < n; i++)
        mr[i] = mdfilt[w2 + i];

    free(m);
    free(mdfilt);
}

void arraydiffd(double * arr, int size, double * diff){
  int i;
  for(i = 0; i < size - 1; i++)
    diff[i] = arr[i+1] - arr[i];
}

void arraydifff(float * arr, int size, float * diff){
  int i;
  for(i = 0; i < size - 1; i++)
    diff[i] = arr[i+1] - arr[i];
}

void arraydiffint(int * arr, int size, int * diff){
  int i;
  for(i = 0; i < size - 1; i++)
    diff[i] = arr[i+1] - arr[i];
}

void arraydiffuint(unsigned int * arr, int size, unsigned int * diff){
  int i;
  for(i = 0; i < size - 1; i++)
    diff[i] = arr[i+1] - arr[i];
}

void flipuduint(int n, unsigned int * arr){
    int i, j;
    unsigned int temp;
    j = n - 1;
    for(i = 0; i < j; i++){
        temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
        j--;
    }
}

void flipudint(int n, int * arr){
    int i, j;
    int temp;
    j = n - 1;
    for(i = 0; i < j; i++){
        temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
        j--;
    }
}

void flipuddouble(int n, double * arr){
    int i, j;
    double temp;
    j = n - 1;
    for(i = 0; i < j; i++){
        temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
        j--;
    }
}

void flipudfloat(int n, float * arr){
    int i, j;
    float temp;
    j = n - 1;
    for(i = 0; i < j; i++){
        temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
        j--;
    }
}

void crossf(float * a, float * b, float * c){
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

void crossd(double * a, double * b, double * c){
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

float absf(float * a){
    return sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2));
}

double absd(double * a){
    return sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2));
}

void normf(float * a, float * b){
    float n;
    n = absf(a);
    b[0] = a[0]/n;
    b[1] = a[1]/n;
    b[2] = a[2]/n;
}

void normd(double * a, double * b){
    double n;
    n = absd(a);
    b[0] = a[0]/n;
    b[1] = a[1]/n;
    b[2] = a[2]/n;
}

double trapz(double * x, double * y, int size){
    int i;
    double res;

    double * a = malloc(sizeof(double) * size);
    double * b = malloc(sizeof(double) * size);
    double * c = malloc(sizeof(double) * size);

    diff(x, size, a);
  
    for(i = 0; i < size; i++)
        b[i] = y[i] + y[i+1];

    multiply(a, b, size, c);

    res = sum(c, size);

    free(a);
    free(b);
    free(c);

    return res;
}

void cumsum(double * a, double * b, int size){
    int i;
    double element = 0;

    for(i = 0; i < size; i++){
        b[i] = element + a[i];
    }
}

double sum(double * a, int size){
    int i;
    double ret = 0;
    for(i = 0; i < size; i++)
        ret = ret + a[i];

    return ret;
}

unsigned int sumuint(unsigned int * a, int size){
    int i;
    unsigned int ret = 0;
    for(i = 0; i < size; i++)
        ret = ret + a[i];
    return ret;
}

void findd(double * a, int size, int * indArray, int * num){
    int i, indx;
    indx = 0;
    for(i = 0; i < size; i++)
        if(a[i] > 0 || a[i] < 0){
            indArray[indx] = i;
            indx++;
        }

    * num = indx;
}

void findf(float * a, int size, int * indArray, int * num){
    int i, indx;
    indx = 0;
    for(i = 0; i < size; i++)
        if(a[i] > 0 || a[i] < 0){
            indArray[indx] = i;
            indx++;
        }

    * num = indx;
}

void find(int *a, int size, int * indArray, int * num){
    int i, indx;
    indx = 0;

    for(i = 0; i < size; i++)
        if(a[i] > 0){
            indArray[indx] = i;
            indx++;
        }

    * num = indx;
}

void diff(double * x, int size, double * y){
  int i;
  for(i = 0; i < size; i++)
    y[i] = x[i+1] - x[i];
}

int any(int * a, int size, int * num){
    int i, n;
    int ret = 0;
    n = 0;
    for(i = 0; i < size; i++)
        if(a[i] == 1){
            ret = 1;
            n++;
        }

    * num = n;

    return ret;
}

int sign(double input){
    if(input > 0)
        return 1;
    else if(input == 0)
        return 0;
    else if(input < 0)
        return -1;

    return -1;
}

double std(double * arr, int l){
    int i;
    double mean, sum;

    for(i = 0; i < l; i++)
        mean += arr[i];

    mean /= l;

    for(i = 0; i < l; i++)
        sum += (arr[i] - mean) * (arr[i] - mean);

    return sqrt(sum/l);  
}

float stdf(float * arr, int l){
    int i;
    float mean, sum;

    for(i = 0; i < l; i++)
        mean += arr[i];

    mean /= l;

    for(i = 0; i < l; i++)
        sum += (arr[i] - mean) * (arr[i] - mean);

    return sqrt(sum/l);  
}

void multiply(double * x, double * y, int size, double * z){
  int i;
  for(i = 0; i < size; i++)
    z[i] = x[i] * y[i];
}

void printArrayShort(char * name, int length, short * arr){
  int i;
  printf("%s: ", name);
  for(i = 0; i < length; i++)
    printf("%d ", arr[i]);
  printf("\n");
}

void printArrayUShort(char * name, int length, unsigned short * arr){
  int i;
  printf("%s: ", name);
  for(i = 0; i < length; i++)
    printf("%x ", arr[i]);
  printf("\n");
}

void printArrayInt(char * name, int length, int * arr){
  int i;
  printf("%s: ", name);
  for(i = 0; i < length; i++)
    printf("%d ", arr[i]);
  printf("\n");
}

void printArrayUInt(char * name, int length, unsigned int * arr){
  int i;
  printf("%s: ", name);
  for(i = 0; i < length; i++)
    printf("%u ", arr[i]);
  printf("\n");
}

void printArrayFloat(char * name, int length, float * arr){
  int i;
  printf("%s: ", name);
  for(i = 0; i < length; i++)
    printf("%f ", arr[i]);
  printf("\n");
}

void printArrayDouble(char * name, int length, double * arr){
  int i;
  printf("%s: ", name);
  for(i = 0; i < length; i++)
    printf("%f ", arr[i]);
  printf("\n");
}

void abs_int(int * arr, int l){
    int i;
    for(i = 0; i < l; i++)
        arr[i] = abs(arr[i]);
}

void abs_float(float * arr, int l){
    int i;
    for(i = 0; i < l; i++)
        arr[i] = fabsf(arr[i]);
}

void abs_double(double * arr, int l){
    int i;
    for(i = 0; i < l; i++)
        arr[i] = fabs(arr[i]);
}

void find_all_inf(double * arr,  int l, int * indexinf, int * i_inf, int * index, int * ind){
    int i, ireg = 0, iinf = 0;

    for(i = 0; i < l; i++)
        if(isinf(arr[i]) == 1){
            indexinf[iinf] = i;
            iinf++;
        }
        else{
            index[ireg] = i;
            ireg++;
        }

    * i_inf = iinf;
    * ind = ireg;
}

void move_elements_end_double(double * arr, int * index, int l, int * beg, int size_beg, int * end, int size_end){
    int i;
    double arr_temp[l];

    if(l == size_beg || l == size_end)
        return;

    for(i = 0; i < size_beg; i++){
        index[i] = beg[i];
        arr_temp[i] = arr[beg[i]];
    }

    for(i = size_beg; i < size_beg + size_end; i++){
        index[i] = end[i - size_beg];
        arr_temp[i] = arr[end[i - size_beg]];
    }
    
    for(i = 0; i < l; i++)
        arr[i] = arr_temp[i];

}

void move_elements_end_float(float * arr, int * index, int l, int * beg, int size_beg, int * end, int size_end){
    int i;
    float arr_temp[l];

    if(l == size_beg || l == size_end)
        return;

    for(i = 0; i < size_beg; i++){
        index[i] = beg[i];
        arr_temp[i] = arr[beg[i]];
    }

    for(i = size_beg; i < size_beg + size_end; i++){
        index[i] = end[i - size_beg];
        arr_temp[i] = arr[end[i - size_beg]];
    }
    
    for(i = 0; i < l; i++)
        arr[i] = arr_temp[i];

}

void find_all_inf_float(float * arr,  int l, int * indexinf, int * i_inf, int * index, int * ind){
    int i, ireg = 0, iinf = 0;

    for(i = 0; i < l; i++)
        if(isinf(arr[i]) == 1){
            indexinf[iinf] = i;
            iinf++;
        }
        else{
            index[ireg] = i;
            ireg++;
        }

    * i_inf = iinf;
    * ind = ireg;
}
