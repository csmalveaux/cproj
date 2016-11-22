#ifndef SORT_H__
#define SORT_H__

void switchuint(unsigned int * arr, int size, int index1, int index2);
void switchint(int * arr, int size, int index1, int index2);
void switchfloat(float * arr, int size, int index1, int index2);
void switchdouble(double * arr, int size, int index1, int index2);

void sortedarraymapuint(unsigned int * unsorted, unsigned int * sorted, int size, int * index);
void sortedarraymapint(int * unsorted, int * sorted, int size, int * index);
void sortedarraymapfloat(float * unsorted, float * sorted, int size, int * index);
void sortedarraymapdouble(double * unsorted, double * sorted, int size, int * index);

void selectionsortuint(unsigned int * arr, int l);
void selectionsortint(int * arr, int l);
void selectionsortfloat(float * arr, int l);
void selectionsortdouble(double * arr, int l);

void bubblesortuint(unsigned int * arr, int l);
void bubblesortint(int * arr, int l);
void bubblesortfloat(float * arr, int l);
void bubblesortdouble(double * arr, int l);

void treesortuint(unsigned int * arr, int size);
void treesortint(int * arr, int size);
void treesortfloat(float * arr, int size);
void treesortdouble(double * arr, int size);

void heapsortuint(unsigned int * arr, int size);
void heapsortint(int * arr, int size);
void heapsortfloat(float * arr, int size);
void heapsortdouble(double * arr, int size);

void quicksortuint(unsigned int * arr, int size);
void quicksortint(int * arr, int size);
void quicksortfloat(float * arr, int size);
void quicksortdouble(double * arr, int size);

void mergesortuint(unsigned int * arr, unsigned int * scratch, int size);
void mergesortint(int * arr, int * scratch, int size);
void mergesortfloat(float * arr, float * scratch, int size);
void mergesortdouble(double * arr, double * scratch, int size);

void countingsortuint(unsigned int * arr, int size, unsigned int max_value);
void countingsortint(int * arr, int size, int max_value);

#endif /* SORT_H__ */
