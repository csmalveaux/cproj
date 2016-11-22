#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "sort.h"
#include "search.h"
#include "binarytree.h"

void switchuint(unsigned int * arr, int size, int index1, int index2);
void switchint(int * arr, int size, int index1, int index2);
void switchfloat(float * arr, int size, int index1, int index2);
void switchdouble(double * arr, int size, int index1, int index2);

void sortedarraymapuint(unsigned int * unsorted, unsigned int * sorted, int size, int * index);
void sortedarraymapint(int * unsorted, int * sorted, int size, int * index);
void sortedarraymapfloat(float * unsorted, float * sorted, int size, int * index);
void sortedarraymapdouble(double * unsorted, double * sorted, int size, int * index);

void insertionsortuint(unsigned int * arr, int l);
void insertionsortint(int * arr, int l);
void insertionsortfloat(float * arr, int l);
void insertionsortdouble(double * arr, int l);

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


void switchuint(unsigned int * arr, int size, int index1, int index2){
    unsigned int temp1, temp2;

    //if(index1 > size - 1 || index2 > size - 1 || index1 < 0 || index2 < 0) {fputs ("Array index out of bounds",stderr); exit (2);}

    temp1 = arr[index1];
    temp2 = arr[index2];

    arr[index1] = temp2;
    arr[index2] = temp1;
}

void switchint(int * arr, int size, int index1, int index2){
    int temp1, temp2;

    //if(index1 > size - 1 || index2 > size - 1 || index1 < 0 || index2 < 0) {fputs ("Array index out of bounds",stderr); exit (2);}

    temp1 = arr[index1];
    temp2 = arr[index2];

    arr[index1] = temp2;
    arr[index2] = temp1;
}

void switchfloat(float * arr, int size, int index1, int index2){
    float temp1, temp2;

    //if(index1 > size - 1 || index2 > size - 1 || index1 < 0 || index2 < 0) {fputs ("Array index out of bounds",stderr); exit (2);}

    temp1 = arr[index1];
    temp2 = arr[index2];

    arr[index1] = temp2;
    arr[index2] = temp1;
}

void switchdouble(double * arr, int size, int index1, int index2){
    double temp1, temp2;

    //if(index1 >= size || index2 >= size || index1 < 0 || index2 < 0) {fputs ("Array index out of bounds",stderr); exit (2);}

    temp1 = arr[index1];
    temp2 = arr[index2];

    arr[index1] = temp2;
    arr[index2] = temp1;
}

void sortedarraymapuint(unsigned int * unsorted, unsigned int * sorted, int size, int * index){
	int i, j;
	unsigned int temp[size];
	memcpy(temp, sorted, sizeof(unsigned int)*size);

	for(i = 0; i < size; i++)
		for(j = 0; j < size; j++)
			if(unsorted[i] == temp[j]){
				temp[j] = UINT_MAX;
				index[i] = j;
				continue;
			}

}

void sortedarraymapint(int * unsorted, int * sorted, int size, int * index){
	int i, j;
	int temp[size];
	memcpy(temp, sorted, sizeof(int)*size);

	for(i = 0; i < size; i++)
		for(j = 0; j < size; j++)
			if(unsorted[i] == temp[j]){
				temp[j] = -INT_MAX;
				index[i] = j;
				break;
			}

}

void sortedarraymapfloat(float * unsorted, float * sorted, int size, int * index){
	int i, j;
	float temp[size];
	memcpy(temp, sorted, sizeof(float)*size);

	for(i = 0; i < size; i++)
		for(j = 0; j < size; j++)
			if(unsorted[i] == temp[j]){
				temp[j] = -HUGE_VAL;
				index[i] = j;
				continue;
			}

}

void sortedarraymapdouble(double * unsorted, double * sorted, int size, int * index){
	int i, j;
	double temp[size];
	memcpy(temp, sorted, sizeof(double)*size);

	for(i = 0; i < size; i++)
		for(j = 0; j < size; j++)
			if(unsorted[i] == temp[j]){
				temp[j] = -HUGE_VAL;
				index[i] = j;
				break;
			}

}

void insertionsortuint(unsigned int * arr, int l){
	int i, j;
	unsigned int temparr[l], temp, val;
	for(i = 0; i < l; i++){
		val = arr[i];
		for(j = 0; j <= i; j++){
			if(j == i) temparr[j] = val;
			if(val < temparr[j]){
				temp = temparr[j];
				temparr[j] = val;
				val = temp;
			}

		}
	}
}

void insertionsortint(int * arr, int l){
	int i, j;
	int temparr[l], temp, val;
	for(i = 0; i < l; i++){
		val = arr[i];
		for(j = 0; j <= i; j++){
			if(j == i) temparr[j] = val;
			if(val < temparr[j]){
				temp = temparr[j];
				temparr[j] = val;
				val = temp;
			}

		}
	}
}

void insertionsortfloat(float * arr, int l){
	int i, j;
	float temparr[l], temp, val;
	for(i = 0; i < l; i++){
		val = arr[i];
		for(j = 0; j <= i; j++){
			if(j == i) temparr[j] = val;
			if(val < temparr[j]){
				temp = temparr[j];
				temparr[j] = val;
				val = temp;
			}

		}
	}
}

void insertionsortdouble(double * arr, int l){
	int i, j;
	double temparr[l], temp, val;
	for(i = 0; i < l; i++){
		val = arr[i];
		for(j = 0; j <= i; j++){
			if(j == i) temparr[j] = val;
			if(val < temparr[j]){
				temp = temparr[j];
				temparr[j] = val;
				val = temp;
			}

		}
	}
}

void selectionsortuint(unsigned int * arr, int l){
	int i, ind = 0;
	unsigned int dummy;

	for(i = 0; i < l - i; i++){
		unsigned int temp[l - i];
		memcpy(temp, arr + i, sizeof(unsigned int)*(l - i));
		findminuint(temp, l - i, &dummy, &ind);
		if(0 != ind)
			switchuint(arr, l, i, ind + i);
	}

}

void selectionsortint(int * arr, int l){
	int i, ind = 0;
	int dummy;

	for(i = 0; i < l - i; i++){
		int temp[l - i];
		memcpy(temp, arr + i, sizeof(int)*(l - i));
		findminint(temp, l - i, &dummy, &ind);
		if(0 != ind)
			switchint(arr, l, i, ind + i);
	}

}

void selectionsortfloat(float * arr, int l){
	int i, ind = 0;
	float dummy;

	for(i = 0; i < l - i; i++){
		float temp[l - i];
		memcpy(temp, arr + i, sizeof(float)*(l - i));
		findminfloat(temp, l - i, &dummy, &ind);
		if(0 != ind)
			switchfloat(arr, l, i, ind + i);
	}

}

void selectionsortdouble(double * arr, int l){
	int i, ind;
	double dummy;

	for(i = 0; i < l - i; i++){
		double temp[l - i];
		memcpy(temp, arr + i, sizeof(double)*(l - i));

		findmindouble(temp, l - i, &dummy, &ind);
		if(0 != ind)
			switchdouble(arr, l, i, ind + i);
	}

}

void bubblesortuint(unsigned int * arr, int l){
	int i, order = 0;

	while(order != 1){
		order = 1;
		for(i = 0; i < l - 1; i++)
			if(arr[i] > arr[i + 1]){
				order = 0;
				switchuint(arr, l, i, i + 1);
			}
	}
}

void bubblesortint(int * arr, int l){
	int i, order = 0;

	while(order != 1){
		order = 1;
		for(i = 0; i < l - 1; i++)
			if(arr[i] > arr[i + 1]){
				order = 0;
				switchint(arr, l, i, i + 1);
			}
	}
}

void bubblesortfloat(float * arr, int l){
	int i, order = 0;

	while(order != 1){
		order = 1;
		for(i = 0; i < l - 1; i++)
			if(arr[i] > arr[i + 1]){
				order = 0;
				switchfloat(arr, l, i, i + 1);
			}
	}
}

void bubblesortdouble(double * arr, int l){
	int i, order = 0;

	while(order != 1){
		order = 1;
		for(i = 0; i < l - 1; i++)
			if(arr[i] > arr[i + 1]){
				order = 0;
				switchdouble(arr, l, i, i + 1);
			}
	}
}

void treesortuint(unsigned int * arr, int size){
	int i, root, parent, child1, child2, check, order = 0; 

	while(order != 0){
		order = 1;
		for(i = 0; i < size; i++){
			root = i;
			check = findchildren(size, root, &child1, &child2);

			// No children checking parent value
			if(check == -2){ 
				check = findparent(size, root, &parent);
				if(check == -1) return; // Only one element in binary tree; root. Already sorted by default.
				if(root % 2 == 0) if(arr[root] < arr[parent]){ switchuint(arr, size, root, parent); order = 0; }
				if(root % 2 == 1) if(arr[root] > arr[parent]){ switchuint(arr, size, root, parent); order = 0; }
			}

			// Only one child, checking only left hand child's value.
			if(check == -1) if(arr[child1] > arr[root]){ switchuint(arr, size, child1, root); order = 0; }

			// Has both left and right hand children, checking both children's values
			if(check == 0){
				if(arr[child1] > arr[root]){ switchuint(arr, size, child1, root); order = 0; }
				if(arr[child2] < arr[root]){ switchuint(arr, size, child2, root); order = 0; }
			}

		}
	}
}

void treesortint(int * arr, int size){
	int i, root, parent, child1, child2, check, order = 0; 

	while(order != 0){
		order = 1;
		for(i = 0; i < size; i++){
			root = i;
			check = findchildren(size, root, &child1, &child2);

			// No children checking parent value
			if(check == -2){ 
				check = findparent(size, root, &parent);
				if(check == -1) return; // Only one element in binary tree; root. Already sorted by default.
				if(root % 2 == 0) if(arr[root] < arr[parent]){ switchint(arr, size, root, parent); order = 0; }
				if(root % 2 == 1) if(arr[root] > arr[parent]){ switchint(arr, size, root, parent); order = 0; }
			}

			// Only one child, checking only left hand child's value.
			if(check == -1) if(arr[child1] > arr[root]){ switchint(arr, size, child1, root); order = 0; }

			// Has both left and right hand children, checking both children's values
			if(check == 0){
				if(arr[child1] > arr[root]){ switchint(arr, size, child1, root); order = 0; }
				if(arr[child2] < arr[root]){ switchint(arr, size, child2, root); order = 0; }
			}

		}
	}
}

void treesortfloat(float * arr, int size){
	int i, root, parent, child1, child2, check, order = 0; 

	while(order != 0){
		order = 1;
		for(i = 0; i < size; i++){
			root = i;
			check = findchildren(size, root, &child1, &child2);

			// No children checking parent value
			if(check == -2){ 
				check = findparent(size, root, &parent); 
				if(check == -1) return; // Only one element in binary tree; root. Already sorted by default.
				if(root % 2 == 0) if(arr[root] < arr[parent]){ switchfloat(arr, size, root, parent); order = 0; }
				if(root % 2 == 1) if(arr[root] > arr[parent]){ switchfloat(arr, size, root, parent); order = 0; }
			}

			// Only one child, checking only left hand child's value.
			if(check == -1) if(arr[child1] > arr[root]){ switchfloat(arr, size, child1, root); order = 0; }

			// Has both left and right hand children, checking both children's values
			if(check == 0){
				if(arr[child1] > arr[root]){ switchfloat(arr, size, child1, root); order = 0; }
				if(arr[child2] < arr[root]){ switchfloat(arr, size, child2, root); order = 0; }
			}

		}
	}
}

void treesortdouble(double * arr, int size){
	int i, root, parent, child1, child2, check, order = 0; 

	while(order != 0){
		order = 1;
		for(i = 0; i < size; i++){
			root = i;
			check = findchildren(size, root, &child1, &child2);

			// No children checking parent value
			if(check == -2){ 
				check = findparent(size, root, &parent);
				if(check == -1) return; // Only one element in binary tree; root. Already sorted by default.
				if(root % 2 == 0) if(arr[root] < arr[parent]){ switchdouble(arr, size, root, parent); order = 0; }
				if(root % 2 == 1) if(arr[root] > arr[parent]){ switchdouble(arr, size, root, parent); order = 0; }
			}

			// Only one child, checking only left hand child's value.
			if(check == -1) if(arr[child1] > arr[root]){ switchdouble(arr, size, child1, root); order = 0; }

			// Has both left and right hand children, checking both children's values
			if(check == 0){
				if(arr[child1] > arr[root]){ switchdouble(arr, size, child1, root); order = 0; }
				if(arr[child2] < arr[root]){ switchdouble(arr, size, child2, root); order = 0; }
			}

		}
	}
}

void heapsortuint(unsigned int * arr, int size){
	int i, root, parent, child1, child2, check, order = 0; 

	while(order != 0){
		order = 1;
		for(i = 0; i < size; i++){
			root = i;
			check = findchildren(size, root, &child1, &child2);

			// No children checking parent value
			if(check == -2){ 
				check = findparent(size, root, &parent);
				if(check == -1) return; // Only one element in binary tree; root. Already sorted by default.
				if(root % 2 == 0) if(arr[root] > arr[parent] || arr[parent] == UINT_MAX){ switchuint(arr, size, root, parent); order = 0; }
				if(root % 2 == 1) if(arr[root] > arr[parent] || arr[parent] == UINT_MAX){ switchuint(arr, size, root, parent); order = 0; }
			}

			// Only one child, checking only left hand child's value.
			if(check == -1) if(arr[child1] > arr[root]){ switchuint(arr, size, child1, root); order = 0; }

			// Has both left and right hand children, checking both children's values
			if(check == 0){
				if(arr[child1] > arr[root] && arr[child1] != UINT_MAX){ switchuint(arr, size, child1, root); order = 0; }
				if(arr[child2] > arr[root] && arr[child2] != UINT_MAX){ switchuint(arr, size, child2, root); order = 0; }
				if(arr[child1] > arr[child2] && arr[child2] != UINT_MAX){ switchuint(arr, size, child1, child2); order = 0; }
			}

		}
	}
}

void heapsortint(int * arr, int size){
	int i, root, parent, child1, child2, check, order = 0; 

	while(order != 0){
		order = 1;
		for(i = 0; i < size; i++){
			root = i;
			check = findchildren(size, root, &child1, &child2);

			// No children checking parent value
			if(check == -2){ 
				check = findparent(size, root, &parent);
				if(check == -1) return; // Only one element in binary tree; root. Already sorted by default.
				if(root % 2 == 0) if(arr[root] > arr[parent] || arr[parent] == INT_MAX || arr[parent] == -INT_MAX){ switchint(arr, size, root, parent); order = 0; }
				if(root % 2 == 1) if(arr[root] > arr[parent] || arr[parent] == INT_MAX || arr[parent] == -INT_MAX){ switchint(arr, size, root, parent); order = 0; }
			}

			// Only one child, checking only left hand child's value.
			if(check == -1) if(arr[child1] > arr[root]){ switchint(arr, size, child1, root); order = 0; }

			// Has both left and right hand children, checking both children's values
			if(check == 0){
				if(arr[child1] > arr[root] &&  arr[parent] != INT_MAX && arr[parent] != -INT_MAX){ switchint(arr, size, child1, root); order = 0; }
				if(arr[child2] > arr[root] &&  arr[parent] != INT_MAX && arr[parent] != -INT_MAX){ switchint(arr, size, child2, root); order = 0; }
				if(arr[child1] > arr[child2] &&  arr[parent] != INT_MAX && arr[parent] != -INT_MAX){ switchint(arr, size, child1, child2); order = 0; }
			}

		}
	}
}

void heapsortfloat(float * arr, int size){
	int i, root, parent, child1, child2, check, order = 0; 

	while(order != 0){
		order = 1;
		for(i = 0; i < size; i++){
			root = i;
			check = findchildren(size, root, &child1, &child2);

			// No children checking parent value
			if(check == -2){ 
				check = findparent(size, root, &parent);
				if(check == -1) return; // Only one element in binary tree; root. Already sorted by default.
				if(root % 2 == 0) if(arr[root] > arr[parent] || fabs(arr[parent]) == HUGE_VAL ){ switchfloat(arr, size, root, parent); order = 0; }
				if(root % 2 == 1) if(arr[root] > arr[parent] || fabs(arr[parent]) == HUGE_VAL ){ switchfloat(arr, size, root, parent); order = 0; }
			}

			// Only one child, checking only left hand child's value.
			if(check == -1) if(arr[child1] > arr[root]){ switchfloat(arr, size, child1, root); order = 0; }

			// Has both left and right hand children, checking both children's values
			if(check == 0){
				if(arr[child1] > arr[root] && fabs(arr[child1]) != HUGE_VAL){ switchfloat(arr, size, child1, root); order = 0; }
				if(arr[child2] > arr[root] && fabs(arr[child2]) != HUGE_VAL){ switchfloat(arr, size, child2, root); order = 0; }
				if(arr[child1] > arr[child2] && fabs(arr[child2]) != HUGE_VAL){ switchfloat(arr, size, child1, child2); order = 0; }
			}

		}
	}
}

void heapsortdouble(double * arr, int size){
	int i, root, parent, child1, child2, check, order = 0; 

	while(order != 0){
		order = 1;
		for(i = 0; i < size; i++){
			root = i;
			check = findchildren(size, root, &child1, &child2);

			// No children checking parent value
			if(check == -2){ 
				check = findparent(size, root, &parent);
				if(check == -1) return; // Only one element in binary tree; root. Already sorted by default.
				if(root % 2 == 0) if(arr[root] > arr[parent] || fabs(arr[parent]) == HUGE_VAL ){ switchdouble(arr, size, root, parent); order = 0; }
				if(root % 2 == 1) if(arr[root] > arr[parent] || fabs(arr[parent]) == HUGE_VAL ){ switchdouble(arr, size, root, parent); order = 0; }
			}

			// Only one child, checking only left hand child's value.
			if(check == -1) if(arr[child1] > arr[root]){ switchdouble(arr, size, child1, root); order = 0; }

			// Has both left and right hand children, checking both children's values
			if(check == 0){
				if(arr[child1] > arr[root] && fabs(arr[child1]) != HUGE_VAL){ switchdouble(arr, size, child1, root); order = 0; }
				if(arr[child2] > arr[root] && fabs(arr[child2]) != HUGE_VAL){ switchdouble(arr, size, child2, root); order = 0; }
				if(arr[child1] > arr[child2] && fabs(arr[child2]) != HUGE_VAL){ switchdouble(arr, size, child1, child2); order = 0; }
			}

		}
	}
}

void quicksortuint(unsigned int * arr, int size){
	int i, i1 = 0, i2 = 0;
	unsigned int divider, temp1[size], temp2[size];

	if(size == 1) return;
	
	divider = arr[size/2];
	
	for(i = 0; i < size; i++){
		if(arr[i] < divider) { temp1[i1] = arr[i]; i1++; }
		if(arr[i] >= divider) { temp2[i2] = arr[i]; i2++; }
	}

	quicksortuint(temp1, i1);
	quicksortuint(temp2, i2);

	memcpy(arr, temp1, sizeof(unsigned int) * i1);
	memcpy(arr + i1, temp2, sizeof(unsigned int) * i2);

}

void quicksortint(int * arr, int size){
	int i, i1 = 0, i2 = 0;
	int divider, temp1[size], temp2[size];

	if(size == 1) return;
	
	divider = arr[size/2];
	
	for(i = 0; i < size; i++){
		if(arr[i] < divider) { temp1[i1] = arr[i]; i1++; }
		if(arr[i] >= divider) { temp2[i2] = arr[i]; i2++; }
	}

	quicksortint(temp1, i1);
	quicksortint(temp2, i2);

	memcpy(arr, temp1, sizeof(int) * i1);
	memcpy(arr + i1, temp2, sizeof(int) * i2);

}

void quicksortfloat(float * arr, int size){
	int i, i1 = 0, i2 = 0;
	float divider, temp1[size], temp2[size];

	if(size == 1) return;
	
	divider = arr[size/2];
	
	for(i = 0; i < size; i++){
		if(arr[i] < divider) { temp1[i1] = arr[i]; i1++; }
		if(arr[i] >= divider) { temp2[i2] = arr[i]; i2++; }
	}

	quicksortfloat(temp1, i1);
	quicksortfloat(temp2, i2);

	memcpy(arr, temp1, sizeof(float) * i1);
	memcpy(arr + i1, temp2, sizeof(float) * i2);

}

void quicksortdouble(double * arr, int size){
	int i, i1 = 0, i2 = 0;
	double divider, temp1[size], temp2[size];

	if(size == 1) return;
	
	divider = arr[size/2];
	
	for(i = 0; i < size; i++){
		if(arr[i] < divider) { temp1[i1] = arr[i]; i1++; }
		if(arr[i] >= divider) { temp2[i2] = arr[i]; i2++; }
	}

	quicksortdouble(temp1, i1);
	quicksortdouble(temp2, i2);

	memcpy(arr, temp1, sizeof(double) * i1);
	memcpy(arr + i1, temp2, sizeof(double) * i2);

}

void mergesortuint(unsigned int * arr, unsigned int * scratch, int size){
	int i, s = 0, i1 = 0, i2 = 0, size1, size2;
	size1 = size/2;
	size2 = size - size/2;
	unsigned int temp1[size1], scratch1[size1];
	unsigned int temp2[size2], scratch2[size2];

	memcpy(temp1, arr, size1);
	memcpy(temp2, arr + size1, size2);

	mergesortuint(temp1, scratch1, sizeof(unsigned int) * size1);
	mergesortuint(temp2, scratch2, sizeof(unsigned int) * size2);

	if( size == 1 ) return;

	while(i1 < size1 && i2 < size2){
		if(temp1[i1] <= temp2[i2]) { scratch[s] = temp1[i1]; i1++; }
		else { scratch[s] = temp2[i2]; i2++; }
		s++;
	}

	for(i = s; i < size; i++){
		if(i1 < size1) { scratch[i] = temp1[i1]; i1++; }
		if(i2 < size2) { scratch[i] = temp2[i2]; i2++; }  
	}

	memcpy(arr, scratch, sizeof(unsigned int) * size);

}

void mergesortint(int * arr, int * scratch, int size){
	int i, s = 0, i1 = 0, i2 = 0, size1, size2;
	if( size == 1 ) return;
	size1 = size/2;
	size2 = size - size/2;
	int temp1[size1], scratch1[size1];
	int temp2[size2], scratch2[size2];

	memcpy(temp1, arr, sizeof(int) * size1);
	memcpy(temp2, arr + size1, sizeof(int) * size2);

	mergesortint(temp1, scratch1, size1);
	mergesortint(temp2, scratch2, size2);

	// if( size == 1 ) return;

	while(i1 < size1 && i2 < size2){
		if(temp1[i1] <= temp2[i2]) { scratch[s] = temp1[i1]; i1++; }
		else { scratch[s] = temp2[i2]; i2++; }
		s++;
	}

	for(i = s; i < size; i++){
		if(i1 < size1) { scratch[i] = temp1[i1]; i1++; }
		if(i2 < size2) { scratch[i] = temp2[i2]; i2++; }  
	}

	memcpy(arr, scratch, sizeof(int) * size);

}

void mergesortfloat(float * arr, float * scratch, int size){
	int i, s = 0, i1 = 0, i2 = 0, size1, size2;
	size1 = size/2;
	size2 = size - size/2;
	float temp1[size1], scratch1[size1];
	float temp2[size2], scratch2[size2];

	memcpy(temp1, arr, size1);
	memcpy(temp2, arr + size1, size2);

	mergesortfloat(temp1, scratch1, sizeof(float) * size1);
	mergesortfloat(temp2, scratch2, sizeof(float) * size2);

	if( size == 1 ) return;

	while(i1 < size1 && i2 < size2){
		if(temp1[i1] <= temp2[i2]) { scratch[s] = temp1[i1]; i1++; }
		else { scratch[s] = temp2[i2]; i2++; }
		s++;
	}

	for(i = s; i < size; i++){
		if(i1 < size1) { scratch[i] = temp1[i1]; i1++; }
		if(i2 < size2) { scratch[i] = temp2[i2]; i2++; }  
	}

	memcpy(arr, scratch, sizeof(float) * size);

}

void mergesortdouble(double * arr, double * scratch, int size){
	int i, s = 0, i1 = 0, i2 = 0, size1, size2;
	size1 = size/2;
	size2 = size - size/2;
	double temp1[size1], scratch1[size1];
	double temp2[size2], scratch2[size2];

	memcpy(temp1, arr, size1);
	memcpy(temp2, arr + size1, size2);

	mergesortdouble(temp1, scratch1, sizeof(double) * size1);
	mergesortdouble(temp2, scratch2, sizeof(double) * size2);

	if( size == 1 ) return;

	while(i1 < size1 && i2 < size2){
		if(temp1[i1] <= temp2[i2]) { scratch[s] = temp1[i1]; i1++; }
		else { scratch[s] = temp2[i2]; i2++; }
		s++;
	}

	for(i = s; i < size; i++){
		if(i1 < size1) { scratch[i] = temp1[i1]; i1++; }
		if(i2 < size2) { scratch[i] = temp2[i2]; i2++; }  
	}

	memcpy(arr, scratch, sizeof(double) * size);

}

void countingsortuint(unsigned int * arr, int size, unsigned int max_value){
	int i, j, index = 0;
	unsigned int * counts = calloc(max_value, sizeof(unsigned int));

	for(i = 0; i < size; i++)
		counts[arr[i]] += 1;

	for(i = 0; i <= max_value; i++)
		for(j = 0; j < counts[i]; j++){
			arr[index] = i;
			i++;
		}
}

void countingsortint(int * arr, int size, int max_value){
	int i, j, index = 0;
	unsigned int * counts = calloc(max_value, sizeof(unsigned int));

	for(i = 0; i < size; i++)
		counts[arr[i]] += 1;

	for(i = 0; i <= max_value; i++)
		for(j = 0; j < counts[i]; j++){
			arr[index] = i;
			i++;
		}
}
