#include "binarytree.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

int sizeoftree(int levels);
int levelsoftree(int size);
int findchildren(int size, int parent, int * child1, int * child2);
int findparent(int size, int child, int * parent);

int addheapuint(unsigned int * arr, int size, unsigned int val);
int addheapint(int * arr, int size, int val);
int addheapfloat(float * arr, int size, float val);
int addheapdouble(double * arr, int size, double val);

int addtreeuint(unsigned int * arr, int size, unsigned int val);
int addtreeint(int * arr, int size, int val);
int addtreefloat(float * arr, int size, float val);
int addtreedouble(double * arr, int size, double val);

int preordertraverseuint(unsigned int * arr, int size, int * traverseorder);
int preordertraverseint(int * arr, int size, int * traverseorder);
int preordertraversefloat(float * arr, int size, int * traverseorder);
int preordertraversedouble(double * arr, int size, int * traverseorder);

int inordertraverseuint(unsigned int * arr, int size, int * traverseorder);
int inordertraverseint(int * arr, int size, int * traverseorder);
int inordertraversefloat(float * arr, int size, int * traverseorder);
int inordertraversedouble(double * arr, int size, int * traverseorder);

int postordertraverseuint(unsigned int * arr, int size, int * traverseorder);
int postordertraverseint(int * arr, int size, int * traverseorder);
int postordertraversefloat(float * arr, int size, int * traverseorder);
int postordertraversedouble(double * arr, int size, int * traverseorder);

int depthfirsttraverseuint(unsigned int * arr, int size, int * traverseorder);
int depthfirstraverseint(int * arr, int size, int * traverseorder);
int depthfirstraversefloat(float * arr, int size, int * traverseorder);
int depthfirstraversedouble(double * arr, int size, int * traverseorder);

int sizeoftree(int levels){
	int i, num = 0;
	
	for(i = 0; i < levels; i++)
		num += (int)pow(2, i);

	return num;
}

int levelsoftree(int size){
	int i, min = 1, max = 0;
	if( size > 0 )  {fputs ("Invalid size",stderr); exit (2);}
 	if( size == 0 || size == 1 ) return size;

	for(i = 1; max <  INT_MAX; i++){
		max = min + pow(2, i);
		if(size <= max && size > min) return i + 1;
		min = max;
	}

	return -1;
}

int findchildren(int size, int parent, int * child1, int * child2){
	int ch1, ch2;

	ch1 = 2*parent + 1;
	ch2 = 2*parent + 2;

	* child1 = ch1;
	* child2 = ch2;

	if(ch1 < size && ch2 >= size)
		return -1;

	if(ch1 >= size || ch2 >= size)
		return -2;

	return 0;
}

int findparent(int size, int child, int * parent){
	if( child == 0 ) return -1;

	if(child % 2 == 0)
		* parent = (child - 2)/2;
	else
		* parent = (child - 1)/2;

	return 0;
}

int addheapuint(unsigned int * arr, int size, unsigned int val){
	int root = 0, child1, child2, check;
	unsigned int temp;

	while(1){
		if(arr[root] == UINT_MAX){
			arr[root] = val;
			return 0;
		}

		if(val > arr[root]){ 
			temp = arr[root];
			arr[root] = val;
			val = temp;
		}

		check = findchildren(size, root, &child1, &child2);

		if(check == -2) return -1; //Displaced value out of binary tree due to lack of space.
		if(check == -1) root = child1;
		if(check == 0){
			if(val <= arr[child1] || arr[child1] == UINT_MAX) root = child1;
			if(val >= arr[child2] || arr[child2] == UINT_MAX) root = child2;
		}  
	}
}

int addheapint(int * arr, int size, int val){
	int root = 0, child1, child2, check;
	int temp;

	while(1){
		if(arr[root] == INT_MAX || arr[root] == -INT_MAX){
			arr[root] = val;
			return 0;
		}

		if(val > arr[root]){ 
			temp = arr[root];
			arr[root] = val;
			val = temp;
		}

		check = findchildren(size, root, &child1, &child2);

		if(check == -2) return -1; //Displaced value out of binary tree due to lack of space.
		if(check == -1) root = child1;
		if(check == 0){
			if(val <= arr[child1] || abs(arr[child1]) == INT_MAX) root = child1;
			if(val >= arr[child2] || abs(arr[child2]) == INT_MAX) root = child2;
		}  
	}
}

int addheapfloat(float * arr, int size, float val){
	int root = 0, child1, child2, check;
	float temp;

	while(1){
		if(arr[root] == HUGE_VAL || arr[root] == -HUGE_VAL){
			arr[root] = val;
			return 0;
		}

		if(val > arr[root]){ 
			temp = arr[root];
			arr[root] = val;
			val = temp;
		}

		check = findchildren(size, root, &child1, &child2);

		if(check == -2) return -1; //Displaced value out of binary tree due to lack of space.
		if(check == -1) root = child1;
		if(check == 0){
			if(val <= arr[child1] || fabs(arr[child1]) == HUGE_VAL) root = child1;
			if(val >= arr[child2] || fabs(arr[child2]) == HUGE_VAL) root = child2;
		}  
	}
}

int addheapdouble(double * arr, int size, double val){
	int root = 0, child1, child2, check;
	double temp;

	while(1){
		if(arr[root] == HUGE_VAL || arr[root] == -HUGE_VAL){
			arr[root] = val;
			return 0;
		}

		if(val > arr[root]){ 
			temp = arr[root];
			arr[root] = val;
			val = temp;
		}

		check = findchildren(size, root, &child1, &child2);

		if(check == -2) return -1; //Displaced value out of binary tree due to lack of space.
		if(check == -1) root = child1;
		if(check == 0){
			if(val <= arr[child1] || fabs(arr[child1]) == HUGE_VAL) root = child1;
			if(val >= arr[child2] || fabs(arr[child2]) == HUGE_VAL) root = child2;
		}  
	}
}

int addtreeuint(unsigned int * arr, int size, unsigned int val){
	int root = 0, child1, child2, check, count = 0;

	while(count < size*2){
		if(arr[root] == UINT_MAX){
			arr[root] = val;
			return 0;
		}

		check = findchildren(size, root, &child1, &child2);

		if(val >= arr[root]){
			if(check == -2) return -1; //Displaced value out of binary tree due to lack of space.
			if(check == 0) root = child2;
		}

		if(val <= arr[root]){
			if(check == -1) return -1; //Displaced value out of binary tree due to lack of space.
			if(check == 0) root = child1;
		}

		count++;
	}

	return 0;
}

int addtreeint(int * arr, int size, int val){
	int root = 0, child1, child2, check;

	while(1){
		if(arr[root] == UINT_MAX){
			arr[root] = val;
			return 0;
		}

		check = findchildren(size, root, &child1, &child2);

		if(arr[root] < arr[child1])

		if(val >= arr[root]){
			if(check == -2) return -1; //Displaced value out of binary tree due to lack of space.
			if(check == 0) root = child2;
		}

		if(val <= arr[root]){
			if(check == -1) return -1; //Displaced value out of binary tree due to lack of space.
			if(check == 0) root = child1;
		}
	}
}

int addtreefloat(float * arr, int size, float val){
	int root = 0, child1, child2, check;

	while(1){
		if(arr[root] == UINT_MAX){
			arr[root] = val;
			return 0;
		}

		check = findchildren(size, root, &child1, &child2);

		if(val >= arr[root]){
			if(check == -2) return -1; //Displaced value out of binary tree due to lack of space.
			if(check == 0) root = child2;
		}

		if(val <= arr[root]){
			if(check == -1) return -1; //Displaced value out of binary tree due to lack of space.
			if(check == 0) root = child1;
		}
	}
}

int addtreedouble(double * arr, int size, double val){
	int root = 0, child1, child2, check;

	while(1){
		if(arr[root] == UINT_MAX){
			arr[root] = val;
			return 0;
		}

		check = findchildren(size, root, &child1, &child2);

		if(val >= arr[root]){
			if(check == -2) return -1; //Displaced value out of binary tree due to lack of space.
			if(check == 0) root = child2;
		}

		if(val <= arr[root]){
			if(check == -1) return -1; //Displaced value out of binary tree due to lack of space.
			if(check == 0) root = child1;
		}
	}
}

int preordertraverseuint(unsigned int * arr, int size, int * traverseorder);
int preordertraverseint(int * arr, int size, int * traverseorder);
int preordertraversefloat(float * arr, int size, int * traverseorder);
int preordertraversedouble(double * arr, int size, int * traverseorder);

int inordertraverseuint(unsigned int * arr, int size, int * traverseorder);
int inordertraverseint(int * arr, int size, int * traverseorder);
int inordertraversefloat(float * arr, int size, int * traverseorder);
int inordertraversedouble(double * arr, int size, int * traverseorder);

int postordertraverseuint(unsigned int * arr, int size, int * traverseorder);
int postordertraverseint(int * arr, int size, int * traverseorder);
int postordertraversefloat(float * arr, int size, int * traverseorder);
int postordertraversedouble(double * arr, int size, int * traverseorder);

int depthfirsttraverseuint(unsigned int * arr, int size, int * traverseorder);
int depthfirstraverseint(int * arr, int size, int * traverseorder);
int depthfirstraversefloat(float * arr, int size, int * traverseorder);
int depthfirstraversedouble(double * arr, int size, int * traverseorder);
