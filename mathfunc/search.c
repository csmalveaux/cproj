#include <math.h>
#include <limits.h>
#include "search.h"
#include "binarytree.h"
#include <stdio.h>

// #define MIN(a,b) (((a)<(b))?(a):(b))
// #define MAX(a,b) (((a)>(b))?(a):(b))

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

int findnodeheapuint(unsigned int * arr, int l, unsigned int val, int * index);
int findnodeheapint(int * arr, int l, int val, int * index);
int findnodeheapfloat(float * arr, int l, float val, float tol, int * index);
int findnodeheapdouble(double * arr, int l, double val, double tol, int * index);

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

int linearsearchuint(unsigned int * arr, int l, unsigned int val, int * index){
	int i;

	for(i = 0; i < l; i++)
		if(arr[i] ==  val){
			* index = i;
			return 0;
		}

	return -1;
}

int linearsearchint(int * arr, int l, int val, int * index){
	int i;

	for(i = 0; i < l; i++)
		if(arr[i] ==  val){
			* index = i;
			return 0;
		}

	return -1;
}

int linearsearchfloat(float * arr, int l, float val, float tol, int * index){
	int i;

	for(i = 0; i < l; i++)
		if(fabs(arr[i] - val) <=  tol){
			* index = i;
			return 0;
		}

	return -1;
}

int linearsearchdouble(double * arr, int l, double val, double tol, int * index){
	int i;

	for(i = 0; i < l; i++)
		if(fabs(arr[i] - val) <=  tol){
			* index = i;
			return 0;
		}

	return -1;
}

int binarysearchuint(unsigned int * arr, int l, unsigned int val, int * index){
	int min = 0, max = l - 1, mid;

	mid = (min + max)/2;

	while(min <= max){
		
		if(val == arr[mid]) { * index = mid; return 0; }
		if(val < arr[mid]) max = mid - 1;
		else if(val > arr[mid]) min = mid + 1;

		mid = (min + max)/2;
	}

	return -1;

}

int binarysearchint(int * arr, int l, int val, int * index){
	int min = 0, max = l - 1, mid;

	while(min < max){
		mid = (min + max)/2;
		if(val == arr[mid]) { * index = mid; return 0; }
		if(val < arr[mid]) max = mid + 1;
		else if(val > arr[mid]) min = mid - 1;
	}

	return -1;

}

int binarysearchfloat(float * arr, int l, float val, float tol, int * index){
	int min = 0, max = l - 1, mid;

	while(min < max){
		mid = (min + max)/2;
		if(fabs(arr[mid] - val) <=  tol) { * index = mid; return 0; }
		if(val < arr[mid]) max = mid + 1;
		else if(val > arr[mid]) min = mid - 1;
	}

	return -1;

}

int binarysearchdouble(double * arr, int l, double val, double tol, int * index){
	int min = 0, max = l - 1, mid;

	while(min < max){
		mid = (min + max)/2;
		if(fabs(arr[mid] - val) <=  tol) { * index = mid; return 0; }
		if(val < arr[mid]) max = mid + 1;
		else if(val > arr[mid]) min = mid - 1;
	}

	return -1;

}

int interpolationsearchuint(unsigned int * arr, int l, unsigned int val, int * index){
	int min, max, mid;

	while (min < max){
		mid = min + (max - min)*(val - arr[min])/(arr[max] - arr[min]);
		if(val == arr[mid]) { * index = mid; return 0; }
		if(val < arr[mid]) max = mid;
		else if(val > arr[mid]) min = mid;
	}

	return -1;

}

int interpolationsearchint(int * arr, int l, int val, int * index){
	int min, max, mid;

	while (min < max){
		mid = min + (max - min)*(val - arr[min])/(arr[max] - arr[min]);
		if(val == arr[mid]) { * index = mid; return 0; }
		if(val < arr[mid]) max = mid;
		else if(val > arr[mid]) min = mid;
	}

	return -1;

}

int interpolationsearchfloat(float * arr, int l, float val, float tol, int * index){
	int min = 0, max = l - 1, mid;

	while(min < max){
		mid = min + (max - min)*(val - arr[min])/(arr[max] - arr[min]);
		if(fabs(arr[mid] - val) <=  tol) { * index = mid; return 0; }
		if(val < arr[mid]) max = mid;
		else if(val > arr[mid]) min = mid;
	}

	return -1;

}

int interpolationsearchdouble(double * arr, int l, double val, double tol, int * index){
	int min = 0, max = l - 1, mid;

	while(min < max){
		mid = min + (max - min)*(val - arr[min])/(arr[max] - arr[min]);
		if(fabs(arr[mid] - val) <=  tol) { * index = mid; return 0; }
		if(val < arr[mid]) max = mid;
		else if(val > arr[mid]) min = mid;
	}

	return -1;

}

int findnodeheapuint(unsigned int * arr, int l, unsigned int val, int * index){
	int root = 0, child1, child2, check;

	while(1){
		if(arr[root] == val) { * index = root; return 0; }

		check = findchildren(l, root, &child1, &child2);
		
		if(arr[root] > val){
			if(check != -1) { root = child1; continue; }
			else 
				return -1;
		}
		else{
			if(check != -2) { root = child2; continue; }
			else 
				return -1;
		}
	} 
}

int findnodeheapint(int * arr, int l, int val, int * index){
	int root = 0, child1, child2, check;

	while(1){
		if(arr[root] == val) { * index = root; return 0; }

		check = findchildren(l, root, &child1, &child2);
		
		if(arr[root] > val){
			if(check != -1) { root = child1; continue; }
			else 
				return -1;
		}
		else{
			if(check != -2) { root = child2; continue; }
			else 
				return -1;
		}
	} 
}

int findnodeheapfloat(float * arr, int l, float val, float tol, int * index){
	int root = 0, child1, child2, check;

	while(1){
		if(fabs(arr[root] - val) < tol) { * index = root; return 0; }

		check = findchildren(l, root, &child1, &child2);
		
		if(arr[root] > val){
			if(check != -1) { root = child1; continue; }
			else 
				return -1;
		}
		else{
			if(check != -2) { root = child2; continue; }
			else 
				return -1;
		}
	} 
}

int findnodeheapdouble(double * arr, int l, double val, double tol, int * index){
	int root = 0, child1, child2, check;

	while(1){
		if(fabs(arr[root] - val) < tol) { * index = root; return 0; }

		check = findchildren(l, root, &child1, &child2);
		
		if(arr[root] > val){
			if(check != -1) { root = child1; continue; }
			else 
				return -1;
		}
		else{
			if(check != -2) { root = child2; continue; }
			else 
				return -1;
		}
	} 
}

int findminuint(unsigned int * arr, int l, unsigned int * val, int * index){
	unsigned int min = UINT_MAX;
	int i, ind;

	for(i = 0; i < l; i++)
		if(arr[i] < min){
			min = arr[i];
			ind = i;
		}

	* index = ind;
	* val = min;

    if(min == UINT_MAX)
        return -1;

    return 0;
}

int findminint(int * arr, int l, int * val, int * index){
	int min = INT_MAX;
	int i, ind;

	for(i = 0; i < l; i++)
		if(arr[i] < min){
			min = arr[i];
			ind = i;
		}

	* index = ind;
	* val = min;

    if(min == INT_MAX)
        return -1;

    return 0;
}

int findminfloat(float * arr, int l, float * val, int * index){
	float min = HUGE_VAL;
	int i, ind;

	for(i = 0; i < l; i++)
		if(arr[i] < min){
			min = arr[i];
			ind = i;
		}

	* index = ind;
	* val = min;

    if(min == HUGE_VAL)
        return -1;

    return 0;
}

int findmindouble(double * arr, int l, double * val, int * index){
	double min = HUGE_VAL;
	int i, ind;

	for(i = 0; i < l; i++)
		if(arr[i] < min){
			min = arr[i];
			ind = i;
		}

	* index = ind;
	* val = min;

    if(min == HUGE_VAL)
        return -1;

    return 0;
}

int findmaxuint(unsigned int * arr, int l, unsigned int * val, int * index){
	unsigned int max = 0;
	int i, ind;

	for(i = 0; i < l; i++)
		if(arr[i] > max){
			max = arr[i];
			ind = i;
		}

	* index = ind;
	* val = max;

    if(max == 0)
        return -1;

    return 0;
}

int findmaxint(int * arr, int l, int * val, int * index){
	int max = INT_MIN;
	int i, ind;

	for(i = 0; i < l; i++)
		if(arr[i] > max){
			max = arr[i];
			ind = i;
		}

	* index = ind;
	* val = max;

    if(max == INT_MIN)
        return -1;

    return 0;
}

int findmaxfloat(float * arr, int l, float * val, int * index){
	float max = -HUGE_VAL;
	int i, ind;

	for(i = 0; i < l; i++)
		if(arr[i] > max){
			max = arr[i];
			ind = i;
		}

	* index = ind;
	* val = max;

    if(max == -HUGE_VAL)
        return -1;

    return 0;
}

int findmaxdouble(double * arr, int l, double * val, int * index){
	double max = -HUGE_VAL;
	int i, ind;


	for(i = 0; i < l; i++)
		if(arr[i] > max){
			max = arr[i];
			ind = i;
		}
	printf("does it shit out here?");
	* index = ind;

	printf("does it shit out here?");
	* val = max;

    if(max == -HUGE_VAL)
        return -1;

    return 0;
}
