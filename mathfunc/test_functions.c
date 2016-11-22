#include <stdio.h>
#include <stdint.h>

#include "search.h"
#include "sort.h"
#include "array.h"

int main(){
	unsigned int arr[10] = {71, 6, 55, 32, 7, 10, 2, 5, 8, 99};
	int min_index[10], max_index[10];
	int imin = 0, imax = 0;
	unsigned int min_val, max_val;
	findminmaxuint(arr, 10, 0, min_index, &imin, max_index, &imax, &min_val, &max_val);
	printArrayUInt("max_index", imax, max_index);
	printArrayUInt("min_index", imin, min_index);
}