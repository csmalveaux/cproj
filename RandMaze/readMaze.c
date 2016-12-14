#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void readMaze(unsigned short * map, int * size, int * width, char * filename){
	FILE *fp;
	int height, wid;
	fp = fopen(filename, "r");
	if(fp == NULL) { printf("could not open file\n"); exit(0); }
	
	char buff[10];
	fscanf(fp, "%s", buff);
	if(strcmp(buff, "dim") == 0){
		fscanf(fp, "%d", &height);
		fscanf(fp, "%d", &wid);
		int sz = wid * height;

		* size = sz;
		* width = wid;
		
		map = malloc(sizeof(unsigned short)*sz);
		fread (map, sizeof(unsigned short), sz, fp);
	}
	else{printf("No dim tag\n"); exit(0);}
	
	fclose(fp);
}
