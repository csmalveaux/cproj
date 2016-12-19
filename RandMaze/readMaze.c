#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void readMazeSize(int * size, int * width, char * filename);
void readMaze(unsigned short * map, int size, int width, char * filename);

void readMazeSize(int * size, int * width, char * filename){

	FILE *fp;
	int sz;
	fp = fopen(filename, "r");
	if(fp == NULL) { printf("could not open file\n"); exit(0); }

	fread(size, sizeof(int), 1, fp);
	fread(width, sizeof(int), 1, fp);

	fclose(fp);
}

void readMaze(unsigned short * map, int size, int width, char * filename){
	FILE *fp;
	int sz, wid;
	fp = fopen(filename, "r");
	if(fp == NULL) { printf("could not open file\n"); exit(0); }
	
	fread(&sz, sizeof(int), 1, fp);
	fread(&wid, sizeof(int), 1, fp);

	if(sz != size || wid != width){ printf("Wrong maze file"); exit(0); }
	
	fread (map, sizeof(unsigned short), sz, fp);
	
	fclose(fp);
}
