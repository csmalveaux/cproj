#include <stdio.h>
#include <time.h>

void writeMaze(unsigned short * map, int size, int width);
void writeMazetest(unsigned short * map, int size, int width, char * filename);

void writeMaze(unsigned short * map, int size, int width){
	FILE *fp;
	time_t timer;
	struct tm * currtime;
	char filename[20];
	
	time(&timer);
	currtime = gmtime ( &timer );
	sprintf(filename, "%d%d%d_%d%d%d.mz", currtime->tm_year, currtime->tm_mon, currtime->tm_mday, currtime->tm_hour, currtime->tm_min, currtime->tm_sec);
	
	fp = fopen(filename, "w");
	if(fp == NULL) exit(0);
	fwrite(&size, sizeof(int), 1, fp);
	fwrite(&width, sizeof(int), 1, fp);
	fwrite(map, sizeof(unsigned short), size, fp);
	fclose(fp);
}

void writeMazetest(unsigned short * map, int size, int width, char * filename){
	FILE *fp;
	time_t timer;
	struct tm * currtime;
	
	time(&timer);
	currtime = gmtime ( &timer );
	sprintf(filename, "%d%d%d_%d%d%d.mz", currtime->tm_year, currtime->tm_mon, currtime->tm_mday, currtime->tm_hour, currtime->tm_min, currtime->tm_sec);
	
	fp = fopen(filename, "w");
	if(fp == NULL) {printf("could not open file\n"); exit(0);}
	fwrite(&size, sizeof(int), 1, fp);
	fwrite(&width, sizeof(int), 1, fp);
	fwrite(map, sizeof(unsigned short), size, fp);
	fclose(fp);
}