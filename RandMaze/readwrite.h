#ifndef READ_WRITE_H
#define

void readMazeSize(int * size, int * width, char * filename);
void readMaze(unsigned short * map, int size, int width, char * filename);

void writeMaze(unsigned short * map, int size, int width);
void writeMazetest(unsigned short * map, int size, int width, char * filename);

#endif
