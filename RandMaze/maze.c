#include <time.h>  
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "maze.h"
#include "map.h"

unsigned int curr_pos[MAX_SEEDS];
unsigned int corner[4];
unsigned int sleep_time = 100000000;
int goal_pos = -1;

// Prototypes
void getCorners(const int size, const int width);

int getNCell(const int c, const int s, const int w);
int getSCell(const int c, const int s, const int w);
int getECell(const int c, const int s, const int w);
int getWCell(const int c, const int s, const int w);

int setNCell(const int c, const int sz, const int w, unsigned short * map);
int setSCell(const int c, const int sz, const int w, unsigned short * map);
int setECell(const int c, const int sz, const int w, unsigned short * map);
int setWCell(const int c, const int sz, const int w, unsigned short * map);

void removeNWall(const int c, const int sz, const int w, unsigned short * map);
void removeSWall(const int c, const int sz, const int w, unsigned short * map);
void removeEWall(const int c, const int sz, const int w, unsigned short * map);
void removeWWall(const int c, const int sz, const int w, unsigned short * map);

void getAdjacentCells(const int sz, const int w, const int c, int * n, int * s, int * e, int * ws);
void setAdjacentCells(unsigned short * map, const int sz, const int w, const int c);

void setAllWalls(unsigned short * map, const int sz);
void setCenterGoal(unsigned short * map, const int size, const int width);

void getCorners(const int size, const int width){
	corner[0] = 0;
	corner[1] = size - 1;
	corner[2] = width - 1;
	corner[3] = size -  width;
}

int getNCell(const int c, const int s, const int w){
	if((c - w) < 0)
		return -1;
	return c - w;
}

int getSCell(const int c, const int s, const int w){
	if((c + w) > s)
		return -1;
	return c + w;
}

int getECell(const int c, const int s, const int w){
	if((c + 1) > s)
		return -1;
	if(c % w == w - 1)
		return -1;
	return c + 1;
}

int getWCell(const int c, const int s, const int w){
	if((c - 1) < 0)
		return -1;
	if(c % w == 0)
		return -1;
	return c - 1;
}

int setNCell(const int c, const int sz, const int w, unsigned short * map){
	int n = getNCell(c, sz, w);
	if(n != -1){ map[n] |= SOUTH; return 0; }
	return -1;
}

int setSCell(const int c, const int sz, const int w, unsigned short * map){
	int s = getSCell(c, sz, w);
	if(s != -1){ map[s] |= NORTH; return 0; }
	return -1;
}

int setECell(const int c, const int sz, const int w, unsigned short * map){
	int e = getECell(c, sz, w);
	if(e != -1){ map[e] |= WEST; return 0; }
	return -1;
}

int setWCell(const int c, const int sz, const int w, unsigned short * map){
	int ws = getWCell(c, sz, w);
	if(ws != -1){ map[ws] |= EAST; return 0; }
	return -1;
}

void removeNWall(const int c, const int sz, const int w, unsigned short * map){
	int n = getNCell(c, sz, w);
	if(n < 0) return;
	if((map[c] & NORTH) == NORTH) map[c] ^= NORTH;
	if((map[n] & SOUTH) == SOUTH) map[n] ^= SOUTH;
}

void removeSWall(const int c, const int sz, const int w, unsigned short * map){
	int s = getSCell(c, sz, w);
	if(s < 0) return;
	if((map[c] & SOUTH) == SOUTH) map[c] ^= SOUTH;
	if((map[s] & NORTH) == NORTH) map[s] ^= NORTH;
}

void removeEWall(const int c, const int sz, const int w, unsigned short * map){
	int e = getECell(c, sz, w);
	if(e < 0) return;
	if((map[c] & EAST) == EAST) map[c] ^= EAST;
	if((map[e] & WEST) == WEST) map[e] ^= WEST;
}

void removeWWall(const int c, const int sz, const int w, unsigned short * map){
	int ws = getWCell(c, sz, w);
	if(ws < 0) return;
	if((map[c] & WEST) == WEST) map[c] ^= WEST;
	if((map[ws] & EAST) == EAST) map[ws] ^= EAST;
}

void getAdjacentCells(const int sz, const int w, int c, int * n, int * s, int * e, int * ws){
	* n = getNCell(c, sz, w);
	* s = getSCell(c, sz, w);
	* e = getECell(c, sz, w);
	* ws = getWCell(c, sz, w);
}

void setAdjacentCells(unsigned short * map, const int sz, const int w, const int c){
	unsigned int cell = map[c];

	if((cell & VISIT) ==  VISIT) return;

	if((cell & NORTH) ==  NORTH) setNCell(c, sz, w, map);
	if((cell & SOUTH) ==  SOUTH) setSCell(c, sz, w, map);
	if((cell & EAST) ==  EAST)   setECell(c, sz, w, map);
	if((cell & WEST) ==  WEST)   setWCell(c, sz, w, map);

	map[c] |= VISIT;
}

void setAllWalls(unsigned short * map, int sz){
	int i;
	for(i = 0; i < sz; i++){
		map[i] = 0x00;
		map[i] |= NORTH + SOUTH + EAST + WEST + VISIT;
	}

}

void setCenterGoal(unsigned short * map, const int size, const int width){
	int entry, wall, height = size/width;
	unsigned int center[4];

	center[0] = (height/2 - 1) * width + width/2 - 1;
	center[1] = (height/2 - 1) * width + width/2;
	center[2] = (height/2) * width + width/2 - 1;
	center[3] = (height/2) * width + width/2;

	map[center[0]] ^= (VISIT + EAST + SOUTH);
	map[center[1]] ^= (VISIT + WEST + SOUTH);
	map[center[2]] ^= (VISIT + EAST + NORTH);
	map[center[3]] ^= (VISIT + WEST + NORTH);

	srand (time(NULL));
	entry = rand() % 4;

	map[center[entry]] |= VISIT;

}

