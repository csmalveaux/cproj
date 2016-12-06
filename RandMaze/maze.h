#ifndef MAZE_H
#define MAZE_H

#define MAX_SEEDS		6

extern unsigned int curr_pos[MAX_SEEDS];
extern unsigned int corner[4];
extern unsigned int sleep_time;
extern int goal_pos;

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

void push(unsigned int * stack, int * ind, unsigned int curr);
void pop(unsigned int * stack, int * ind);

int randCell(unsigned short * map, const int size, const int width, const int pseudorandom, int disregard, unsigned int curr, unsigned int * rand_cell);
void removeWall(unsigned short * map, const int size, const int width, unsigned int curr, unsigned int next);

void recursiveBacktracker(unsigned short * map, const int size, const int width, int num, int perfect, int pseudorandom, int print);

void calculateWeightsforSet(unsigned short * map, const int size, const int width, int set_size, unsigned int * set, int * weights);
int getWeight(unsigned short * map, const int size, const int width, unsigned int cell);
void randomizedTraversal(unsigned short * map, const int size, const int width, int pseudorandom, int num, int print);
void randomizedPrim(unsigned short * map, const int size, const int width, int perfect, int pseudorandom, int print);

void wilson(unsigned short * map, const int size, const int width, int pseudorandom, int loops, int print);

void displayMaze(unsigned short * map, const int sz, const int w);
void refreshMaze(unsigned short * map, const int sz, const int w);

void clearCursors();

#endif //MAZE_H