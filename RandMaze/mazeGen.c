#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>  
#include <unistd.h>
#include <math.h>
#include "maze.h"
#include "map.h"
#include "sort.h"
#include "search.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

struct thread{
	int istack1, istack2; 
	unsigned int prev, curr, next;
	unsigned int * stack1;
	int * stack2;
	int finish, diff, momentum, threshold;
};

void initalize_thread(struct thread * thrd, unsigned int start, int perfect, int pseudorandom, const int width);

void push(unsigned int * stack, int * ind, unsigned int curr);
void pop(unsigned int * stack, int * ind);
void delete(unsigned int * stack, int * size, int index);

int randCell(unsigned short * map, const int size, const int width, int pseudorandom, int disregard, unsigned int curr, unsigned int * rand_cell);
void removeWall(unsigned short * map, const int size, const int width, unsigned int curr, unsigned int next);

void recursiveBacktracker(unsigned short * map, const int size, const int width, int num, int perfect, int pseudorandom, int print);

void calculateWeightsforSet(unsigned short * map, const int size, const int width, int set_size, unsigned int * set, int * weights);
int getWeight(unsigned short * map, const int size, const int width, unsigned int cell);
void randomizedTraversal(unsigned short * map, const int size, const int width, int pseudorandom, int num, int print);

int direction(struct thread * thrd, const int width);
int checkCell(unsigned int * set, int setsize, const int size, const int width, unsigned int cell, unsigned int * adj);

int addRandCells(unsigned short * map, const int size, const int width, unsigned int * stack, int * istack, int * stack2, unsigned int curr);
void removeVisited(unsigned short * map, unsigned int * stack, int * istack);
void randomizedPrim(unsigned short * map, const int size, const int width, int perfect, int pseudorandom, int print);

void wilson(unsigned short * map, const int size, const int width, int pseudorandom, int loops, int print);

void initalize_thread(struct thread * thrd, unsigned int start, int perfect, int pseudorandom, const int width){
	thrd->curr = start;
	thrd->istack1 = 0;
	thrd->istack2 = 0;
	thrd->finish = 0;
	thrd->diff = 0;
	thrd->momentum = 0;

	if(perfect == 0){
		srand (time(NULL));
		if(pseudorandom)
			thrd->threshold = rand() % width/8 + 1; //randomly select new threshold
		else 
			thrd->threshold = rand() % 6 + 1; //randomly select new threshold
	}

 }

void push(unsigned int * stack, int * ind, unsigned int curr){
	stack[* ind] = curr;
	* ind = * ind + 1;
}

void pop ( unsigned int * stack, int * ind){
	stack[* ind] = 0;
	* ind = * ind - 1;
}

void delete(unsigned int * stack, int * size, int index){
	memcpy(stack + index, stack + (index + 1), (* size - index) * sizeof(unsigned int));
	* size = * size - 1;
}

int randCell(unsigned short * map, const int size, const int width, int pseudorandom, int disregard, unsigned int curr, unsigned int * rand_cell){
	int rnum, north, south, east, west;
	int ind = 0, available[4];

	//Getting all of the adjacent cells
	getAdjacentCells(size, width, curr, &north, &south, &east, &west);

	//Checking if the adjacent cells exist and if they have been visisted
	if(north > -1){
		if(disregard) { available[ind] = north; ind++; }
		else if((map[north] & VISIT) == VISIT) { available[ind] = north; ind++; }
	}
	if(south > -1){
		if(disregard) { available[ind] = south; ind++; }
		else if((map[south] & VISIT) == VISIT) { available[ind] = south; ind++; }
	}
	if(east > -1){
		if(disregard) { available[ind] =  east; ind++; }
		else if((map[east]  & VISIT) == VISIT) { available[ind] =  east; ind++; }
	}
	if(west > -1){
		if(disregard) { available[ind] =  west; ind++; }
		else if((map[west]  & VISIT) == VISIT) { available[ind] =  west; ind++; }
	}


	//Checking if there are no available and return 
	if(ind == 0) return 0;

	if(ind == 1) {* rand_cell = available[0]; return 1;}

	if(pseudorandom){
		srand (time(NULL));
	}

	rnum = rand() % ind;

	//Selecting random cell from available neighbor options
	* rand_cell = available[rnum];

	return ind;
}

void removeWall(unsigned short * map, const int size, const int width, unsigned int curr, unsigned int next){
	int north, south, east, west;
	getAdjacentCells(size, width, curr, &north, &south, &east, &west);

	if(next == north) removeNWall(curr, size, width, map);
	else if(next == south) removeSWall(curr, size, width, map);
	else if(next == east) removeEWall(curr, size, width, map);
	else if(next == west) removeWWall(curr, size, width, map);

}

int direction(struct thread * thrd, const int width){
	int df, mo;
	df = thrd->curr - thrd->next;
	mo = thrd->momentum;

	//checking the direction of the last two moves in comparison to the direction of the next move
	if(df == thrd->diff){//same direction, add to momentum
		mo++;
		thrd->momentum = mo;
		return 0;
	}
	else{//different direction
		thrd->diff = df;
		if(mo >= thrd->threshold){//check if enough momentum from the previous direction exceeds threshold
			thrd->momentum = 0;
			usleep(sleep_time/10);
			// srand (time(NULL));
			thrd->threshold = rand() % width/8 + 1; //randomly select new threshold
			return 1;
		} 
		thrd->momentum = 0;
	}
	return 0;
}

void recursiveBacktracker(unsigned short * map, const int size, const int width, int num, int perfect, int pseudorandom, int print){
	int i, cell, done = 1;
	struct thread thrd[num];

	setAllWalls(map, size);
	setCenterGoal(map, size, width);
	srand (time(NULL));

	for(i = 0; i < num; i++){
		if(i < 4)
			initalize_thread(&thrd[i], corner[i], perfect, pseudorandom, MIN(width, size/width));
		else{
			
			usleep(sleep_time);
			
			cell = rand() % (size - i) + i;
			initalize_thread(&thrd[i], cell, perfect, pseudorandom, MIN(width, size/width));
		}
		thrd[i].stack1 = calloc(size, sizeof(unsigned int));
		map[thrd[i].curr] ^= VISIT;
		curr_pos[i] = thrd[i].curr;
	}

	if(print)
		displayMaze(map, size, width);
	
	do{
		done  = 1;
		for(i = 0; i < num; i++){
			if(thrd[i].finish == 0){
				if(randCell(map, size, width, pseudorandom, 0, thrd[i].curr, &thrd[i].next)){
					push(thrd[i].stack1, &thrd[i].istack1, thrd[i].curr);
					removeWall(map, size, width, thrd[i].curr, thrd[i].next);
					if(perfect == 0){ 
						if(direction(&thrd[i], MIN(width, size/width))){
							map[thrd[i].curr] |= VISIT;
						}
					}
					thrd[i].curr = thrd[i].next;
					curr_pos[i] = thrd[i].curr;
					map[thrd[i].curr] ^= VISIT;
				}
				else{
					pop(thrd[i].stack1, &thrd[i].istack1);
					thrd[i].curr = thrd[i].stack1[thrd[i].istack1];
					curr_pos[i] = thrd[i].curr;
				}

				if(thrd[i].istack1 == 0){
					thrd[i].finish = 1;
					free(thrd[i].stack1);
				}
				else
					done = 0;
			}
		}
		
		if(print)
			refreshMaze(map, size, width);

	}while(done == 0);

	if(print == 0)
		displayMaze(map, size, width);

}

void calculateWeightsforSet(unsigned short * map, const int size, const int width, int set_size, unsigned int * set, int * weights){
	int i;
	for(i = 0; i < set_size; i++){
		// Recalcualtes only for cells that aren't already confirmed to be surrounded by visted neighbors or edges.
		if(weights[i] > 0)
			weights[i] = getWeight(map, size, width, set[i]);
	}
}

int getWeight(unsigned short * map, const int size, const int width, unsigned int cell){
	int north, south, east, west;
	int weight = 0;

	getAdjacentCells(size, width, cell, &north, &south, &east, &west);

	//Checking if the adjacent cells exist and if they have been visisted
	if((north > -1) && ((map[north] & VISIT) == VISIT)) weight += (int) map[north];
	else weight--;
	if((south > -1) && ((map[south] & VISIT) == VISIT)) weight += (int) map[south];
	else weight--;
	if((east > -1)  && ((map[east]  & VISIT) == VISIT)) weight += (int) map[east];
	else weight--;
	if((west > -1)  && ((map[west]  & VISIT) == VISIT)) weight += (int) map[west];
	else weight--;

	return weight;

}

int checkCell(unsigned int * set, int setsize, const int size, const int width, unsigned int cell, unsigned int * adj){
	int i, rnum, check, north, south, east, west;
	int ind = 0, available[4];

	unsigned int tempset[setsize];
	mergesortint(set, tempset, setsize);

	//Getting all of the adjacent cells
	getAdjacentCells(size, width, cell, &north, &south, &east, &west);

	//Checking if the adjacent cells exist and if they have been visisted
	if(binarysearchuint(set, setsize, north, &check) == -1) { available[ind] = north; ind++; }
	if(binarysearchuint(set, setsize, south, &check) == -1) { available[ind] = south; ind++; }
	if(binarysearchuint(set, setsize, east, &check) == -1) { available[ind] =  east; ind++; }
	if(binarysearchuint(set, setsize, west, &check) == -1) { available[ind] =  west; ind++; }

	//Checking if there are no available and return 
	if(ind == 0) return 0;

	if(ind == 1) {* adj = available[0]; return 1;}

	// srand (time(NULL));
	rnum = rand() % ind;

	//Selecting random cell from available neighbor options
	* adj = available[rnum];

	return 1;
}

void randomizedTraversal(unsigned short * map, const int size, const int width, int pseudorandom, int num, int print){
	int i, j, cell, rnid, range, shift, weight, done = 0;
	struct thread thrd[num];

	setAllWalls(map, size);
	setCenterGoal(map, size, width);

	for(i = 0; i < num; i++){
		if(i < 4)
			initalize_thread(&thrd[i], corner[i], 0, pseudorandom, width);
		else{
			
			usleep(sleep_time);
			srand (time(NULL));
			cell = rand() % (size - i) + i;
			initalize_thread(&thrd[i], cell, 0, pseudorandom, width);
		}
		thrd[i].stack1 = calloc(size, sizeof(unsigned int));
		thrd[i].stack2 = calloc(size, sizeof(int));
		map[thrd[i].curr] ^= VISIT;
		curr_pos[i] = thrd[i].curr;

		push(thrd[i].stack1, &thrd[i].istack1, thrd[i].curr);
		weight = getWeight(map, size, width, thrd[i].curr);
		push(thrd[i].stack2, &thrd[i].istack2, weight);
	}

	if(print)
		displayMaze(map, size, width);

	srand (time(NULL));
	while(done == 0){
		done = 1;
		for(i = 0; i < num; i++){

			if(thrd[i].finish == 0){
				int tempweight[thrd[i].istack2];
				int temp[thrd[i].istack2];
				int index[thrd[i].istack2];

				calculateWeightsforSet(map, size, width, thrd[i].istack2, thrd[i].stack1, thrd[i].stack2);
				memcpy(tempweight, thrd[i].stack2, sizeof(int) * thrd[i].istack2);

				mergesortint(tempweight, temp, thrd[i].istack2);
				sortedarraymapint(tempweight, thrd[i].stack2, thrd[i].istack2, index);

				if(tempweight[thrd[i].istack2 - 1] == -4){
					thrd[i].finish = 1;
					// free(thrd[i].stack1);
					free(thrd[i].stack2);
				}
				else{
					done = 0;
					if((thrd[i].istack2 > 1) && (thrd[i].istack2 < size)){
						for(j = 0; j < thrd[i].istack2; j++){
							if(tempweight[j] > 0){
								shift = j;
								break;
							}
						}

						range = thrd[i].istack2 - shift;

						//Select random cell
						rnid = rand() % range;
					}
					else 
						{ rnid = thrd[i].istack2 - 1; shift = 0; }

					thrd[i].curr = thrd[i].stack1[index[rnid + shift]];

					if(randCell(map, size, width, pseudorandom, 0, thrd[i].curr, &thrd[i].next) > 0){
						push(thrd[i].stack1, &thrd[i].istack1, thrd[i].next);
						removeWall(map, size, width, thrd[i].curr, thrd[i].next);
						weight = getWeight(map, size, width, thrd[i].next);
						push(thrd[i].stack2, &thrd[i].istack2, weight);
						map[thrd[i].next] ^= VISIT;

						curr_pos[i] = thrd[i].next;

						if(print)
							refreshMaze(map, size, width);
					}
				}
			}
		}

	}

	for(i = 0; i < num; i++){
		int nsamples = (int)(thrd[i].istack1/log(thrd[i].istack1));
		int cell, adj;
		for(j = 0; j < nsamples; j++){
			cell = thrd[i].stack1[rand() % nsamples + (thrd[i].istack1 - nsamples - 1)];

			if(checkCell(thrd[i].stack1, thrd[i].istack1, size, width, cell, &adj)){
				removeWall(map, size, width, cell, adj);
				if(print)
					refreshMaze(map, size, width);
			}

		}

		free(thrd[i].stack1);

	}


	if(print == 0)
		displayMaze(map, size, width);

}

int addRandCells(unsigned short * map, const int size, const int width, unsigned int * stack, int * istack, int * stack2, unsigned int curr){
	int i, rnum, north, south, east, west;
	int ind = 0, available[4];

	//Getting all of the adjacent cells
	getAdjacentCells(size, width, curr, &north, &south, &east, &west);

	//Checking if the adjacent cells exist and if they have been visisted
	if(north > -1) if((map[north] & VISIT) == VISIT) { available[ind] = north; ind++; }
	if(south > -1) if((map[south] & VISIT) == VISIT) { available[ind] = south; ind++; }
	if(east > -1)  if((map[east]  & VISIT) == VISIT) { available[ind] =  east; ind++; }
	if(west > -1)  if((map[west]  & VISIT) == VISIT) { available[ind] =  west; ind++; }

	//Checking if there are no available and return 
	if(ind == 0) return 0;

	for(i = 0; i < ind; i++){
		push(stack, istack, available[i]);
		stack2[available[i]] = curr;
	}

	//Selecting random cell from available neighbor options
	

	return 1;
}

void removeVisited(unsigned short * map, unsigned int * stack, int * istack){
	int i;

	for(i = 0; i < * istack; i++)
		if((map[stack[i]] & VISIT) == 0)
			delete(stack, istack, i);
}

void randomizedPrim(unsigned short * map, const int size, const int width, int perfect, int pseudorandom, int print){
	int i, ncells, rnid, ist, done  = 0;
	struct thread thrd;

	setAllWalls(map, size);
	setCenterGoal(map, size, width);

	initalize_thread(&thrd, 0, 0, 0, width);
	srand (time(NULL));

	thrd.stack1 = calloc(size, sizeof(unsigned int));

	map[thrd.curr] ^= VISIT;
	curr_pos[0] = thrd.curr;

	if(print)
		displayMaze(map, size, width);

	srand (time(NULL));
	while(done == 0){

		ncells = randCell(map, size, width, pseudorandom, 0, thrd.curr, &thrd.next);

		if(ncells  > 0){
			removeWall(map, size, width, thrd.curr, thrd.next);
			thrd.prev = thrd.curr;
			thrd.curr = thrd.next;
			curr_pos[0] = thrd.curr;
			map[thrd.curr] =  map[thrd.curr] & ~VISIT;

			if(ncells > 1)
				push(thrd.stack1, &thrd.istack1, thrd.prev);
				//addRandCells(map, size, width, thrd.stack1, &thrd.istack1, thrd.stack2, thrd.prev);
		}
		else{
			if(linearsearchuint(thrd.stack1, thrd.istack1, thrd.curr, &ist) == 0){
				delete(thrd.stack1, &thrd.istack1, ist);
				if(thrd.istack1 == 0) done = 1;
			}

			if(thrd.istack1){
				rnid = rand() % thrd.istack1;
				thrd.curr = thrd.stack1[rnid];
				curr_pos[0] = thrd.curr;
			}	
		}
		
		

		if(print)
			refreshMaze(map, size, width);
	}

	free(thrd.stack1);

	if(print == 0)
		displayMaze(map, size, width);

}

void wilson(unsigned short * map, const int size, const int width, int pseudorandom, int loops, int print){
	int i, j, c, n, s, start, end, original, done = 0;
	struct thread thrd;

	setAllWalls(map, size);
	setCenterGoal(map, size, width);

	initalize_thread(&thrd, 0, 0, 0, width);
	srand (time(NULL));

	thrd.stack1 = calloc(size, sizeof(unsigned int)); //Used to determine if cell touched is in 
	thrd.stack2 = calloc(size, sizeof(unsigned int)); //Used for random-loop traversal

	push(thrd.stack1, &thrd.istack1, thrd.curr);

	map[thrd.curr] ^= VISIT;

	while(done == 0){
		for(i = 0; i < size; i++){
			s = linearsearchuint(thrd.stack1, thrd.istack1, i, &j);
			if(s == -1){
				start = i;
				original = start - 1;
				c = start;
				break;
			}
		}

		if((i == size) && (s != -1))
			done = 1;

		while(linearsearchuint(thrd.stack1, thrd.istack1, c, &i) == -1){
			map[c] = map[c] & ~VISIT;
			push(thrd.stack2, &thrd.istack2, c);
			randCell(map, size, width, pseudorandom, 1, c, &n);
			
			if(linearsearchint(thrd.stack2, thrd.istack2, n, &i) == 0){
				for(j = thrd.istack2 - 1; j > i; j--){
					pop(thrd.stack2, &thrd.istack2);
					map[thrd.stack2[j]] |= VISIT;
					
					if(print)
						refreshMaze(map, size, width);
				}

				c = thrd.stack2[thrd.istack2 - 1];
			}
			else{
				c = n;
				if(print)
					refreshMaze(map, size, width);
			}

		}

		removeWall(map, size, width, c, thrd.stack2[thrd.istack2 - 1]);
		for(j = thrd.istack2 - 1; j > 1; j--){
			removeWall(map, size, width, thrd.stack2[j], thrd.stack2[j - 1]);
			push(thrd.stack1, &thrd.istack1, thrd.stack2[j]);
			pop(thrd.stack2, &thrd.istack2);
		}

		if(loops) removeWall(map, size, width, original, start);

		if(print)
			refreshMaze(map, size, width);

	}

	free(thrd.stack1);
	free(thrd.stack2);

	if(print == 0)
		displayMaze(map, size, width);

}
