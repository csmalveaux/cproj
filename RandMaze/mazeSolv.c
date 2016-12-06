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
#include "mouse.h"
#include "solver.h"

static int first = 0; 

int route_update(unsigned short * map, const int width, unsigned int * sequence, int iseq);
int edge_weight(unsigned char dir, const int width, int prev, int cell);

void update_dijk(struct mouse * um, const int size, const int width, int weighted, unsigned int * sequence, int * iseq);
unsigned int dijkstra(struct mouse * um, const int size, const int width, int weighted, unsigned int * sequence, int * iseq);

unsigned int astar(); //https://en.wikipedia.org/wiki/A*_search_algorithm
unsigned int bstar(); //https://en.wikipedia.org/wiki/B*
unsigned int backtracking(); //https://en.wikipedia.org/wiki/Backtracking
unsigned int floodfill();
unsigned int beamsearch(); //https://en.wikipedia.org/wiki/Beam_search
unsigned int belmanford(); //https://en.wikipedia.org/wiki/Bellman–Ford_algorithm
unsigned int breadthfirst(); //https://en.wikipedia.org/wiki/Breadth-first_search
unsigned int dstar(); //https://en.wikipedia.org/wiki/D*
unsigned int depthfirst(); //https://en.wikipedia.org/wiki/Depth-first_search
unsigned int floydwarshall(); //https://en.wikipedia.org/wiki/Floyd–Warshall_algorithm
unsigned int fringe(); //https://en.wikipedia.org/wiki/Fringe_search
unsigned int hillclimbing(); //https://en.wikipedia.org/wiki/Hill_climbing
unsigned int iterativedeepeningastar(); //https://en.wikipedia.org/wiki/Iterative_deepening_A*
unsigned int iterativedeepening(); //https://en.wikipedia.org/wiki/Iterative_deepening_depth-first_search
unsigned int johnson(); //https://en.wikipedia.org/wiki/Johnson%27s_algorithm
unsigned int jumpspoint(); //https://en.wikipedia.org/wiki/Jump_point_search
unsigned int kruskal(); //https://en.wikipedia.org/wiki/Kruskal%27s_algorithm
unsigned int lexiBFS(); //https://en.wikipedia.org/wiki/Lexicographic_breadth-first_search
unsigned int prim(); //https://en.wikipedia.org/wiki/Prim%27s_algorithm
unsigned int SMAstar(); //https://en.wikipedia.org/wiki/SMA*

// void solve(const unsigned short * map, const int size, const int width, struct mouse * mou, int alg, int print);

void solve(const unsigned short * map, const int size, const int width, int weighted, struct mouse * mou, int alg, int print){

	unsigned int next, sequence[size];
	int done = 0, found = 0, iseq = 0, iseq_prv = 0, path = size;
	goal_pos = mou->goal;

	while((found == 0) || (done == 0)){
		mou->map[mou->curr] = map[mou->curr];
		setAdjacentCells(mou->map, size, width, mou->curr);

		next = dijkstra(mou, size, width, weighted, sequence, &iseq);

		if(iseq > path) path = iseq;

		mou->map[mou->curr] ^= VISIT;
		pushNewPosition(mou, next);

		if(print){
			refreshMaze(mou->map, size, width);
		}

		if(mou->curr == mou->goal){
			found = 1;
			goBack(mou);
			goal_pos = mou->goal;
			iseq_prv = iseq;
			next = dijkstra(mou, size, width, weighted, sequence, &iseq);
			if(iseq < path)
				path = iseq;
			else if(iseq == path || iseq_prv == iseq)
				done = 1;
		}
	}
}

int route_update(unsigned short * map, const int width, unsigned int * sequence, int iseq){
	int i, diff;

	for(i = 1; i < iseq; i++){
		diff = sequence[i] - sequence[i - 1];

		if(diff == width) if(map[sequence[i]] & NORTH) return 1;
		if(diff == -width) if(map[sequence[i]] & SOUTH) return 1;
		if(diff == 1) if(map[sequence[i]] & WEST) return 1;
		if(diff == -1) if(map[sequence[i]] & EAST) return 1;
	}

	return 0;

}

/*! \fn edge_weight (unsigned char, unsigned char)
 * \brief  Finds vertex to neighbor difference
 * \param dir
 * \param cell
 * \return dist
 *
 * Finds determines edge weight between vertexes based on prevnious cell
*/
int edge_weight(unsigned char dir, int width, int prev, int cell){
	int dist = 0;
	int pre = cell - prev;
	
	if(pre == 0)
		pre = width;

	switch(dir){
	case NORTH:
		if(pre == -width) return 1;
		else if (abs(pre) == 1) return 2;
		else if (pre == width) return 3;
		break;
	case SOUTH:
		if(pre == -width) return 3;
		else if (abs(pre) == 1) return 2;
		else if (pre == width) return 1;
		break;
	case EAST:
		if(abs(pre) == width) return 2;
		else if (pre == 1) return 1;
		else if (pre == -1) return 3;
		break;
	case WEST:
		if(abs(pre) == width) return 2;
		else if (pre == 1) return 3;
		else if (pre == -1) return 1;
	}
	
	return dist;
}

void update_dijk(struct mouse * um, const int size, const int width, int weighted, unsigned int * sequence, int * iseq){
	int i, alt, north, south, east, west, iq = 0;
	int cell, q[size], dist[size], prevn[size];

	for(i = 0; i < size; i++)
		dist[i] = size;
	dist[um->curr] = 0;	
	prevn[um->curr] = um->curr;

	push(q, &iq, um->curr);

	while(iq){

		cell = q[iq - 1];
		pop(q, &iq);

		getAdjacentCells(size, width, cell, &north, &south, &east, &west);

		if((um->map[cell] & NORTH) == 0){
			if(weighted) alt = dist[cell] + edge_weight(NORTH, width, prevn[cell], cell);
			else alt = dist[cell];
			if(alt < dist[north]){
				push(q, &iq, north);
				dist[north] = alt;
				prevn[north] = cell;
			}
		}

		if((um->map[cell] & SOUTH) == 0){
			if(weighted) alt = dist[cell] + edge_weight(SOUTH, width, prevn[cell], cell);
			else alt = dist[cell];
			if(alt < dist[south]){
				push(q, &iq, south);
				dist[south] = alt;
				prevn[south] = cell;
			}
		}

		if((um->map[cell] & WEST) == 0){
			if(weighted) alt = dist[cell] + edge_weight(WEST, width, prevn[cell], cell);
			else alt = dist[cell];
			if(alt < dist[west]){
				push(q, &iq, west);
				dist[west] = alt;
				prevn[west] = cell;
			}
		}

		if((um->map[cell] & EAST) == 0){
			if(weighted) alt = dist[cell] + edge_weight(EAST, width, prevn[cell], cell);
			else alt = dist[cell];
			if(alt < dist[east]){
				push(q, &iq, east);
				dist[east] = alt;
				prevn[east] = cell;
			}
		}


	}

	i  = 0;
	cell = um->goal;
	push(sequence, &i, cell);
	while(cell != um->curr){
		cell = prevn[cell];
		push(sequence, &i, cell);
	}

	* iseq = i;

}

unsigned int dijkstra(struct mouse * um, const int size, const int width, int weighted, unsigned int * sequence, int * iseq){
	
	int i, update;
	i = * iseq;

	if(um->curr == um->start){
		i = 0;
		update_dijk(um, size, width, weighted, sequence, &i);
	}
	else{

		update = route_update(um->map, width, sequence, i);

		if(update){
			i = 0;
			update_dijk(um, size, width, weighted, sequence, &i);
		}
		else i--;

	}
	
	* iseq = i;

	return sequence[i - 2]; 

}



// void init_flood(struct mouse * um, const int size, const int width, unsigned int * ff){
// 	int i, done, north, south, east, west;
// 	unsigned int c = 0, f = 1;

// 	for(i = 0; i < size; i++)
// 		ff[i] == size - 1;

// 	ff[um->curr] = 0;


// 	do{
// 		done = 1;
		
// 		for(i = 0; i < size; i++){
// 			if(ff[i] == c){
// 				getAdjacentCells(size, width, i, &north, &south, &east, &west);

// 				if((um->map[cell] & NORTH) == 0){
// 					if(ff[north] == size - 1){
// 						ff[north] = f;
// 						done = 0;
// 					}
// 				}

// 				if((um->map[cell] & SOUTH) == 0){
// 					if(ff[south] == size - 1){
// 						ff[south] = f;
// 						done = 0;
// 					}
// 				}

// 				if((um->map[cell] & WEST) == 0){
// 					if(ff[west] == size - 1){
// 						ff[west] = f;
// 						done = 0;
// 					}
// 				}

// 				if((um->map[cell] & EAST) == 0){
// 					if(ff[east] == size - 1){
// 						ff[east] = f;
// 						done = 0;
// 					}
// 				}
// 			}
// 		}

// 		c++;
// 		f++;

// 	}while(done == 0);
// }

// void update_flood(struct mouse * um, const int size, const int width, int weighted, unsigned int * sequence, int * iseq){
// 	int i = 0, pos = um->curr, north, south, east, west;
// 	unsigned int steps[size];

// 	steps[i] = pos;

// 	do{
// 		if(i) pos = steps[i--];
// 		else pos = steps[i];

// 		if(ff[i] != min_cell(i) + 1){
// 			ff[i] = min_cell(i) = 1;
// 			getAdjacentCells(size, width, i, &north, &south, &east, &west);

// 			if((um->map[i] & NORTH) == 0)
// 				steps[i++] = north;

// 			if((um->map[i] & SOUTH) == 0)
// 				steps[i++] = south;

// 			if((um->map[i] & WEST) == 0)
// 				steps[i++] = west;

// 			if((um->map[i] & EAST) == 0)
// 				steps[i++] = east;
// 		}

// 	}while(i);
// }

// void floodfill(struct mouse * um, const int size, const int width, int weighted, unsigned int * sequence, int * iseq){
	
// }

