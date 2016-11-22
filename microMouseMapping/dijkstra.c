/*!\file dijkstra.cpp
 * \brief Dijkstra algorithm C++ File
 *
 * Contains dijkstra algorithm functions
*/

#include "dijkstra.h"

/*! \var unsigned char dijk[256]
	\brief  Array which contains each cell dijkstra distance from goal.*/
unsigned char dijk[256];

/*! \var unsigned char prevn[256]
	\brief  Array which contains prevnious node for each node.*/
unsigned char prevn[256];

/*! \var unsigned char seq[256]
	\brief  Array which contains sequence from source to target.*/
unsigned char seq[256];

/*! \var unsigned char seq_it
	\brief  Sequence arrary iterator*/
unsigned char seq_it;

/* \var int r_update
   \brief  Indicates if route has changed.*/
int r_update = 0;

/*!\fn init_dijk (void)
 * \brief  Intializes the flood fill values
 * \return void
 *
 * Initializes dijk[] array to 256.
*/
void init_dijk (void){
	
	for(int i = 0; i < 16; i++)
		for(int k = 0; k < 16; k++)
			dijk[k + i*16] = 0xFF;

	prevn[0] = 0;
	dijk[x + y*16] = 0;

}

/*!\fn dijk_goto (void)
 * \brief  Dijkstra next move function
 * \return void
*/
void dijk_goto (void){
	
	int pos = x + y*16;
	int u;
	static int it = 0;

	mmap[pos] |= VISIT;
	
	if(r_update)
		it = 1;
	else
		it++;

	u = seq[seq_it - it]; 

	prevn[u] = pos;

	x = u % 16;
	y = u / 16;

}

/*! \fn update_dijk
 * \brief  Updates dijkstra values of adjacent cells.
 * \return void
*/
void update_dijk (void){
	
	unsigned char u = x + y*16;
	unsigned char temp;
	int alt, q_it, found;
	unsigned char q[256];		// All nodes in the graph are unoptimized - thus are in Q
	unsigned char dijk[256];
	
	for(int i = 0; i < 256; i++)
		dijk[i] = 0xFF;

	dijk[x + y*16] = 0;

	found = 0;
	q_it = 0;
	
	q[q_it] = u;

	do{
		
		if(q_it)
			q_it = q_it - 1;
		
		u = q[q_it];

		//u = min_vertex(u);		// Vertex in Q with the smallest dist[]
			
		//if(dijk[u] == 0xFF)		// All remaining vertices are inaccessible from source
		//	break;
						
		// If u == target
		//if(u == (gx + gy*16)){
		//	found = 1;
		//	break;
		//}

		temp = mmap[u];

		// For each neighbor u
		if(!((temp&NORTH)==NORTH)){
			
			alt = dijk[u] + edge_weight(NORTH, u);
			if(alt < dijk[u-16]){
				q[q_it++] = u - 16;
				dijk[u-16] = alt;
				prevn[u-16] = u;
			}
		}
		if(!((temp&EAST)==EAST)){
			
			alt = dijk[u] + edge_weight(EAST, u);
			if(alt < dijk[u+1]){
				q[q_it++] = u + 1;
				dijk[u+1] = alt;
				prevn[u+1] = u;
			}
		}
		if(!((temp&SOUTH)==SOUTH)){
			
			alt = dijk[u] + edge_weight(SOUTH, u);
			if(alt < dijk[u+16]){
				q[q_it++] = u + 16;
				dijk[u+16] = alt;
				prevn[u+16] = u;
			}
		}
		if(!((temp&WEST)==WEST)){
			
			alt = dijk[u] + edge_weight(WEST, u);
			if(alt < dijk[u-1]){
				q[q_it++] = u - 1;
				dijk[u-1] = alt;
				prevn[u-1] = u;
			}
		}
		
	}while(q_it);		

	//if(found){
		seq_it = 0;		// Empty sequence
		u = gx + gy*16;	// Target
		while (u != (x + y*16)){
			seq[seq_it] = u;
			u = prevn[u];
			seq_it++;
		}
	//}

}


/*! \fn min_vertex (int pos)
 * \brief  Finds lowest valued vertex.
 * \param pos
 * \return int next_node
 *
 * Finds lowest value edge and returns the vertex.
*/
int min_vertex(int pos){
	
	unsigned char temp = mmap[pos];
	int min_vert = 256;
	int vertex;
	int next_node;

	if(!((temp&NORTH)==NORTH)){
		vertex = dijk[pos] + edge_weight(NORTH, pos);
		if(min_vert > vertex){
			min_vert = vertex;
			next_node = pos - 16;
		}
	}
	if(!((temp&SOUTH)==SOUTH)){
		vertex = dijk[pos] + edge_weight(SOUTH, pos);
		if(min_vert > dijk[pos + 16]){
			min_vert = vertex;
			next_node = pos + 16;
		}
	}
	if(!((temp&EAST)==EAST)){
		vertex = dijk[pos] + edge_weight(EAST, pos);
		if(min_vert > dijk[pos + 1]){
			min_vert = vertex;
			next_node = pos + 1;
		}
	}
	if(!((temp&WEST)==WEST)){
		vertex = dijk[pos] + edge_weight(WEST, pos);
		if(min_vert > dijk[pos - 1]){
			min_vert = vertex;
			next_node = pos - 1;
		}
	}

	dijk[next_node] = min_vert;
	prevn[next_node] = pos;

	return next_node;
}

/*! \fn edge_weight (unsigned char, unsigned char)
 * \brief  Finds vertex to neighbor difference
 * \param dir
 * \param cell
 * \return dist
 *
 * Finds determines edge weight between vertexes based on prevnious cell
*/
int edge_weight(unsigned char dir, unsigned char cell){
	int dist = 0;
	int pre = cell - prevn[cell];
	
	if(pre == 0)
		pre = 16;

	switch(dir){
	case NORTH:
		if(pre == -16)
			dist = 1;
		else if (pre == 1 || pre == -1)
			dist = 2;
		else if (pre == 16)
			dist = 3;
		break;
	case SOUTH:
		if(pre == -16)
			dist = 3;
		else if (pre == 1 || pre == -1)
			dist = 2;
		else if (pre == 16)
			dist = 1;
		break;
	case EAST:
		if(pre == -16||pre == 16)
			dist = 2;
		else if (pre == 1)
			dist = 1;
		else if (pre == -1)
			dist = 3;
		break;
	case WEST:
		if(pre == -16||pre == 16)
			dist = 2;
		else if (pre == 1)
			dist = 3;
		else if (pre == -1)
			dist = 1;
	}
	return dist;
}

/*! \fn route_update (void)
 * \brief  Determins if the maze needs to be updated.
 * \return int
 *
 * Taking into account the current route to target cell, 
 * if there are new walls which compromise the set path, 
 * the maze is recalculated.
*/
int route_update (void){
	unsigned char temp;
	int u, d;
		
	for(int i = 0; i < seq_it; i++){
		u = seq[i];
		d = seq[i] - seq[i-1];
		temp = mmap[u];
		if(d == 16)
			if(((temp&NORTH) == NORTH))
				return 1;
		if(d == -16)
			if(((temp&SOUTH) == SOUTH))
				return 1;
		if(d == 1)
			if(((temp&WEST) == WEST))
				return 1;
		if(d == -1)
			if(((temp&EAST) == EAST))
				return 1;
	}
	return 0;
}