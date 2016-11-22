/*!\headerfile dijkstra.h "dijkstra.h
   \brief Dijkstra algorithm header file

   Contains dijkstra functions
*/

#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include "data.h"

extern unsigned char dijk[256];
extern unsigned char prevn[256];
extern unsigned char seq[256];
extern unsigned char q[256];
extern unsigned char seq_it;

extern int r_update;

void init_dijk (void);
void update_dijk (void);
void dijk_goto (void);
int min_vertex (int);
int edge_weight (unsigned char, unsigned char);
int route_update (void);

#endif // DIJKSTRA_H