/*!\file floodfill.cpp
 * \brief Flood-fill algorithm C++ File
 *
 * Contains flood-fill algorithm functions.
 */

#include "floodfill.h"

/*! \var extern unsigned char ff[256]
	\brief  Array which contains each cell flood-fill distance from goal.*/
unsigned char ff[256];

/*! \fn init_flood (void)
 * \brief  Intializes the flood fill values
 * \return void
 *
 * Initializes ff[] array to 255.
*/
void init_flood (void){
	
	int m = 0;
	int f = 1;
	int done = 0;
	int tc;
	unsigned char temp;

	for(int i = 0; i < 16; i++)
		for(int k = 0; k < 16; k++)
			ff[k + i*16] = 255;
	
	ff[gx + gy*16] = 0;
	do{
		done = 1;
		for(int i = 0; i < 16; i++)
			for(int k = 0; k < 16; k++)
				if(ff[k + i*16] == m){
					tc = k + i*16;
					temp = mmap[tc];
					if(!((temp&NORTH)==NORTH))
						if(ff[tc - 16] == 255){
							ff[tc - 16] = f;
							done = 0;
						}
					if(!((temp&WEST)==WEST))
						if(ff[tc - 1] == 255){
							ff[tc - 1] = f;
							done = 0;
						}
					if(!((temp&SOUTH)==SOUTH))
						if(ff[tc + 16] == 255){
							ff[tc + 16] = f;
							done = 0;
						}
					if(!((temp&EAST)==EAST))
						if(ff[tc + 1] == 255){
							ff[tc + 1] = f;
							done = 0;
						}
				}
		m++;
		f++;
	}while(!done);
}

/*! \fn min_cell (int pos)
 * \brief  Finds lowest flood-filled valued adjacent cell.
 * \param pos
 * \return unsigned char min_cell
 *
 * Finds lowest flood-filled valued adjacent cell.
*/
unsigned char min_cell(int pos){
	 
	unsigned char min_cell = 255;
	unsigned char temp = mmap[pos];

	if(!((temp&NORTH)==NORTH))
		if(ff[pos - 16] < min_cell)
			min_cell = ff[pos - 16];
	if(!((temp&SOUTH)==SOUTH))
		if(ff[pos + 16] < min_cell)
			min_cell = ff[pos + 16];
	if(!((temp&WEST)==WEST))
		if(ff[pos - 1] < min_cell)
			min_cell = ff[pos - 1];
	if(!((temp&EAST)==EAST))
		if(ff[pos + 1] < min_cell)
			min_cell = ff[pos + 1];
	
	return min_cell;
}

/*! \fn update_flood (void)
 * \brief  Updates flood-filled values of adjacent cells.
 * \return void
*/
void update_flood (void){
	
	int p = 0;
	int pos = x + y*16;
	unsigned char st[255];
	unsigned char temp;

	st[p] = pos;
	
	do{
		if(p)
			pos = st[p--];
		else
			pos = st[p];
		temp = mmap[pos];
		if(ff[pos] != (min_cell(pos) + 1)){
			ff[pos] = min_cell(pos) + 1;
			if(!((temp&NORTH)==NORTH))
				st[p++] = pos - 16;
			if(!((temp&SOUTH)==SOUTH))
				st[p++] = pos + 16;
			if(!((temp&WEST)==WEST))
				st[p++] = pos - 1;
			if(!((temp&EAST)==EAST))
				st[p++] = pos + 1;
		}
	}while(p);
}

/*! \fn ff_goto (void)
 * \brief  Flood-fill next move function
 * \return void
*/
void ff_goto (void){
	
	unsigned char temp = mmap[x + y*16];
	int n, e, w, s;
	
	s = temp&SOUTH;
	n = temp&NORTH;
	e = temp&EAST;
	w = temp&WEST;

	mmap[x+y*16] |= VISIT;

	int pos = x + y*16;
	unsigned char mc = min_cell(pos);
	if(!(n == NORTH))
		if(mc==ff[pos-16]){
			y--;
			return;
		}
	if(!(w == WEST))
		if(mc==ff[pos-1]){
			x--;
			return;
		}
	if(!(s == SOUTH))
		if(mc==ff[pos+16]){
			y++;
			return;
		}
	if(!(e == EAST))
		if(mc==ff[pos+1]){
			x++;
		}
}
