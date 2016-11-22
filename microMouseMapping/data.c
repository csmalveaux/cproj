/*!\file data.cpp
 * \brief Data C++ File
 *
 * Defining global variables.
*/

#include "data.h"

/*! \var extern unsigned char mmap[256]
	\brief  Array which contains each maze cell data.
 *	  
 *	Array which represents all the cells in teh 16 by 16 cell micromoues maze. Each cell is updated as a bit field.
 *		Example: 
			mmap[x] |= NORTH + WEST + VISIT;
			inicates that at cell x, there is a north and west wall and it has been visited.
 *	In order to remove a wall or to mark the cell as unvisited:
			mmap[x]&= ~(WEST + VISIT);
*/
unsigned char mmap[256];

/*! \var extern int x
	\brief  Current x-coordinate of mouse.*/
int x = 0;
/*! \var extern int y
  * \brief  Current x-coordinate of mouse.*/
int y = 0;
/*! \var extern int gx
	\brief  Current x-coordinate of mouse.
* Default gy = 15;*/
int gx = 7;
/*! \var extern int gy
  * \brief  Current x-coordinate of mouse.
  * Default gy = 15;*/
int gy = 7;
/*! \var extern int alg
  * \brief  Current algorithm used.
  * Default flood fill.;*/
int alg = FLOOD;