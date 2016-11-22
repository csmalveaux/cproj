/*!\headerfile floodfill.h "floodfill.h
   \brief Flood-fill algorithm header file

   Contains flood-fill functions
*/

#ifndef FLOODFILL_H
#define FLOODFILL_H

#include "data.h"

extern unsigned char ff[256];

void init_flood (void);
unsigned char min_cell(int);
void update_flood (void);
void ff_goto (void);

#endif // FLOODFILL_H