/*!\headerfile floodfill.h "floodfill.h
   \brief Flood-fill algorithm header file

   Contains flood-fill functions
*/

#ifndef MAPPING_H
#define MAPPING_H

#include "data.h"

using namespace std;

void init (void);
int wsensors (void);
void mwupdate (void);
void mgoto(void);
unsigned char sideassign (string);
void update (void);
int roamg (void);
void init_alg (void);
void disp_route (void);

#endif // MAPPING_H