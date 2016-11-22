/*!\headerfile menu.h "menu.h
   \brief Menu Header File
*/

#ifndef MENU_H
#define MENU_H

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <windows.h>
#include "data.h"

using namespace std;

extern int mode;
extern int first;

int menu (void);
int commandsel (string);
void commdisp (void);
void mazedisp (void);
void helpdisp (string);

#endif // MENU_H