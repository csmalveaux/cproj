/*!\headerfile data.h "data.h
   \brief Global data and definitions header file

   Contains global data and definitions
*/

#ifndef DATA_H
#define DATA_H

#include "menu.h"
//#include "mapping.h"
//#include "floodfill.h"

using namespace std;


/************
!!CELL DATA!!
*************/
/*!\def NORTH 
 * \brief Bit represents  north wall in cell.*/
#define NORTH  0x01 	
/*!\def SOUTH
 * \brief Bit represents  south wall in cell.*/
#define SOUTH  0x04
/*!\def EAST 
 * \brief Bit represents  east wall in cell.*/
#define EAST   0x02		
/*!\def WEST 
 * \brief Bit represents  west wall in cell.*/
#define WEST   0x08
/*!\def VISIT 
 * \brief Bit represents  if the cell was visited.*/	
#define VISIT  0x10

#define TOPLEFT 119
#define TOPRGHT 120
#define BTMLEFT 135
#define BTMRGHT 136

/*************
!!MENU MODES!!
**************/
/*!\def EXIT 
 * \brief Exit mode.*/
#define EXIT	0
/*!\def INIT 
 * \brief Initialize mode.*/
#define INIT	1
/*!\def DISP 
 * \brief Display maze mode.*/
#define DISP	2
/*!\def EXPLR 
 * \brief Explore maze mode.*/
#define EXPLR	3
/*!\def SOLVE 
 * \brief Solve maze mode.*/
#define SOLVE	4
/*!\def CLEAR 
 * \brief Clear screen mode.*/
#define CLEAR	5
/*!\def SETT 
 * \brief Settings mode.*/
#define SETT	6
/*!\def HELP 
 * \brief Help mode.*/
#define HELP	7
/*!\def SAVE 
 * \brief Save file.*/
#define SAVE	8
/*!\def LOAD 
 * \brief Load file.*/
#define LOAD	9
/*!\def STAT 
 * \brief Status mode.*/
#define STAT	10
/*!\def STATIS 
 * \brief Statistic mode.*/
#define STATIS	11
/*!\def NOINST 
 * \brief Not an instruction.*/
#define NOINST	12

/*************
!!ALGORITHMS!!
**************/
/*!\def FLOOD 
 * \brief Flood fill algorithm.*/
#define FLOOD	1
/*!\def DIJK 
 * \brief Dijkstra algorithm.*/
#define DIJK	2

extern unsigned char mmap[256];

extern int x;
extern int y;
extern int gx;
extern int gy;
extern int alg;

#endif // DATA_H