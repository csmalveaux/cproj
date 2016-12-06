#ifndef MAP_H
#define MAP_H

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

#endif //MAP_H