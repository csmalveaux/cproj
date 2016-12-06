#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "maze.h"
#include "map.h"
#include "mouse.h"
#include "solver.h"

int main(int argc, char *argv[] )  {
	int i, size, height, width, animate, type;

	if (argc == 2) {
       size = atoi(argv[1]);
       height = size;
       width =  size;
       type = 0;
    }
    else if (argc == 3) {
       height = atoi(argv[1]);
       width = atoi(argv[2]);
       animate = 0;
       type = 0;
    }
    else if (argc == 4) {
       height = atoi(argv[1]);
       width = atoi(argv[2]);
       type = atoi(argv[3]);
       animate = 0;
    }
    else if (argc == 5) {
       height = atoi(argv[1]);
       width = atoi(argv[2]);
       type = atoi(argv[3]);
       animate = atoi(argv[4]);
       sleep_time /= (height*width);

    }
    else if (argc >= 6) { 
       height = atoi(argv[1]);
       width = atoi(argv[2]);
       type = atoi(argv[3]);
       animate = atoi(argv[4]);
       sleep_time = atoi(argv[5]);
    }
    
    size = height*width;

	unsigned short map[size];
  getCorners(size, width);

	switch (type){
		case 0: 
			recursiveBacktracker(map, size, width, 1, 0, 0, animate);
			break;
		case 1: 
			recursiveBacktracker(map, size, width, 6, 0, 0, animate);
			break;
		case 2:
			randomizedTraversal(map, size, width, 0, 1, animate);
      break;
    case 3:
      randomizedTraversal(map, size, width, 0, 6, animate);
      break;
    case 4:
      randomizedPrim(map, size, width, 0, 0, animate);
      break;
    case 5:
      wilson(map, size, width, 0, 0, animate);

	}

  clearCursors();
  
  struct mouse umouse;

  intialize_mouse(&umouse, 0, (height/2 - 1) * width + width/2 - 1);
  umouse.map = calloc(size, sizeof(unsigned short));
  setPerimeter(umouse.map, size, width);

  refreshMaze(umouse.map, size, width);
  sleep_time *= 10;
  solve(map, size, width, 1, &umouse, 1, animate);

  free(umouse.map);

}