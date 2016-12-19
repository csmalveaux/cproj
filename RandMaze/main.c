#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "maze.h"
#include "map.h"
#include "mouse.h"
#include "solver.h"
#include "readwrite.h"

int main(int argc, char *argv[] )  {
	int i, size, height, width, animate, type, solve, save;

	if (argc == 2) {
       size = atoi(argv[1]);
       height = size;
       width =  size;
       type = 0;
       animate = 0;
       solve = 0;
       save = 0;
    }
  else if (argc == 3) {
    height = atoi(argv[1]);
    width = atoi(argv[2]);
    animate = 0;
    type = 0;
    solve = 0;
    save = 0;
  }
  else if (argc == 4) {
    height = atoi(argv[1]);
    width = atoi(argv[2]);
    type = atoi(argv[3]);
    animate = 0;
    solve = 0;
    save = 0;
  }
  else if (argc == 5) {
    height = atoi(argv[1]);
    width = atoi(argv[2]);
    type = atoi(argv[3]);
    animate = atoi(argv[4]);
    sleep_time /= (height*width);
    solve = 0;
    save = 0;
  }
  else if (argc >= 6) { 
    height = atoi(argv[1]);
    width = atoi(argv[2]);
    type = atoi(argv[3]);
    animate = atoi(argv[4]);
    sleep_time = atoi(argv[5]);
    solve = 0;
    save = 0;
  }
  else if (argc >= 7) { 
    height = atoi(argv[1]);
    width = atoi(argv[2]);
    type = atoi(argv[3]);
    animate = atoi(argv[4]);
    sleep_time = atoi(argv[5]);
    save = atoi(argv[6]);
    solve = 0;
  }
  else if (argc >= 8) { 
    height = atoi(argv[1]);
    width = atoi(argv[2]);
    type = atoi(argv[3]);
    animate = atoi(argv[4]);
    sleep_time = atoi(argv[5]);
    save = atoi(argv[6]);
    solve = atoi(argv[7]);
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

  if(solve){
    struct mouse umouse;
    intialize_mouse(&umouse, 0, (height/2 - 1) * width + width/2 - 1);
    umouse.map = calloc(size, sizeof(unsigned short));
    setPerimeter(umouse.map, size, width);

    refreshMaze(umouse.map, size, width);
    sleep_time *= 10;
    solve(map, size, width, 1, &umouse, 1, 1);
  }

  if(save){
    writeMaze(map, size, width);
  }

  int size2, width2;

  readMazeSize(&size2, &width2, "1161118_25147.mz");

  unsigned short map2[size2];

  readMaze(map2, size2, width2, "1161118_25147.mz");
  displayMaze(map2, size2, width2);

}