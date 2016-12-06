#include <stdio.h>
#include "map.h"
#include "maze.h"
#include "mouse.h"


void displayMaze(unsigned short * map, int sz, int w);
void displayRow(unsigned short * map, int sz, int w, int row);

void displayBoard();
void displayNoBoard();
void displayCurr();
void displayVis();
void displayEmpty();
void displayGoal();

void clearMaze(int sz, int w);
void refreshMaze(unsigned short * map, const int sz, const int w);

void clearCursors();
void displayMouse(const int width);

void displayMaze(unsigned short * map, const int sz, const int w){
	int i;

	for(i = 0; i < sz/w; i++)
		displayRow(map, sz, w, i);

		printf("\n");

}

void displayRow(unsigned short * map, int sz, int w, int row){
	int i, j, k, pos = 0;
	for(j = 0; j < 3; j++){
		for(i = 0; i < w; i++){
			switch(j){
				case 0:
					if(row == 0){ 
						if(i == 0) printf("┌"); 
							  
						if ((map[row*w + i] & NORTH) == NORTH) displayBoard();
						else displayNoBoard();

						if(i == w - 1) printf("┐");
						else  printf("┬");  

					}
					else{
						if(i == 0) printf("├");  

						if ((map[row*w + i] & NORTH) == NORTH) displayBoard();
						else displayNoBoard();

						if(i == w - 1) printf("┤");
						else  printf("┼");  
					}
					break;
				
				case 1: 
					if(i == 0) { 
						if ((map[row*w + i] & WEST) == WEST) printf("│"); 
						else printf(" ");
					}
					
					for(k = 1; k <= MAX_SEEDS; k++){
						if(curr_pos[k -  1] == (row * w + i)){
							pos = k;
							break;
						}
					}

					if(pos > 0){ displayCurr(pos); pos = 0; }
					else if(mouse_pos == (row * w + i)) displayMouse(w);
					else if(goal_pos == (row * w + i)) displayGoal(w);
					else if((map[row*w + i] & VISIT) == VISIT) displayVis();
					else displayEmpty();

					if((map[row*w + i] & EAST) == EAST) printf("│"); 
					else printf(" ");
					break;

				case 2:
					if(row == sz/w - 1){ 
						if(i == 0) printf("└"); 
							  
						if ((map[row*w + i] & SOUTH) == SOUTH) displayBoard();
						else displayNoBoard();

						if(i == w - 1) printf("┘");
						else  printf("┴");  

					}
			}
		}
		if((j < 2)) printf("\n");
	}

}

void displayBoard(){
	printf("───");
}

void displayNoBoard(){
	printf("   ");
}

void displayMouse(int width){
	if      (mouse_dir == -width) printf(" ▼ ");
	else if (mouse_dir == width) printf(" ▲ ");
	else if (mouse_dir == 1) printf(" ◀︎ ");
	else if (mouse_dir == -1) printf(" ▶︎ ");
	else printf(" ◼︎ ");
}

void displayCurr(int pos){
	switch(pos){
		case 1:
			printf(" ■ ");
			break;
		case 2:
			printf(" ● ");
			break;
		case 3:
			printf(" ✖︎ ");
			break;
		case 4:
			printf(" ❖ ");
			break;
		case 5:
			printf(" ❤︎ ");
			break;
		case 6:
			printf(" ◎ ");
	}
	
}

void displayGoal(){
	printf(" ★ ");
}

void displayVis(){
	printf(" ▫ ");
}

void displayEmpty(){
	printf("   ");
}

void clearMaze(int sz, int w){
	int i;
	for(i = 0; i <= 2*sz/w; i++)
		printf("\033[A\r");
}

void refreshMaze(unsigned short * map, const int sz, const int w){
	clearMaze(sz, w);
	displayMaze(map, sz, w);
	usleep(sleep_time);
}

void clearCursors(){
	int i;
	for(i = 0; i < MAX_SEEDS; i++)
		curr_pos[i] = -1;
}
