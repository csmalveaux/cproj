/*!\file mapping.cpp
 * \brief Mapping C++ File
 *
 * Contains mapping functions
*/

#include "mapping.h"
#include "data.h"
#include "floodfill.h"
#include "dijkstra.h"

/*! \fn init (void)
 * \brief  Function to initialize data in mmap array
 * \return void
 *
 * Function initializes each cell in mmap array to unvisited, and establishes walls on the north, south, west, and east side of  the maze.
*/
void init (void){

	// clearing mmap with 0x00
	for(int i = 0; i < 256; i++)
		mmap[i] = 0x00;
	
	// setting up boundary walls
	for(int i = 0; i < 16; i++){
		mmap[i]|= NORTH;
		mmap[i*16]|=WEST;
		mmap[240+i]|=SOUTH;
		mmap[15+i*16]|=EAST;
	}
	x = 0;
	y = 0;

	init_alg();
}

/*! \fn init_alg (void)
 * \brief  Function to initialize algorithm data
 * \return void
 *
 * Function initializes data associated with specified algorithm.
*/
void init_alg (void){
	//if (alg == FLOOD)
	//	init_flood();
	//else if (alg == DIJK)
		init_dijk();
}

/*! \fn wsensors (void)
 * \brief  Takes user input to simulate mouse sensor input
 * \return int
 *
 * Function allows the user to input the walls of the current cell, in the form of a string, in lieu of the physical wall sensing circuit. The string parses out each word, separated by a whitespace. 
Five acceptable inputs: north, south, east, west, and back. Inputs not need to be in any particular order. Separate only by whitespace. Lowercase only. 
*/
int wsensors (void){
	string walls, str;
	size_t it;
	unsigned char side;
	cout << "Current route: ";
	disp_route ();
	cout.width(10); cout << left << "enter walls of current cell (north, south, east, west) \nenter 'back' to exit explore mode\ncurrent position: ("<< x <<" ,"<< y <<") \n";
	getline (cin, walls);
	while(!walls.empty()){	
		it = walls.find_first_of(" ");
		str.assign(walls, 0, it);
			
		side = sideassign(str);
		if (side != 0)
			mmap[x+y*16] |= side;
		else if (side == 0){
			mode = DISP;
			return 0;
			break;
		}
		walls.erase (0, it + 1);
		if(!str.compare(walls))
			walls.clear();
		}
		str.clear();
		return 1;
}

/*! \fn mwupdate (void)
 * \brief  Updates maze walls
 * \return void
 *
 * Function takes adjacent walls of cells and updates their walls without the cells having been visited.
 Example: mmap[2] has a east wall, thus mmap[3] will be updated to indicate a west wall, while remaing unvisited.
*/
void mwupdate (void){
	unsigned char temp;
	for(int i = 0; i < 16; i++)
		for(int k = 0; k <16; k++){
			temp = mmap[i + k*16];
			if(i != 0)
				if((temp&NORTH) == NORTH)
					mmap[i + (k-1)*16]|=SOUTH;
			if(i != 15)
				if((temp&SOUTH) == SOUTH)
					mmap[i + (k+1)*16]|=NORTH;
			if(k != 0)
				if((temp&WEST) == WEST)
					mmap[i-1 + k*16]|=EAST;
			if(k != 15)
				if((temp&EAST) == EAST)
					mmap[i+1 + k*16]|=WEST;
		}
}

/*! \fn mgoto (void)
 * \brief  Simple movement algorithm
 * \return void
 *
 * Function to decide the mouse's next movement based on the surrounding walls and visited cells around the current cell. Mouse moves based on preferences 1) south, 2) east, 3) west, 4) north. Back tracing over visisted squares still needs to be implemented.
 * Change this function to play with the algorithm. 
*/
void mgoto(void){
	
	//if (alg == FLOOD)
	//	ff_goto ();
	//else if (alg == DIJK)
		dijk_goto ();

}

/*! \fn update (void)
 * \brief  Function to update algorithm data
 * \return void
 *
 * Function updates data associated with specified algorithm.
*/
void update (void){
	//if (alg == FLOOD)
	//	update_flood();
	//else if (alg == DIJK)
		if(x == 0 && y == 0){
			update_dijk();
		}
		else{
			r_update = route_update();
			if(r_update || roamg())
				update_dijk();
		}
}

/*! \fn sideassign (string)
 * \brief  Takes a string input and outputs the associated defined direction
 * \param str a string
 * \return unsigned char
 *
 * Takes a string and compares it to direction strings in order to output the direction assiociated unsigned char
*/
unsigned char sideassign (string str){
	
	if(!str.compare("north"))
		return NORTH;
	else if (!str.compare("south"))
		return SOUTH;
	else if (!str.compare("east"))
		return EAST;
	else if (!str.compare("west"))
		return WEST;
	else if (!str.compare("back"))
		return 0;
	return 1;
}

/*! \fn roamg (void)
 * \brief  Changes the goal cell based on current location.
 * \return int
 *
 * Uses current location to decied on which goal cell to find first.
 * Returns 1 if the goal cell has been updated. Returns 0 if no change.
*/
int roamg (void){
	int old_g, new_g;
	
	old_g = gx + gy*16;

	if( x < 8 && y < 8){
		gx = 7;
		gy = 7;
	}
	else if( x >= 8 && y < 8){
		gx = 8;
		gy = 7;
	}
	else if( x < 8 && y >= 8){
		gx = 7;
		gy = 8;
	}
	else if( x >= 8 && y >= 8){
		gx = 8;
		gy = 8;
	}

	new_g = gx + gy*16;

	if(old_g == new_g)
		return 0;

	return 1;
}

/*! \fn disp_route (void)
 * \brief  Displays the route of the maze to the target cell.
 * \return void
*/
void disp_route (void){
	int u;

	if(seq[0] == 0x00 && seq[1] == 0x00)
		return;

	for(int i = seq_it; i >= 0; i--){
		u = seq[i];
		cout << hex << u;
		if(i)
			cout << ", ";
	}
	cout << endl;
}


