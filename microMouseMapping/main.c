/*! \file main.cpp
	\brief Micromouse mapping GUI
	\author Chloe Malveaux
	\version 1.5
	\date 2011*/	 

/*! \mainpage μMouse Mapping Documentation
 * Version 1.5
 * \n 07 Mar. 2011
 * \n by Chloe Malveaux
 * \section intro Introduction
 * 
 * An environment to test maze mapping and solving algorithms, before testing it physically on the mouse. Written in C++.
 * \subsection update Updates
 *    - Update 17/02/11: Separate c++ and header files for menu, algorithms, mapping, and data definitions and functions.
 *    - Update 18/02/11: Updated menu.
 *    - Update 04/03/11: Added stack to dijkstra update.
 *    - Update 07/03/11: Added route change check.
 * 
 * \subsection funct Functionality
 *    - Initializes maze and mouse global variables
 *    - Displays maze and mouse information.
 *    - Explore maze using default algorithm, Dijkstra.
 *
 * \subsection functFF Future Functionality
 *    - Display history of commands entered by user.
 *    - Display current status of mouse and maze.
 *    - Change settings (i.e algorithm, goal cell, start cell, search fastest vs shortest path)
 *    - Create and save maze for the mouse to solve
 *    - Save mouse and maze data as a file.
 *    - Load mouse and maze data from file.
 *    - Maintain maze statistics from current solution, such as algorthim used and number of cells visited before target cell found, number of times algorithm updates maze, and "time" till final solution.
 * \section open Open and Running
 * Open zipped file and run the executable file with the C Runtime Libraries in the current directory to ensure the program runs. Run with Windows.
 * \n Zipped file should contain three files:
 *    - uMouse_Mapping.exe
 *    - msvcp100d.dll
 *    - msvcr100d.dll
 *
 * \section pages More Information
 * - \subpage menu Using Program
*/

/*! \page menu Using Program
  \section sec Menu
  - \ref display 
  - \ref clear 
  - \ref help 
  - \ref exit 
  - \ref initialize
  - \ref explore  
  \n Program opens with the options listed above. Commands are entered in the commandline: [ratking]$ (command)
  \n Misspelled or non-existant commands entered will result in the output of "Invalid Option"
  \subsection display
     "[ratking]$ display" or "[ratking]$ disp"
	 \n Displays the maze with the walls explored by the virtual mouse. No additional input necessary.
  \subsection clear
     "[ratking]$ clear" or "[ratking]$ clr"
	 \n Clears the screen, however does not clear any data. No additional input necessary.
  \subsection help
     "[ratking]$ help"
	 \n Displays in-program description of the all commands available. No additional input necessary.
	 \n "[ratking$ help (command)"
	 \n Displays description command specified.
  \subsection exit
     "[ratking]$ exit"
	 \n Exits program.
  \subsection initialize
     "[ratking]$ initialize" or "[ratking]$ init"
	 \n Initializes the maze and mouse data to default values. No additional input necessary.
	 \n Default values for maze are no walls except for exterior walls.
	 \n Default values for mouse are current location is at the start cell and facing east.
  \subsection explore
	 "[ratking]$ explore" or "[ratking]$ explr"
	 \n Allows the user to start traversing the maze by simulating the sensors which would detect walls. Displays maze and current mouse location.
	 \n Acceptable inputs are "north", "south", "east", "west", and "back". No particular order of wall inputs is necessary. Separate by space. Lower case and correct spelling is required. Uppercase or incorrect spelling may result incorrect maze data, and undesired mouse movement.
	 \n In order to return back from the "explore" function, enter "back". 
	 \n When using Dijkstra algorithm, the route the mouse intends to take to the target is displayed as well.
	 \n Ensure that maze is initialized before attempting to traverse the maze.
*/
	  
#include "menu.h"

/*!\fn main (void)
 *	\brief Progam main.
 * Calls menu function
*/
int main (){
	while(menu());
	return 0;
}