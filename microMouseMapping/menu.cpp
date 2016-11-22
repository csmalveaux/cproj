/*!\file menu.cpp
 * \brief Menu C++ File
 *
 * Contains functions for menu display
*/

#include "menu.h"
#include "mapping.h"

/*! \var extern int mode
	\brief  Menu mode selection.*/
int mode = CLEAR;
/*! \var extern int first
 *  \brief  Used to indicate the first time running through explore mode.*/
int first = 0;

/*!\fn menu (void)
 * \brief  Menu function
 * \return int 1 or int 0 
 *
 * Continuous menu loop returns 1 as defualt and returns 0 when user exits the program. 
*/
int menu (void){
	string str, str1, str2;
	size_t it;
	int tx, ty, tm;

	// if the screen is cleared the menu is displayed
	if(mode != EXIT && mode != INIT && mode != HELP && mode != NOINST)
		commdisp ();

	// if maze is to be displayed
	if(mode != EXIT && mode != INIT && mode != CLEAR && mode != HELP && mode != NOINST)
		mazedisp ();
	
	// if not in explore mode then menu options can be entered
	if(mode != EXPLR){
		cout << "[ratking]$ ";
		str.clear();
		str1.clear();
		str2.clear();
		getline (cin, str2);
		while(!str2.empty()){
			it = str2.find_first_of(" ");
			str1.assign(str2, 0, it);
			if(!str1.compare("help"))
				str.assign(str1);
			str2.erase (0, it + 1);
			if(!str2.compare(str1))
				str2.clear();
		}

		if(str1.compare("help") && str.compare("help") && str.empty()){
			str.assign(str1);
			str1.clear();
		}

		mode = commandsel (str);
	}

	// menu switch
	switch(mode){
	
	// exits program 
	case EXIT:
		return mode;
		break;
	
	// initalizes the maze walls 
	case INIT:
		init();
		cout.width(20); cout << left << "Maze Initialized." << endl;
		first = 0;
		break;

	// enters explore mode
	case EXPLR:
		mode = EXPLR;
		
		if (first){
			if(wsensors ()){
				tx = x;
				ty = y;
				tm = mode;
				mwupdate();
				x = tx;
				y = ty;
				mode = tm;
				update();
				mgoto ();
			}
		}

		first = 1;
		break;
	
	// enters solve mode, not yet implemented
	case SOLVE:
		cout.width(20); cout << left << "Function not yet available." << endl;
		break;
	case DISP: break;
	case CLEAR: break;
	case SETT:
		cout.width(20); cout << left << "Function not yet available." << endl;
		break;
	case HELP:
		helpdisp (str1);
		break;
	case SAVE:
		cout.width(20); cout << left << "Function not yet available." << endl;
		break;
	case LOAD:
		cout.width(20); cout << left << "Function not yet available." << endl;
		break;
	case STAT:
		cout.width(20); cout << left << "Function not yet available." << endl;
		break;
	case STATIS:
		cout.width(20); cout << left << "Function not yet available." << endl;
		break;
	default:
		cout.width(20); cout << left << "Invalid Option" << endl;
	}

	if (mode != EXIT && mode != INIT && mode != HELP && mode != NOINST)
		system ("cls");
	return 1;
}

/*! \fn commandsel (string str)
 * \brief  Takes a string command input and outputs the int associated with the command
 * \param str
 * \return int
 *
 * Takes a string and compares it to possible command input strings in order to output associated int command. Lower case only.
*/
int commandsel (string str){
	
	if(!str.compare("init")||!str.compare("initialize")||!str.compare("1"))
		return INIT;
	else if (!str.compare("disp")||!str.compare("display")||!str.compare("2"))
		return DISP;
	else if (!str.compare("explr")||!str.compare("explore")||!str.compare("3"))
		return EXPLR;
	else if (!str.compare("solve")||!str.compare("4"))
		return SOLVE;
	else if (!str.compare("clr")||!str.compare("clear")||!str.compare("5"))
		return CLEAR;
	else if (!str.compare("set")||!str.compare("settings")||!str.compare("6"))
		return SETT;
	else if (!str.compare("help"))
		return HELP;
	else if (!str.compare("exit")||!str.compare("0"))
		return EXIT;
	else
		return NOINST;
}

/*! \fn commdisp (void)
 * \brief  Displays the commands
 * \return void
 *
 * Function displays the options available to the user.
 * 1) Initialize 2) Display maze 3) Explore Maze 4) Solve Maze 5) Clear screen 0) Exit
*/
void commdisp (void){
	cout << endl;
	cout.width(48); cout << internal << "Micromouse Maze Mapping GUI\n\n" << endl;
	cout.width(20); cout << internal << "Menu Functions";	cout.width(20); cout << internal << "Maze Functions";	cout.width(20); cout << internal << "Data Functions\n" << endl;
	cout.width(20); cout << left << "\t display";			cout.width(20); cout <<  left << "  initialize";		cout.width(20); cout <<  left << "  settings*" << endl;
	cout.width(20); cout << left << "\t history*";			cout.width(20); cout <<  left << "  explore";			cout.width(20); cout <<  left << "  save*" << endl;
	cout.width(20); cout << left << "\t clear";				cout.width(20); cout <<  left << "  solve*";			cout.width(20); cout <<  left << "  load*" << endl;
	cout.width(20); cout << left << "\t help";				cout.width(20); cout <<  left << "  status*";			cout.width(20); cout <<  left << "  statistics*" << endl;
	cout.width(20); cout << left << "\t exit\n" << endl;
	cout << " 'help [command]' displays information about command.\n If command is not specified, information on all supported commands are listed.\n" << endl;
	cout << " *Indicates function not yet avaliable.\n" << endl;
}

/*! \fn mazedisp (void)
 * \brief  Displays the maze
 * \return void
 *
 * Function displays the maze. "*" denotes visited, @ denotes current location of mouse
*/
void mazedisp (void){

	string vline, hline, cnhline, cshline, bnhline;
	unsigned char currc, leftc, rightc, topc, bottc;
	
	for(int i = 0; i < 16; i++){
			
			vline.clear();
			hline.clear();
			cnhline.clear();
			cshline.clear();
			bnhline.clear();

			hline.append(" ");
			cnhline.append(" ");
			cshline.append(" ");
			bnhline.append(" ");
			
			for(int k = 0; k < 16; k++){
				
				//currpos = x + y*16;
				int pos = k+i*16;
				
				// initializing temp cell data
				currc = mmap[pos];
				leftc = 0x00;
				rightc= 0x00;
				topc = 0x00;
				bottc = 0x00;
				
				// saving cells around current cell
				if(k != 0)				// if west side cells
					leftc = mmap[pos-1];
				else
					leftc |= EAST;
				if(k != 15)				// if east side cells
					rightc = mmap[pos+1];
				else
					rightc |= WEST;
				if(i == 0)				// if north side cells
					topc = mmap[pos-16];
				else
					topc |= SOUTH;
				if(i == 15)				// if south side cells
					bottc = mmap[pos+16];
				else
					bottc |= SOUTH;

				// setting cnhline string
				if((currc&NORTH)==NORTH)
					cnhline.append("--- ");
				else
					cnhline.append("    ");

				// setting cshline string
				if((currc&SOUTH)==SOUTH)
					cshline.append("--- ");
				else
					cshline.append("    ");

				// setting bnhline string
				if((bottc&NORTH)==NORTH)
					bnhline.append("--- ");
				else
					bnhline.append("    ");

				// prints | for west side of maze
				if(k == 0)
					if((currc&WEST)==WEST)
						vline.append("|");
					else
						vline.append(" ");

				// prints asterisk if cell visited, else whitespace, @ if mouse is currently in the cell
				if(pos==(x + y*16))
					vline.append(" @ ");
				else if((currc&VISIT)==VISIT)
					vline.append(" * ");
				else
					vline.append("   ");

				// prints | if cell has west wall or right cell has east wall, else whitespace
				if((currc&EAST)==EAST||(rightc&WEST)==WEST)
					vline.append("|");
				else
					vline.append(" ");
			}
		
			// checking the current row and bottom row for shared walls
			for(unsigned int g = 0; g < cshline.length(); g++)
				if((cshline.compare(g,1,"-") == 0)||(bnhline.compare(g,1,"-") == 0))
					hline.append("-");
				else
					hline.append(" ");
			// displaying strings
			if(i==0){
				cout.width(70); cout << internal << cnhline << endl; // displays north wall
			}
			cout.width(70); cout << internal << vline << endl;		 // displays vertical walls
			cout.width(70); cout << internal << hline << endl;		 // displays horizonal walls
	}
}

/*! \fn helpdisp (string)
 * \brief  Displays the help
 * \param str
 * \return void
 *
 * Function displays the command help information
*/
void helpdisp (string str){
	if (!str.compare("help")){
		cout << "HELP commands\n" << endl;
		cout << "Menu functions." << endl;
	}

	if (!str.compare("disp")||!str.compare("display")||!str.compare("help")){
		cout.width(10); cout << left << "\tdisplay";
		cout.width(30); cout << left << "\t- Prints maze and current location of mouse." << endl;
	}
	if (!str.compare("history")||!str.compare("help")){
		cout.width(10); cout << left << "\thistory";
		cout.width(30); cout << left << "\t- Displays history of commands entered." << endl;
	}
	if (!str.compare("clr")||!str.compare("clear")||!str.compare("help")){
		cout.width(10); cout << left << "\tclear";
		cout.width(30); cout << left << "\t- Clears the screen execpt for menu." << endl;
	}
	if (!str.compare("help")||!str.compare("help")){
		cout.width(10); cout << left << "\thelp";
		cout.width(30); cout << left << "\t- Displays information about specified command." << endl;
		cout.width(10); cout << left << "\t";
		cout.width(30); cout << left << "\t  Entering 'help' within a mode, will present more" << endl;
		cout.width(10); cout << left << "\t";
		cout.width(30); cout << left << "\t  detailed information of mode." << endl;
	}
	if (!str.compare("exit")||!str.compare("help")){
		cout.width(10); cout << left << "\texit";
		cout.width(30); cout << left << "\t- Closes and exit program." << endl;
	}
	
	if (!str.compare("help")){
		cout << endl;
		cout << "Maze functions." << endl;
	}
	if (!str.compare("initialize")||!str.compare("help")){
		cout.width(10); cout << left << "\tinitialize";
		cout.width(30); cout << left << "\t- Resets map data and start and goal cells." << endl;
		cout.width(10); cout << left << "\t";
		cout.width(30); cout << left << "\t  To change start and goal cells from default, enter " << endl;
		cout.width(10); cout << left << "\t";
		cout.width(30); cout << left << "\t  'settings'." << endl;
	}
	if (!str.compare("explr")||!str.compare("explore")||!str.compare("help")){
		cout.width(10); cout << left << "\texplore";
		cout.width(30); cout << left << "\t- Moves mouse through maze." << endl;
	}
	if (!str.compare("solve")||!str.compare("help")){
		cout.width(10); cout << left << "\tsolve";
		cout.width(30); cout << left << "\t- Calculates routes from start to go cell." << endl;
	}
	if (!str.compare("status")||!str.compare("help")){
		cout.width(10); cout << left << "\tstatus";
		cout.width(30); cout << left << "\t- Displays current status of maze." << endl;
		cout.width(10); cout << left << "\t";
		cout.width(30); cout << left << "\t  i.e. % maze explored, if goal found, algorithm used..." << endl;
	}
	
	if (!str.compare("help")){
		cout << endl;
		cout << "Data functions." << endl;
	}
	if (!str.compare("set")||!str.compare("settings")||!str.compare("help")){
		cout.width(10); cout << left << "\tsettings";
		cout.width(30); cout << left << "\t- Allows user to change or reset default settings." << endl;
		cout.width(10); cout << left << "\t";
		cout.width(30); cout << left << "\t  Settings include start and goal cell, and algorithm." << endl;
	}
	if (!str.compare("save")||!str.compare("help")){
		cout.width(10); cout << left << "\tsave";
		cout.width(30); cout << left << "\t- Save current maze and status." << endl;
	}
	if (!str.compare("load")||!str.compare("help")){
		cout.width(10); cout << left << "\tload";
		cout.width(30); cout << left << "\t- Loads previously saved mazes." << endl;
	}
	if (!str.compare("statistics")||!str.compare("help")){
		cout.width(10); cout << left << "\tstatistics";
		cout.width(30); cout << left << "\t- For loaded maze, compares algorithm solutions." << endl;
	}
}

