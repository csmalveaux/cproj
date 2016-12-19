/****************************************************************************
 * Copyright (C) 2016 by Chloe Malveaux                                     *
 *                                                                          *
 * This file is part of Maze Explorer.                                      *
 *                                                                          *
 *   Maze Explorer is free software: you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as         *
 *   published by the Free Software Foundation, either version 3 of the     *
 *   License, or (at your option) any later version.                        *
 *                                                                          *
 *   Maze Explorer is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *   GNU Lesser General Public License for more details.                    *
 *                                                                          *
 *   You should have received a copy of the GNU Lesser General Public       *
 *   License along with Maze Explorer.                                      *
 *   If not, see <http://www.gnu.org/licenses/>.                            *
 ****************************************************************************/

/**
* @file buildMaze.c
* @author Chloe Malveaux
* @date Dec 18 2016
* @brief Displays and saves randomly generated mazes.
* 
* Through terminal user interface, selects maze size (height, width), 
* generation algorithm, display options, and start/goal cells to create
* randomly generated maze. User can choose to save and name maze created.
**/

#include <ncurses.h>

#define HEIGHT 50
#define WIDTH  50

/**
* @fn main
* @brief Main function
**/
int main(int argc, char *argv[]){
	
	WINDOW * mainWin;
	int offsetx, offsety;

	initscr();
    refresh();

    offsetx = (COLS - WIDTH) / 2;
    offsety = (LINES - HEIGHT) / 2;

    mainWin = newwin(HEIGHT, WIDTH, offsety, offsetx);
    box(mainWin, 0, 0);
    wrefresh(mainWin);
    getch();
    delwin(mainWin);
    endwin();

    return 0;

}