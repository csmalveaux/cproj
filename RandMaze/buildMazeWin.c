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
#include <menu.h>

#define ARRAY_SIZE(a) (sizeof(a) / sizeof(a[0]))


/** Prototypes **/
void init_curses(WINDOW * win);
WINDOW * init_mainWindow();
MENU * init_mainMenu(ITEM ** items, WINDOW * win);
ITEM ** init_MenuBarItems();
void free_menu_all(MENU * menu, ITEM ** items);

/** Constants **/

/**
* @fn main
* @brief Main function
**/
void main(int argc, char *argv[]){

	WINDOW * main_win;
	WINDOW * menu_win;

	ITEM ** menu_items;
	MENU * main_menu;

	init_curses(main_win);

	main_win = init_mainWindow();
	menu_win = derwin(main_win, 2, COLS - 2, 1, 1);
	box(menu_win, 0, 0);
	wrefresh(main_win);
	wrefresh(menu_win);
	//menu_items = init_MenuBarItems();
	//main_menu = init_mainMenu(menu_items, main_win);

	getch();

	//free_menu_all(main_menu, menu_items);
	endwin();
}

/**
* @fn init_curses
* @brief Initializes Ncurses
**/
void init_curses(WINDOW * win){
	initscr();
	noecho();
	refresh();

	keypad(win, TRUE);
}

/**
* @fn init_mainWindow
* @brief Initializes main window
**/
WINDOW * init_mainWindow(){
	WINDOW * win;
	win = newwin(LINES, COLS, 0, 0);
	box(win, 0, 0);
	touchwin(win);
	wrefresh(win);

	return win;
}

/**
* @fn init_mainMenu
* @brief Initializes main menu bar
**/
MENU * init_mainMenu(ITEM ** items, WINDOW * win){
	MENU * menu;
	
	menu = new_menu((ITEM **)items);
	set_menu_format(menu, 1, 4);
	set_menu_win(menu, win);
    //set_menu_sub(menu, subwin);
	post_menu(menu);
	wrefresh(win);

	return menu;
}


/**
* @fn init_MenuBarItems
* @brief Initializes menu bar options
**/
ITEM ** init_MenuBarItems(){
	ITEM ** items;

	items[0] = new_item("File", "(F1)");
	items[1] = new_item("Generate", "(F2)");
	items[2] = new_item("Edit", "(F3)");
	items[3] = new_item("View", "(F4)");
	items[4] = (ITEM *) NULL;

	return items;
}

/**
* @fn free_menu_all
* @brief Frees menu and all items
**/
void free_menu_all(MENU * menu, ITEM ** items){
	int i;
	free_menu(menu);

	for(i = 0; i < ARRAY_SIZE(items); i++){
		free_item(items[i]);
	}
}