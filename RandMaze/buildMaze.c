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

#define HEIGHT 20
#define WIDTH  

#define ENTER 10
#define ESCAPE 27

void init_curses();
void draw_menubar(WINDOW *menubar);
void delete_menu(ITEM **items,int count);
int scroll_menu(ITEM **items,int count,int menu_start_col);
ITEM **draw_menu(int start_col);

/**
* @fn main
* @brief Main function
**/
int main(int argc, char *argv[]){

    int key;
    WINDOW *menubar, *messagebar;

    init_curses();

    menubar = subwin(stdscr, 1, COLS, 0, 0);
    messagebar = subwin(stdscr, 1, COLS - 1,23,1);
    draw_menubar(menubar);

    refresh();
    
    do {
        int selected_item;
        ITEM **menu_items;
        key = getch();
        werase(messagebar);
        wrefresh(messagebar);
        if (key==KEY_F(1)) {
            menu_items=draw_menu(0);
            selected_item=scroll_menu(menu_items,8,0);
            delete_menu(menu_items,9);
            if (selected_item<0)
                wprintw(messagebar,"You haven't selected any item.");
            else
                wprintw(messagebar,
                  "You have selected menu item %d.",selected_item+1);
            touchwin(stdscr);
            refresh();
        } else if (key==KEY_F(2)) {
            menu_items=draw_menu(20);
            selected_item=scroll_menu(menu_items,8,20);
            delete_menu(menu_items,9);
            if (selected_item<0)
                wprintw(messagebar,"You haven't selected any item.");
            else
                wprintw(messagebar,
                  "You have selected menu item %d.",selected_item+1);
            touchwin(stdscr);
            refresh();
        }
    } while (key!=ESCAPE);

    delwin(menubar);
    delwin(messagebar);
    endwin();
    return 0;

}

void init_curses(){
    initscr();
    curs_set(0);
    noecho();
    keypad(stdscr,TRUE);
}

void draw_menubar(WINDOW *menubar)
{
    waddstr(menubar,"File");
    waddstr(menubar,"(F1)");
    wmove(menubar,0,20);
    wattron(menubar,WA_DIM);
    waddstr(menubar,"Generate");
    waddstr(menubar,"(F2)");
    wattroff(menubar,WA_DIM);
    wmove(menubar,0,40);
    wattron(menubar,WA_DIM);
    waddstr(menubar,"View");
    waddstr(menubar,"(F3)");
    wattroff(menubar,WA_DIM);
    wmove(menubar,0,60);
    wattron(menubar,WA_DIM);
    waddstr(menubar,"Edit");
    waddstr(menubar,"(F4)");
    wattroff(menubar,WA_DIM);
}

ITEM ** draw_filebar(){
    ITEM **items;
    items = (ITEM **) malloc(6 * sizeof(ITEM *));

    items[0] = new_item("","");
}

ITEM **draw_menu(int start_col)
{
    int i;
    ITEM **items;
    items=(ITEM **)malloc(9*sizeof(ITEM *));

    items[0]=newwin(10,19,1,start_col);
    box(items[0],ACS_VLINE,ACS_HLINE);
    items[1] = new_item("New","New");//subwin(items[0],1,17,2,start_col+1);
    //wprintw(items[1],"New");
    items[2] = new_item("Open","Open"); //subwin(items[0],1,17,3,start_col+1);
    //wprintw(items[2],"Open");
    items[3] = new_item("Save","Save"); //subwin(items[0],1,17,4,start_col+1);
    //wprintw(items[3],"Save");
    items[4] = new_item("Save As","Save As"); //subwin(items[0],1,17,5,start_col+1);
    //wprintw(items[4],"Save As");
    items[5] = new_item("Print","Print"); //subwin(items[0],1,17,6,start_col+1);
    //wprintw(items[5],"Print");

    wrefresh(items[0]);
    return items;
}

void delete_menu(ITEM **items,int count)
{
    int i;
    for (i=0;i<count;i++)
        delwin(items[i]);
    free(items);
}

int scroll_menu(ITEM **items,int count,int menu_start_col)
{
    int key;
    int selected=0;
    while (1) {
        key=getch();
        if (key==KEY_DOWN || key==KEY_UP) {
            //wbkgd(items[selected+1],COLOR_PAIR(2));
            wnoutrefresh(items[selected+1]);
            if (key==KEY_DOWN) {
                selected=(selected+1) % count;
            } else {
                selected=(selected+count-1) % count;
            }
            wbkgd(items[selected+1],COLOR_PAIR(1));
            wnoutrefresh(items[selected+1]);
            doupdate();
        } else if (key==KEY_LEFT || key==KEY_RIGHT) {
            delete_menu(items,count+1);
            touchwin(stdscr);
            refresh();
            items=draw_menu(20-menu_start_col);
            return scroll_menu(items,8,20-menu_start_col);
        } else if (key==ESCAPE) {
            return -1;
        } else if (key==ENTER) {
            return selected;
        }
    }
}