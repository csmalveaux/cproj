#include "mouse.h"
#include "map.h"

unsigned int mouse_pos = -1;
unsigned int mouse_dir = -1;

void intialize_mouse(struct mouse * ms, unsigned int start, unsigned int goal);
void setPerimeter(unsigned short * map, const int sz, const int width);
void pushNewPosition(struct mouse * ms, unsigned int next);
void goBack(struct mouse * ms);

void intialize_mouse(struct mouse * ms, unsigned int start, unsigned int goal){
	ms->start = start;
	ms->goal = goal;
	ms->curr = start;
	ms->prev = 0;
	ms->next = 0;
	ms->dir = 0;
}

void setPerimeter(unsigned short * map, const int sz, const int width){
	int i;
	for(i = 0; i < sz; i++){
		map[i] = 0x00;
		map[i] |= VISIT;

		if(i < width)
			map[i] |= NORTH;

		if(i % width == 0)
			map[i] |= WEST;

		if(i % width == (width - 1))
			map[i] |= EAST;

		if(i >= (sz - width))
			map[i] |= SOUTH;
	}
}

void pushNewPosition(struct mouse * ms, unsigned int next){
	ms->prev = ms->curr;
	ms->curr = next;
	mouse_pos = ms->curr;
	ms->dir = ms->prev - ms->curr;
	mouse_dir = ms->dir;
}

void goBack(struct mouse * ms){
	unsigned int temp;
	temp = ms->start;
	ms->start = ms->goal;
	ms->goal = temp;
}
