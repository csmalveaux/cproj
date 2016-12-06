#ifndef MOUSE_H
#define MOUSE_H

struct mouse{
	unsigned int start, next, curr, prev, goal;
	int dir;
	unsigned short * map;
};

extern unsigned int mouse_pos;
extern unsigned int mouse_dir;

void intialize_mouse(struct mouse * ms, unsigned int start, unsigned int goal);
void pushNewPosition(struct mouse * ms, unsigned int next);
void setPerimeter(unsigned short * map, const int sz, const int width);
void goBack(struct mouse * ms);

#endif // MOUSE_H