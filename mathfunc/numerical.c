#include "numerical.h"

#define A	7
#define B	5
#define M 	11

static int Xn;

int random(int seed);

int random(int seed){
	int Xn1 = (A*Xn + B)%M;
	Xn = Xn1;
	return Xn1;
}