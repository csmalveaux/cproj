#ifndef BINARYTREE_H__
#define BINARYTREE_H__

int sizeoftree(int levels);
int levelsoftree(int size);
int findchildren(int size, int parent, int * child1, int * child2);
int findparent(int size, int child, int * parent);

int addheapuint(unsigned int * arr, int size, unsigned int val);
int addheapint(int * arr, int size, int val);
int addheapfloat(float * arr, int size, float val);
int addheapdouble(double * arr, int size, double val);

int addtreeuint(unsigned int * arr, int size, unsigned int val);
int addtreeint(int * arr, int size, int val);
int addtreefloat(float * arr, int size, float val);
int addtreedouble(double * arr, int size, double val);

#endif /* BINARYTREE_H__ */
