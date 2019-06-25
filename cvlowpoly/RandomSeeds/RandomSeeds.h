#ifndef __RANDOMSEEDS_H__
#define __RANDOMSEEDS_H__

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<iostream>
#include "edlines.h"	

typedef vector<Pixel> EdgeChain;
typedef EdgeChain Seeds;
using namespace std;
#define random(x) (rand()%x)


Seeds creatRandomSeeds(int xsize, int ysize, int NumPoints);


#endif