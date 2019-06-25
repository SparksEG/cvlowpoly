/*=================================================================================================
Author: Renato Farias (renatomdf@gmail.com)
Created on: April 13th, 2012
Updated on: June 12th, 2012
About: This is a simple CPU implementation of the Jump Flooding algorithm from the paper "Jump
Flooding in GPU With Applications to Voronoi Diagram and Distance Transform" [Rong 2006]. The
result is a Voronoi diagram generated from a number of seeds which the user provides with mouse
clicks. You can also click on and drag around a seed to reposition it, if that's your thing.
=================================================================================================*/

#ifndef __JFA_H__
#define __JFA_H__
/*=================================================================================================
INCLUDES
=================================================================================================*/
#include "ImgLib.h"
#include "edlines.h"	
#include <vector>
#include <fstream>
#include <unordered_map>


typedef vector<Pixel> EdgeChain;
typedef EdgeChain Seeds;
//typedef vector<vector<Pixel> >  VecEdgeChains;

typedef unordered_multimap<int, int> VoiPixelMap;
//namespace JFA {
//
//	//这里因为命名相同无法编译通过的原因，简单通过设置命名空间来解决，后续需要完善
//
//	//using namespace std;
//	typedef struct Pixel {
//		int x, y;
//
//		Pixel(int tx, int ty) {
//			x = tx;
//			y = ty;
//		}
//		Pixel() {};
//	} Pixel;
//
//	typedef std::vector<Pixel> EdgeChain;
//}




// Represents a Pixel with (x,y) coordinates



/*=================================================================================================
FUNCTIONS
=================================================================================================*/
//void ExecuteJumpFlooding(vector<Pixel> Seeds);

void ExecuteJumpFloodingDis(float* flowDisMap, Seeds Seeds, unsigned int xsize, unsigned int ysize);

void ExecuteJumpFloodingVoi(int* VoiMap, const Seeds Seeds, unsigned int xsize, unsigned int ysize);

//EdgeChain VertexOptimization(const EdgeChain constrainedPoints, const EdgeChain cvtSeeds,
//	unsigned int xsize, unsigned int ysize, int optTimes, int* distMap);

Seeds VertexOptimization(const Seeds cvtSeeds, unsigned int xsize, unsigned int ysize, int optTimes, float* flowDistMap);

//EdgeChains VertexOptimization(const EdgeChains constrainedPoints, const EdgeChains cvtSeeds,
//	unsigned int xsize, unsigned int ysize, int optTimes, float* distMap);
#endif