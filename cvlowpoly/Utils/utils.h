#ifndef __UTILS__
#define __UTILS__

#include "PathCompacter.h"
#include "edlines.h"
#include "ImgLib.h"
#include "JFA.h"

typedef vector<Pixel> EdgeChain;
typedef vector<vector<Pixel> > VecEdgeChains;

/****************************
功能：将 算法edlines的数据结构 【EdgeChains】 转换成 散点 的数据结构【EdgeChain】
输入：链式vec存储的所有图像边缘 edges
输出：单vec存储的所有图像边缘
返回：是否成功
****************************/
int EdgeChainsToSeeds(EdgeChain &seeds, const EdgeChains edges);

/****************************
功能：将 算法edlines的数据结构 【EdgeChains】 转换成  vec<vec<Pixel> > 的数据结构【VecEdgeChains】
输入：链式vec存储的所有图像边缘 edges
输出：双vec存储的所有图像边缘
返回：是否成功
****************************/
int EdgeChainsToVecEdgeChains(VecEdgeChains &vedges, const EdgeChains edges);


/****************************
功能：将 vec<vec> Pixel 的数据结构【VecEdgeChains】 转换成  算法edlines的数据结构 【EdgeChains】
输入：双vec存储的所有图像边缘
输出：链式vec存储的所有图像边缘 edges
返回：是否成功
****************************/
int VecEdgeChainsToEdgeChains(EdgeChains &edges, const VecEdgeChains vedges);


/****************************
功能：将 vec<Pixel>  的数据结构【EdgeChain】  转换成  算法PathCompact的数据结构 【DVector2D*】
输入：vec存储的一条图像边缘
返回：数组指针存储的一条图像边缘
****************************/
DVector2D* EDToPath(EdgeChain edge);


/****************************
功能：将 算法PathCompact的数据结构【DVector2D*】 转换成  vec<Pixel> 的数据结构【EdgeChain】
输入：数组指针存储的一条图像边缘edge，该边缘的像素点个数
返回：vec存储的一条图像边缘
****************************/
EdgeChain PathToED(DVector2D* segment, unsigned int uPointsInCurrentPath);


/****************************
功能：将两点及其经过他们中间线段上的点找到，并更新添加到给的容器【EdgeChains】中
输入：两点p1，p2
返回：更新后的容器【EdgeChains】
****************************/
void PixelsLinkToEdgeChain(EdgeChain& edge, Pixel p1, Pixel p2);


/****************************
功能：将 算法edlines的数据结构 【EdgeChains】 转换成  算法JFA的数据结构 【EdgeChain】
输入：链式vec存储的所有图像边缘 edges
返回：vec存储的 所有边缘上的散点
****************************/
void EdgeChainsToJFASeeds(EdgeChain& Seeds, const EdgeChains& edges);


/****************************
功能：将所有的边缘上的点画成点
输入：所有图像边缘 edges，图像长xsize 宽ysize
****************************/
int DrawPlots(EdgeChains edges, unsigned int xsize,
	unsigned int ysize);


/****************************
功能：将所有的边缘画成线
输入：所有图像边缘 edges，图像宽xsize 长ysize
****************************/
int DrawLines(EdgeChains edges, unsigned int xsize,
	unsigned int ysize);


#endif
