#include "utils.h"
#include "PathCompacter.h"
#include <vector>
#include <iostream>
using namespace std;

//struct EdgeChains {
//	std::vector<unsigned int> xCors;//all the x coordinates of edge points
//	std::vector<unsigned int> yCors;//all the y coordinates of edge points
//	std::vector<unsigned int> sId;  //the start index of each edge in the coordinate arrays
//	unsigned int numOfEdges;//the number of edges whose length are larger than minLineLen; numOfEdges < sId.size;
//};

// The following struct represents a unsigned int precision 2D point.
//struct Pixel {
//	unsigned int x;//X coordinate
//	unsigned int y;//Y coordinate
//};
//typedef vector<Pixel> EdgeChain;
//typedef vector<vector<Pixel> > VecEdgeChains;

//int EdgeChainsToVecEdgeChains(VecEdgeChains &vedges, const EdgeChains edges) {
//
//	Pixel dot;
//	vedges.resize(edges.numOfEdges);
//	unsigned int index = 0;
//	for (int i = 0; i < edges.numOfEdges; i++) {
//		for (int k = 0; k < (edges.sId[i + 1] - edges.sId[i]) && index < edges.xCors.size(); k++) {
//			dot.x = edges.xCors[index];
//			dot.y = edges.yCors[index];
//			vedges[i].push_back(dot);
//			index++;
//		}
//	}
//	return 1;
//
//}
//
//int VecEdgeChainsToEdgeChains(EdgeChains &edges, const VecEdgeChains vedges) {
//
//	edges.sId.clear();
//	edges.xCors.clear();
//	edges.yCors.clear();
//	edges.numOfEdges = vedges.size();
//	unsigned int index = 0;
//	edges.sId.push_back(index);
//	for (int i = 0; i < vedges.size(); i++) {
//		for (int j = 0; j < vedges[i].size(); j++) {
//			edges.xCors.push_back(vedges[i][j].x);
//			edges.yCors.push_back(vedges[i][j].y);
//			index++;
//		}
//		edges.sId.push_back(index);
//	}
//	return 1;
//
//}
//
//DVector2D* EDToPath(EdgeChain edge) {
//
//	size_t num = edge.size();
//	DVector2D* segment = new DVector2D[num];
//	for (int i = 0; i < num; i++) {
//		segment[i].dX = edge[i].x;
//		segment[i].dY = edge[i].y;
//	}
//	return segment;
//
//}
//
//EdgeChain PathToED(DVector2D* segment, unsigned int uPointsInCurrentPath) {
//
//	unsigned int num = uPointsInCurrentPath;
//	EdgeChain edge(num);
//	for (int i = 0; i < num; i++) {
//		edge[i].x = segment[i].dX;
//		edge[i].y = segment[i].dY;
//	}
//	return edge;
//
//}

//static double perpendicularDistance(DVector2D start, DVector2D end, DVector2D mid,
//	double dSquareSegmentLength)
//{
//	double dArea;
//
//	dArea = start.dX * (mid.dY - end.dY) +
//		mid.dX * (end.dY - start.dY) +
//		end.dX * (start.dY - mid.dY);
//	return dArea * dArea / dSquareSegmentLength;
//}



//int main() {
//	//DeviationMetric perpendicularDistanceDeviationMetric = &perpendicularDistance;
//	DeviationMetric perpendicularDistanceDeviationMetric = &shortestDistanceToSegment;
//
//
//
//	EdgeChains edges;
//	edges.yCors = { 14,14,14,14,14,14,14,14,14,14,14,14,14,14 };
//	edges.xCors = { 15,14,13,12,11,10,9,8,7,6,5,4,3,2 };
//	edges.numOfEdges = 1;
//	edges.sId = { 0,16 };
//
//	VecEdgeChains vedges;
//	EdgeChainsToVecEdgeChains(vedges, edges);
//	VecEdgeChainsToEdgeChains(edges, vedges);
//	for (int i = 0; i < vedges.size(); i++) {
//		DVector2D *segment = EDToPath(vedges[i]);
//		unsigned int upointsincurrentpath = vedges[i].size();	
//		compactPath(segment, upointsincurrentpath, segment, &upointsincurrentpath, 0.02, perpendicularDistanceDeviationMetric);
//		vedges[i] = PathToED(segment, upointsincurrentpath);
//		delete segment;
//	}
//	VecEdgeChainsToEdgeChains(edges, vedges);
//	int a = 1;
//
//}
