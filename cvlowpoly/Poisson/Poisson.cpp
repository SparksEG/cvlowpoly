#include "Poisson.h"

float GetDistance(const sPoint& P1, const sPoint& P2) {

	return sqrt((P1.x - P2.x) * (P1.x - P2.x) + (P1.y - P2.y) * (P1.y - P2.y));
}

sGridPoint ImageToGrid(const sPoint& P, float CellSize)			//将离散的散点放到网格点上，CellSize是step？
{
	return sGridPoint((int)(P.x / CellSize), (int)(P.y / CellSize));
}

//初始化，生成一个 m_H * m_W 的网格，其实是一幅图像
sGrid::sGrid(int W, int H, float CellSize)	
	: m_W(W), m_H(H), m_CellSize(CellSize) {

	m_Grid.resize(m_H);
	for (auto i = m_Grid.begin(); i != m_Grid.end(); i++) {

		i->resize(m_W);
	}
}

void sGrid::Insert(const sPoint& P) {	//插入一个散点

	sGridPoint G = ImageToGrid(P, m_CellSize);
	m_Grid[G.x][G.y] = P;
}


//在此点周围为D的正方形内存在一个已经被访问的点A 且 到点A到该点的距离小于MinDist，返回true 
bool sGrid::IsInNeighbourhood(sPoint Point, float MinDist, float CellSize) {		 // point.m_valid 永远为 true

	sGridPoint G = ImageToGrid(Point, CellSize);

	// number of adjucent cells to look for neighbour points      
	const int D = 5;

	// scan the neighbourhood of the point in the grid
	for (int i = G.x - D; i < G.x + D; i++) {

		for (int j = G.y - D; j < G.y + D; j++) {

			if (i >= 0 && i < m_W && j >= 0 && j < m_H) {

				sPoint P = m_Grid[i][j];

				if (P.m_Valid && GetDistance(P, Point) < MinDist) {

					return true;
				}
			}
		}
	}

	return false;
}


sPoint PopRandom(std::vector<sPoint>& Points, DefaultPRNG& Generator)	//从Points中随机剔除一个点
{
	const int Idx = Generator.RandomInt(Points.size() - 1);
	const sPoint P = Points[Idx];
	Points.erase(Points.begin() + Idx);

	return P;
}


sPoint GenerateRandomPointAround(const sPoint& P, float MinDist, DefaultPRNG& Generator)		//环形内随机产生一个附近的点
{
	// start with non-uniform distribution
	float R1 = Generator.RandomFloat();
	float R2 = Generator.RandomFloat();

	// radius should be between MinDist and 2 * MinDist
	float Radius = MinDist * (R1 + 1.0f);

	// random angle
	float Angle = 2 * 3.141592653589f * R2;

	// the new point is generated around the point (x, y)
	float X = P.x + Radius * cos(Angle);
	float Y = P.y + Radius * sin(Angle);

	return sPoint(X, Y);
}

/**
Return a vector of generated points

NewPointsCount - refer to bridson-siggraph07-poissondisk.pdf for details (the value 'k')
Circle  - 'true' to fill a circle, 'false' to fill a rectangle
MinDist - minimal distance estimator, use negative value for default
**/
	//NewPointsCount: 每次增加新的点时，需要向附近查询的次数


