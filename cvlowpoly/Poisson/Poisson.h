
#ifndef __POISSON_H__
#define __POISSON_H__

#include <vector>
#include <random>
#include <stdint.h>
#include <time.h>

class DefaultPRNG		//α�������
{
public:
	DefaultPRNG()
		: m_Gen(std::random_device()()), m_Dis(0.0f, 1.0f)	//random_device: �������
	{
		// prepare PRNG
		m_Gen.seed(time(nullptr));
	}

	explicit DefaultPRNG(uint32_t seed)
		: m_Gen(seed), m_Dis(0.0f, 1.0f)
	{
	}

	float RandomFloat()
	{
		return static_cast<float>(m_Dis(m_Gen));
	}

	int RandomInt(int Max)
	{
		std::uniform_int_distribution<> DisInt(0, Max);
		return DisInt(m_Gen);
	}

private:
	std::mt19937 m_Gen;								//����һ��α�����
	std::uniform_real_distribution<float> m_Dis;	//���������ȷֲ�
};

struct sPoint			//���ɵĵ�
{
	sPoint(): x(0), y(0), m_Valid(false) {}
	sPoint(float X, float Y): x(X), y(Y), m_Valid(true) {}
	float x;
	float y;
	bool m_Valid;

	bool IsInRectangle() const {

		return x >= 0 && y >= 0 && x <= 1 && y <= 1;
	}

	bool IsInCircle() const {

		float fx = x - 0.5f;
		float fy = y - 0.5f;
		return (fx*fx + fy * fy) <= 0.25f;
	}
};

struct sGridPoint		//�����
{
	sGridPoint(int X, int Y): x(X), y(Y) {}
	int x;
	int y;
};

struct sGrid
{
	sGrid(int W, int H, float CellSize);		//��ʼ��������һ�� m_H * m_W ��������ʵ��һ��ͼ��

	void Insert(const sPoint& P);

	//�ڴ˵���ΧΪD���������ڴ���һ���Ѿ������ʵĵ�A �� ����A���õ�ľ���С��MinDist������true 
	bool IsInNeighbourhood(sPoint Point, float MinDist, float CellSize);


private:
	int m_W;
	int m_H;
	float m_CellSize;

	std::vector< std::vector<sPoint> > m_Grid;
};

template <typename PRNG>
sPoint PopRandom(std::vector<sPoint>& Points, PRNG& Generator)	//��Points������޳�һ����
{
	const int Idx = Generator.RandomInt(Points.size() - 1);
	const sPoint P = Points[Idx];
	Points.erase(Points.begin() + Idx);

	return P;
}

template <typename PRNG>
sPoint GenerateRandomPointAround(const sPoint& P, float MinDist, PRNG& Generator)		//�������������һ�������ĵ�
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

template <typename PRNG>
sPoint PopRandom(std::vector<sPoint>& Points, PRNG& Generator);	//��Points������޳�һ����

template <typename PRNG>
sPoint GenerateRandomPointAround(const sPoint& P, float MinDist, PRNG& Generator);		//�������������һ�������ĵ�

/**
Return a vector of generated points

NewPointsCount - refer to bridson-siggraph07-poissondisk.pdf for details (the value 'k')
Circle  - 'true' to fill a circle, 'false' to fill a rectangle
MinDist - minimal distance estimator, use negative value for default
**/

template <typename PRNG = DefaultPRNG>		//NewPointsCount: ÿ�������µĵ�ʱ����Ҫ�򸽽���ѯ�Ĵ���
std::vector<sPoint> GeneratePoissonPoints(size_t NumPoints, DefaultPRNG& Generator, int NewPointsCount = 10,
	bool Circle = false, float MinDist = -1.0f) {
	if (MinDist < 0.0f) {

		MinDist = sqrt(float(NumPoints)) / float(NumPoints);
	}

	std::vector<sPoint> SamplePoints;
	std::vector<sPoint> ProcessList;

	// create the grid
	float CellSize = MinDist / sqrt(2.0f);

	int GridW = (int)ceil(1.0f / CellSize);
	int GridH = (int)ceil(1.0f / CellSize);

	//sGrid Grid(GridW, GridH, CellSize);
	GridW = 32;
	GridH = 4096;
	sGrid Grid(GridW, GridH, CellSize);

	sPoint FirstPoint;
	do {
		FirstPoint = sPoint(Generator.RandomFloat(), Generator.RandomFloat());
	} while (!(Circle ? FirstPoint.IsInCircle() : FirstPoint.IsInRectangle()));

	// update containers
	ProcessList.push_back(FirstPoint);
	SamplePoints.push_back(FirstPoint);
	Grid.Insert(FirstPoint);

	// generate new points for each point in the queue
	while (!ProcessList.empty() && SamplePoints.size() < NumPoints) {

		sPoint Point = PopRandom(ProcessList, Generator);

		for (int i = 0; i < NewPointsCount; i++) {

			sPoint NewPoint = GenerateRandomPointAround(Point, MinDist, Generator);
			bool Fits = Circle ? NewPoint.IsInCircle() : NewPoint.IsInRectangle();
			if (Fits && !Grid.IsInNeighbourhood(NewPoint, MinDist, CellSize)) {

				ProcessList.push_back(NewPoint);
				SamplePoints.push_back(NewPoint);
				Grid.Insert(NewPoint);
				continue;
			}
		}
	}

	return SamplePoints;
}
float GetDistance(const sPoint& P1, const sPoint& P2);
sGridPoint ImageToGrid(const sPoint& P, float CellSize);

#endif