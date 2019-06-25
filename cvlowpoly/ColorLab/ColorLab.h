#ifndef __COLOR_LAB_H_
#define __COLOR_LAB_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <unordered_map>
using namespace std;


#include <opencv2/core/utility.hpp>  
#include <opencv2/opencv.hpp>

using namespace cv;

//�ඨ�壺��ά����
class Vector2d
{
public:
	double x_;
	double y_;

public:
	Vector2d(double x, double y) :x_(x), y_(y) {}
	Vector2d() :x_(0), y_(0) {}

	//��ά�������, ��˵Ľ����ʵ������������ֱ������������ɵ�ƽ�棬��������ֻ��Ҫ���С�ͷ���
	double CrossProduct(const Vector2d vec)
	{
		return x_ * vec.y_ - y_ * vec.x_;
	}

	//��ά�������
	double DotProduct(const Vector2d vec)
	{
		return x_ * vec.x_ + y_ * vec.y_;
	}

	//��ά��������
	Vector2d Minus(const Vector2d vec) const
	{
		return Vector2d(x_ - vec.x_, y_ - vec.y_);
	}
};

//��������
class Triangle
{
public:
	Vector2d pointA_, pointB_, pointC_;

public:
	Triangle(Vector2d point1, Vector2d point2, Vector2d point3)
		:pointA_(point1), pointB_(point2), pointC_(point3)
	{
		//todo �ж������Ƿ���
	}

	bool IsPointInTriangle(const Vector2d pointP)
	{
		Vector2d PA = pointA_.Minus(pointP);
		Vector2d PB = pointB_.Minus(pointP);
		Vector2d PC = pointC_.Minus(pointP);
		double t1 = PA.CrossProduct(PB);
		double t2 = PB.CrossProduct(PC);
		double t3 = PC.CrossProduct(PA);
		return t1 * t2 >= 0 && t1*t3 >= 0;
	}
};

//void reloadsort();

bool IsPointInTriangle(const Vector2d p, const Triangle t);
void ColorPointInTriangle(vector<vector<bool>> IsPointScanedMap, const Mat& LabMap, vector<vector<Vec3f>>& ColorMap, const Vec6f& triangle);
void ReadImgColorOnOpencv(const Mat img, vector<vector<Vec3f>>& ColorMap);
void WriteImgColorOnOpencv(const vector<vector<Vec3f>> ColorMap, Mat& img);


#endif