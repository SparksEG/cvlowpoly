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

//类定义：二维向量
class Vector2d
{
public:
	double x_;
	double y_;

public:
	Vector2d(double x, double y) :x_(x), y_(y) {}
	Vector2d() :x_(0), y_(0) {}

	//二维向量叉乘, 叉乘的结果其实是向量，方向垂直于两个向量组成的平面，这里我们只需要其大小和方向
	double CrossProduct(const Vector2d vec)
	{
		return x_ * vec.y_ - y_ * vec.x_;
	}

	//二维向量点积
	double DotProduct(const Vector2d vec)
	{
		return x_ * vec.x_ + y_ * vec.y_;
	}

	//二维向量减法
	Vector2d Minus(const Vector2d vec) const
	{
		return Vector2d(x_ - vec.x_, y_ - vec.y_);
	}
};

//三角形类
class Triangle
{
public:
	Vector2d pointA_, pointB_, pointC_;

public:
	Triangle(Vector2d point1, Vector2d point2, Vector2d point3)
		:pointA_(point1), pointB_(point2), pointC_(point3)
	{
		//todo 判断三点是否共线
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