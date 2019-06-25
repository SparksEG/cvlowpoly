#include "ColorLab.h"


bool IsSmallOnVec3b(Vec3b v1, Vec3b v2) {

	return v1[0] < v2[0];
}


float min(float a, float b, float c) {

	int min = a < b ? a : b;
	return min < c ? min : c;
}

float max(float a, float b, float c) {

	int max = a > b ? a : b;
	return max > c ? max : c;
}

bool IsPointInTriangle(const Vector2d p, const Triangle t) {

	Vector2d PA = t.pointA_.Minus(p);
	Vector2d PB = t.pointB_.Minus(p);
	Vector2d PC = t.pointC_.Minus(p);
	double t1 = PA.CrossProduct(PB);
	double t2 = PB.CrossProduct(PC);
	double t3 = PC.CrossProduct(PA);

	return (t1 * t2 >= 0) && (t1 * t3 >= 0) && (t2 * t3 >= 0);
}


// 方法三：依照L分量的40%-60%，LAB对应取平均;
void ColorPointInTriangle(vector<vector<bool>> IsPointScanedMap, const Mat& LabMap, vector<vector<Vec3f>>& ColorMap, const Vec6f& triangle) {

	Triangle t(Vector2d(triangle[1], triangle[0]), Vector2d(triangle[3], triangle[2]), Vector2d(triangle[5], triangle[4]));
	vector<Vec3b> VecLAB;

	int xmin = min(triangle[0], triangle[2], triangle[4]);
	int xmax = max(triangle[0], triangle[2], triangle[4]);
	int ymin = min(triangle[1], triangle[3], triangle[5]);
	int ymax = max(triangle[1], triangle[3], triangle[5]);

	if (xmin < 0 || xmin >= LabMap.cols || xmax < 0 || xmax >= LabMap.cols || ymin < 0 || ymin >= LabMap.rows || ymax < 0 || ymax >= LabMap.rows)
		return;

	// 因为 同一L分量 可能对应多个点即多个RGB，需要用multimap
	for (int j = ymin; j <= ymax; j++) {
		for (int i = xmin; i <= xmax; i++) {

			if (!IsPointScanedMap[j][i]) {
				Vector2d p(j, i);
				if (IsPointInTriangle(p, t)) {;
					VecLAB.push_back(LabMap.at<Vec3b>(j, i));
					IsPointScanedMap[j][i] = true;
				}
			}
		}
	}
	if (VecLAB.empty())
		return;
	sort(VecLAB.begin(), VecLAB.end(), IsSmallOnVec3b);

	//求取40%-60%的L分量的平均颜色
	Vec3f LAB(0, 0, 0);
	int count = 0;
	int colorNum = VecLAB.size();
	int begin = colorNum * 0.4;
	int end = colorNum * 0.6;
	//for (int i = 0; i < VecLAB.size(); ++i) {
	for (int i = begin; i <= end; ++i) {
		Vec3f tpLAB = VecLAB[i];
		//LAB = LAB + tpLAB;
		LAB = LAB + static_cast<Vec3f>(VecLAB[i]);
	}
	//LAB = LAB / int(VecLAB.size());
	LAB = LAB / (end - begin + 1);
	

	//复制平均颜色
	for (int j = ymin; j <= ymax; j++) {
		for (int i = xmin; i <= xmax; i++) {

			if (IsPointScanedMap[j][i]) {
				ColorMap[j][i] = LAB;
			}
		}
	}
}

/*
// 方法四：全部去一样的随机值,正确的
void ColorPointInTriangle(vector<vector<bool>>& IsPointScanedMap, const Mat& LabMap, vector<vector<Vec3f>>& ColorMap, const Vec6f& triangle) {

	Triangle t(Vector2d(triangle[1], triangle[0]), Vector2d(triangle[3], triangle[2]), Vector2d(triangle[5], triangle[4]));
	vector<Vec3b> VecLAB;

	int xmin = min(triangle[0], triangle[2], triangle[4]);
	int xmax = max(triangle[0], triangle[2], triangle[4]);
	int ymin = min(triangle[1], triangle[3], triangle[5]);
	int ymax = max(triangle[1], triangle[3], triangle[5]);

	if (xmin < 0 || xmin >= LabMap.cols || xmax < 0 || xmax >= LabMap.cols || ymin < 0 || ymin >= LabMap.rows || ymax < 0 || ymax >= LabMap.rows)
		return;
	uint red = (rand() % (256 - 0)) + 0;
	uint green = (rand() % (256 - 0)) + 0;
	uint blue = (rand() % (256 - 0)) + 0;
	Vec3b color(blue, green, red);
	for (int j = ymin; j <= ymax; j++) {
		for (int i = xmin; i <= xmax; i++) {
			Vector2d p(j, i);
			if (!IsPointScanedMap[j][i] && IsPointInTriangle(p, t)) {
				ColorMap[j][i] = color;
				IsPointScanedMap[j][i] = true;
			}
		}
	}
}
*/

//// 方法一：L分量取40%-60%，RGB全部取平均
//void ColorPointInTriangle(vector<vector<bool>>& IsPointScanedMap, const vector<vector<float>> L_component, vector<vector<Scalar>>& ColorMap, const Vec6f triangle) {
//
//	Triangle t(Vector2d(triangle[0], triangle[1]), Vector2d(triangle[2], triangle[3]), Vector2d(triangle[4], triangle[5]));
//	vector<float> VecLcomponent;
//
//	int xmin = min(triangle[0], triangle[2], triangle[4]);
//	int xmax = max(triangle[0], triangle[2], triangle[4]);
//	int ymin = min(triangle[1], triangle[3], triangle[5]);
//	int ymax = max(triangle[1], triangle[3], triangle[5]);
//
//	for (int j = ymin; j <= ymax; j++) {
//		for (int i = xmin; i <= xmax; i++) {
//
//			if (!IsPointScanedMap[j][i]) {
//				Vector2d p(j, i);
//				if (IsPointInTriangle(p, t)) {
//					VecLcomponent.push_back(L_component[j][i]);
//					IsPointScanedMap[j][i] = true;
//				}
//			}
//		}
//	}
//	sort(VecLcomponent.begin(), VecLcomponent.end());
//
//	//求取40%-60%的L分量
//	float AverL = 0;
//	int colorNum = VecLcomponent.size();
//	int begin = colorNum * 0.4;
//	int end = colorNum * 0.6;
//	for (int i = begin; i < end; i++) {
//		AverL += VecLcomponent[i];
//	}
//	AverL = AverL / (end - begin + 1);
//
//	//求取平均颜色
//	Scalar AverLAB(0, 0, 0);
//	for (int j = ymin; j <= ymax; j++) {
//		for (int i = xmin; i <= xmax; i++) {
//
//			if (IsPointScanedMap[j][i]) {
//				AverLAB += ColorMap[j][i];
//			}
//		}
//	}
//	AverLAB = AverLAB / (end - begin + 1);
//
//	//复制平均颜色
//	for (int j = ymin; j <= ymax; j++) {
//		for (int i = xmin; i <= xmax; i++) {
//
//			if (IsPointScanedMap[j][i]) {
//				ColorMap[j][i] = AverLAB;
//				ColorMap[j][i][0] = AverL;  //是不是这么来的啊
//			}
//		}
//	}
//}

/*
// 方法二：依照L分量的40%-60%，RGB对应取平均; 用unordered_multimap建立对应关系
void ColorPointInTriangle(vector<vector<bool>>& IsPointScanedMap, const Mat& LabMap, vector<vector<Vec3f>>& ColorMap, const Vec6f& triangle) {

Triangle t(Vector2d(triangle[1], triangle[0]), Vector2d(triangle[3], triangle[2]), Vector2d(triangle[5], triangle[4]));
vector<int> VecLcomponent;
unordered_multimap<int, Vec3b> L2Color;
int xmin = min(triangle[0], triangle[2], triangle[4]);
int xmax = max(triangle[0], triangle[2], triangle[4]);
int ymin = min(triangle[1], triangle[3], triangle[5]);
int ymax = max(triangle[1], triangle[3], triangle[5]);

if (xmin < 0 || xmin >= LabMap.cols || xmax < 0 || xmax >= LabMap.cols || ymin < 0 || ymin >= LabMap.rows || ymax < 0 || ymax >= LabMap.rows)
return;

// 因为 同一L分量 可能对应多个点即多个RGB，需要用multimap
for (int j = ymin; j <= ymax; j++) {
for (int i = xmin; i <= xmax; i++) {

if (!IsPointScanedMap[j][i]) {
Vector2d p(j, i);
if (IsPointInTriangle(p, t)) {
int L_component = LabMap.at<Vec3b>(j,i)[0];
L2Color.insert(make_pair(L_component, ColorMap[j][i]));
VecLcomponent.push_back(L_component);
IsPointScanedMap[j][i] = true;
}
}
}
}
if (VecLcomponent.empty())
return;
sort(VecLcomponent.begin(), VecLcomponent.end());

//求取40%-60%的L分量的平均颜色
Vec3f RGB(0, 0, 0);
int count = 0;
int colorNum = VecLcomponent.size();
int begin = colorNum * 0.4;
int end = colorNum * 0.6;
for (int i = begin; i < end; i++) {
unordered_multimap<int, Vec3b>::iterator iter;
for (iter = L2Color.lower_bound(VecLcomponent[i]); iter != L2Color.upper_bound(VecLcomponent[i]); ++iter) {
RGB += iter->second;
++count;
}

}
Vec3b AverRGB = RGB / count;

//复制平均颜色
for (int j = ymin; j <= ymax; j++) {
for (int i = xmin; i <= xmax; i++) {

if (IsPointScanedMap[j][i]) {
ColorMap[j][i] = AverRGB;
}
}
}
}
*/

void ReadImgColorOnOpencv(const Mat img, vector<vector<Vec3f>>& ColorMap) {

	int ysize = img.rows;
	int xsize = img.cols;

	for (int j = 0; j < ysize; j++) {
		for (int i = 0; i < xsize; i++) {
			
			ColorMap[j][i] = img.at<Vec3b>(j, i) ;
		}
	}
}

void WriteImgColorOnOpencv(const vector<vector<Vec3f>> ColorMap, Mat& img) {

	int ysize = img.rows;
	int xsize = img.cols;

	for (int j = 0; j < ysize; j++) {
		for (int i = 0; i < xsize; i++) {

			img.at<Vec3b>(j, i) = ColorMap[j][i];
		}
	}
}



