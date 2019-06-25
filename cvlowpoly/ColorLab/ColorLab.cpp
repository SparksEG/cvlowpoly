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


// ������������L������40%-60%��LAB��Ӧȡƽ��;
void ColorPointInTriangle(vector<vector<bool>> IsPointScanedMap, const Mat& LabMap, vector<vector<Vec3f>>& ColorMap, const Vec6f& triangle) {

	Triangle t(Vector2d(triangle[1], triangle[0]), Vector2d(triangle[3], triangle[2]), Vector2d(triangle[5], triangle[4]));
	vector<Vec3b> VecLAB;

	int xmin = min(triangle[0], triangle[2], triangle[4]);
	int xmax = max(triangle[0], triangle[2], triangle[4]);
	int ymin = min(triangle[1], triangle[3], triangle[5]);
	int ymax = max(triangle[1], triangle[3], triangle[5]);

	if (xmin < 0 || xmin >= LabMap.cols || xmax < 0 || xmax >= LabMap.cols || ymin < 0 || ymin >= LabMap.rows || ymax < 0 || ymax >= LabMap.rows)
		return;

	// ��Ϊ ͬһL���� ���ܶ�Ӧ����㼴���RGB����Ҫ��multimap
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

	//��ȡ40%-60%��L������ƽ����ɫ
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
	

	//����ƽ����ɫ
	for (int j = ymin; j <= ymax; j++) {
		for (int i = xmin; i <= xmax; i++) {

			if (IsPointScanedMap[j][i]) {
				ColorMap[j][i] = LAB;
			}
		}
	}
}

/*
// �����ģ�ȫ��ȥһ�������ֵ,��ȷ��
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

//// ����һ��L����ȡ40%-60%��RGBȫ��ȡƽ��
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
//	//��ȡ40%-60%��L����
//	float AverL = 0;
//	int colorNum = VecLcomponent.size();
//	int begin = colorNum * 0.4;
//	int end = colorNum * 0.6;
//	for (int i = begin; i < end; i++) {
//		AverL += VecLcomponent[i];
//	}
//	AverL = AverL / (end - begin + 1);
//
//	//��ȡƽ����ɫ
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
//	//����ƽ����ɫ
//	for (int j = ymin; j <= ymax; j++) {
//		for (int i = xmin; i <= xmax; i++) {
//
//			if (IsPointScanedMap[j][i]) {
//				ColorMap[j][i] = AverLAB;
//				ColorMap[j][i][0] = AverL;  //�ǲ�����ô���İ�
//			}
//		}
//	}
//}

/*
// ������������L������40%-60%��RGB��Ӧȡƽ��; ��unordered_multimap������Ӧ��ϵ
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

// ��Ϊ ͬһL���� ���ܶ�Ӧ����㼴���RGB����Ҫ��multimap
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

//��ȡ40%-60%��L������ƽ����ɫ
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

//����ƽ����ɫ
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



