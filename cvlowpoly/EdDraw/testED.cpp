#include "edlines.h"
#include "utils.h"
#include <ctime>

#include "triangle.h"
#include "RandomSeeds.h"
#include "ColorLab.h"
//// The following struct represents a unsigned int precision 2D point.
////struct Pixel {
////	unsigned int x;//X coordinate
////	unsigned int y;//Y coordinate
////};
//typedef vector<Pixel> EdgeChain;
//typedef vector<vector<Pixel> > VecEdgeChains;
//
//int EdgeChainsToVecEdgeChains(VecEdgeChains &vedges, const EdgeChains edges) {
//
//	Pixel dot;
//	vedges.resize(edges.numOfEdges);
//	unsigned int index = 0;
//	for (int i = 0; i < edges.numOfEdges; i++) {
//		for (int k = 0; k < (edges.sId[i+1] - edges.sId[i]) && index < edges.xCors.size(); k++) {
//			dot.x = edges.xCors[index];
//			dot.y = edges.yCors[index];
//			vedges[i].push_back(dot);
//			index++;
//		}
//	}
//	return 1;
//
//}

#define yita 0.02;
//float MINLength;

int main() {

	ImgSet Set;
	Set.InitImg = GetInitImg("img1.png");
	Set.GrayImg = GetGrayImg(Set.InitImg);
	Set.FilterImg = GetFilterImg(Set.GrayImg);
	uint xsize = Getxsize(Set.FilterImg);
	uint ysize = Getysize(Set.FilterImg);
	uchar *data = GetImgData(Set.FilterImg);
	EdgeChains edges = GetEdgeChains(data, xsize, ysize);
	DeviationMetric DistanceDeviationMetric = &perpendicularDistance;
	//DeviationMetric DistanceToSegmentDeviationMetric = &shortestDistanceToSegment;

	float MINLength = (xsize + ysize) * yita;
	VecEdgeChains vedges;
	EdgeChainsToVecEdgeChains(vedges, edges);
	VecEdgeChainsToEdgeChains(edges, vedges);   //����תΪedges��Ŀ���Ǽ�ȥβ���ģ�0,0��
	for (int i = 0; i < vedges.size(); i++) {
		DVector2D *segment = EDToPath(vedges[i]);
		unsigned int upointsincurrentpath = vedges[i].size();
		compactPath(segment, upointsincurrentpath, segment, 
			&upointsincurrentpath, 1, DistanceDeviationMetric);
		vedges[i] = PathToED(segment, upointsincurrentpath);
		delete segment;
	}


	EdgeChains edgesPathCompact;
	VecEdgeChainsToEdgeChains(edgesPathCompact, vedges);
	


	
	//clock_t begin = clock();
	//clock_t end = clock();
	//double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//cout << elapsed_secs << endl;
	//EdgeChain JFASeeds;
	//EdgeChainsToJFASeeds(JFASeeds, edges);					 //���������е�edges�����ɵ��ߣ���һ�㣬�о���ȷһ��
	EdgeChain JFASeeds;
	JFASeeds.resize(edges.xCors.size());
	for (uint i = 0; i < edges.xCors.size(); i++) {

		JFASeeds[i].x = edges.xCors[i];
		JFASeeds[i].y = edges.yCors[i];
	}			//��һ���ֿ���д��һ������


	//EdgeChainsToJFASeeds(Seeds, edgesPathCompact);	         //�����ü򻯺��edges�����ɵ��ߣ���һ��
	//EdgeChain JFASeeds;
	//JFASeeds.resize(edgesPathCompact.xCors.size());
	//for (uint i = 0; i < edgesPathCompact.xCors.size(); i++) {

	//	JFASeeds[i].x = edgesPathCompact.xCors[i];
	//	JFASeeds[i].y = edgesPathCompact.yCors[i];
	//}			//��һ���ֿ���д��һ������

	
	float* flowDisMap = (float*)malloc(sizeof(float) * xsize * ysize);
	ExecuteJumpFloodingDis(flowDisMap, JFASeeds, xsize, ysize);

	int optTimes = 15;
	int NumPoints = 250;
	Seeds RandomSeeds = creatRandomSeeds(xsize, ysize, NumPoints);
	Seeds optEdges = VertexOptimization(RandomSeeds, xsize, ysize, optTimes, flowDisMap);

	DrawLines(edges, xsize, ysize);
	DrawLines(edgesPathCompact, xsize, ysize);

	//DrawPlots(edges, xsize, ysize);

	EdgeChains optSeeds;
	optSeeds.xCors.resize(optEdges.size());
	optSeeds.yCors.resize(optEdges.size());
	for (size_t i = 0; i < optEdges.size(); i++) {
		optSeeds.xCors[i] = optEdges[i].x;
		optSeeds.yCors[i] = optEdges[i].y;
	}

	DrawPlots(optSeeds, xsize, ysize);
	cout << xsize << "\t" << ysize << endl;


	//�ϲ� ������ �� �ؼ���
	//���������ʷ�ģ��
	Size size = Set.InitImg.size();
	Rect rect(0, 0, size.width, size.height);
	Subdiv2D subdiv(rect);

	 
	for (int i = 0; i <optSeeds.xCors.size(); i++) {
		cv::Point2f p;
		p.x = optSeeds.xCors[i];
		p.y = optSeeds.yCors[i];
		subdiv.insert(p);
	}

	//�����ʷּ򻯵Ĺؼ���
	for (int i = 0; i <edgesPathCompact.xCors.size(); i++) {
		cv::Point2f p;
		p.x = edgesPathCompact.xCors[i];
		p.y = edgesPathCompact.yCors[i];
		subdiv.insert(p);
	}
	subdiv.insert(cv::Point2f (0, 0));
	subdiv.insert(cv::Point2f (0, ysize-1));
	subdiv.insert(cv::Point2f (xsize - 1, 0));
	subdiv.insert(cv::Point2f (xsize - 1, ysize - 1));
	////�����ʷֿ��ӻ�ȫ���Ĺؼ���
	//for (int i = 0; i <edges.xCors.size(); i++) {
	//	cv::Point2f p;
	//	p.x = edges.xCors[i];
	//	p.y = edges.yCors[i];
	//	subdiv.insert(p);
	//}



	Mat img(ysize, xsize, CV_8UC3, Vec3b(0, 0, 0));
	Vec3b delaunay_color(255, 0, 255), points_color(0, 0, 255);

	draw_delaunay(img, subdiv, delaunay_color);
	imshow("Delaunay", img);
	waitKey();
	
	vector<vector<Vec3f>> ColorMap(ysize);
	fill(ColorMap.begin(), ColorMap.end(), vector<Vec3f>(xsize, (0,0,0)));
	vector<vector<bool> > IsPointScanedMap(ysize);
	fill(IsPointScanedMap.begin(), IsPointScanedMap.end(), vector<bool>(xsize, false));

	

	//ColorPointInTriangle(IsPointScanedMap, L_component, ColorMap, const Vec6f triangle)
	Mat LabMap(ysize, xsize, CV_8UC3, Vec3b(0, 0, 0));
	Mat RGBMap(ysize, xsize, CV_8UC3, Vec3b(0, 0, 0));
	cvtColor(Set.InitImg, LabMap, COLOR_RGB2Lab);
	ReadImgColorOnOpencv(LabMap, ColorMap);

	std::vector<Vec6f> triangleList;
	subdiv.getTriangleList(triangleList);	//��������ʷֵ�������
	for (size_t i = 0; i < triangleList.size(); ++i) {
		ColorPointInTriangle(IsPointScanedMap, LabMap, ColorMap, triangleList[i]);
	}
	WriteImgColorOnOpencv(ColorMap, LabMap);
	cvtColor(LabMap, RGBMap, COLOR_Lab2RGB);
	imshow("Res", RGBMap);
	waitKey();
	int a = 1;
}

/*
��һ�β������ܣ����ɵ���������٣�

��һ�β�����ϡ�����ɵ������Ҫ�࣬����������
��Ҳ����Ҫ̫�࣬������������



�ӣ��Ŵ���Խ��Խ����    ������edges���㣩
�ӣ��Ŵ���Խ�࣬���ɵ�Խ��edges����
���ɲ�����Ϊ���ֲ�edges�㲻��


�����ж����ӣ��źã�����edges���㣩
�����ж����ɵ��ã�����edges���㣩
�ӣ��Ŵ���Խ�࣬���ɵ�Խ��edges����

Ϊʲô���ӣ��ţ�
��Ϊû�в��ɵ����ʱ��ֻ�б߽��ϵĵ㣬��ʱ���д�������������γ��֣�
�����벴�ɵ���������Ӿ��ȵķֲ��ڷǱ߽��ĵط��������˼���������Ч��


*/