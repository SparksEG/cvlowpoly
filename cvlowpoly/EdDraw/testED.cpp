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
	VecEdgeChainsToEdgeChains(edges, vedges);   //重新转为edges，目的是减去尾部的（0,0）
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
	//EdgeChainsToJFASeeds(JFASeeds, edges);					 //是利用所有的edges点连成的线，慢一点，感觉正确一点
	EdgeChain JFASeeds;
	JFASeeds.resize(edges.xCors.size());
	for (uint i = 0; i < edges.xCors.size(); i++) {

		JFASeeds[i].x = edges.xCors[i];
		JFASeeds[i].y = edges.yCors[i];
	}			//这一部分可以写成一个函数


	//EdgeChainsToJFASeeds(Seeds, edgesPathCompact);	         //是利用简化后的edges点连成的线，快一点
	//EdgeChain JFASeeds;
	//JFASeeds.resize(edgesPathCompact.xCors.size());
	//for (uint i = 0; i < edgesPathCompact.xCors.size(); i++) {

	//	JFASeeds[i].x = edgesPathCompact.xCors[i];
	//	JFASeeds[i].y = edgesPathCompact.yCors[i];
	//}			//这一部分可以写成一个函数

	
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


	//合并 采样点 与 关键点
	//绘制三角剖分模型
	Size size = Set.InitImg.size();
	Rect rect(0, 0, size.width, size.height);
	Subdiv2D subdiv(rect);

	 
	for (int i = 0; i <optSeeds.xCors.size(); i++) {
		cv::Point2f p;
		p.x = optSeeds.xCors[i];
		p.y = optSeeds.yCors[i];
		subdiv.insert(p);
	}

	//三角剖分简化的关键点
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
	////三角剖分可视化全部的关键点
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
	subdiv.getTriangleList(triangleList);	//获得三角剖分的三角形
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
第一次采样较密，泊松点采样即可少；

第一次采样较稀，泊松点采样需要多，迭代次数？
多也不需要太多，看迭代次数？



加！号次数越多越不好    （对于edges不足）
加！号次数越多，泊松点越在edges附近
泊松采样是为了弥补edges点不足


初步判定不加！号好（对于edges不足）
初步判定泊松点多好（对于edges不足）
加！号次数越多，泊松点越在edges附近

为什么不加！号：
因为没有泊松点采样时，只有边界上的点，此时会有大量尖锐的三角形出现；
而加入泊松点采样，更加均匀的分布在非边界点的地方，降低了尖锐三角形效果


*/