#include "utils.h"


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

int EdgeChainsToSeeds(EdgeChain &seeds, const EdgeChains edges) {

	if (!seeds.empty())
		seeds.clear();
	seeds.resize(edges.xCors.size());
	for (uint i = 0; i < edges.xCors.size(); i++) {

		seeds[i].x = edges.xCors[i];
		seeds[i].y = edges.yCors[i];
	}

	return 1;
}

int EdgeChainsToVecEdgeChains(VecEdgeChains &vedges, const EdgeChains edges) {

	Pixel dot;
	vedges.resize(edges.numOfEdges);
	for (uint i = 0; i < edges.numOfEdges; i++) {
		for (uint k = edges.sId[i]; k < (edges.sId[i + 1]); k++) {
			dot.x = edges.xCors[k];
			dot.y = edges.yCors[k];
			vedges[i].push_back(dot);
		}
	}

	return 1;
}

int VecEdgeChainsToEdgeChains(EdgeChains &edges, const VecEdgeChains vedges) {

	edges.sId.clear();
	edges.xCors.clear();
	edges.yCors.clear();
	edges.numOfEdges = vedges.size();
	unsigned int index = 0;
	edges.sId.push_back(index);
	for (int i = 0; i < vedges.size(); i++) {
		for (int j = 0; j < vedges[i].size(); j++) {
			edges.xCors.push_back(vedges[i][j].x);
			edges.yCors.push_back(vedges[i][j].y);
			index++;
		}
		edges.sId.push_back(index);
	}

	return 1;
}


DVector2D* EDToPath(EdgeChain edge) {

	size_t num = edge.size();
	DVector2D* segment = new DVector2D[num];
	for (int i = 0; i < num; i++) {
		segment[i].dX = edge[i].x;
		segment[i].dY = edge[i].y;
	}
	return segment;

}


EdgeChain PathToED(DVector2D* segment, unsigned int uPointsInCurrentPath) {

	unsigned int num = uPointsInCurrentPath;
	EdgeChain edge(num);
	for (uint i = 0; i < num; i++) {
		edge[i].x = unsigned int(segment[i].dX);
		edge[i].y = unsigned int(segment[i].dY);
	}

	return edge;
}

/*
功能：向edge中push p1,p1与p2之间,p2 这些Pixels
*/
void PixelsLinkToEdgeChain(EdgeChain& edge, Pixel p1, Pixel p2) {

	Pixel tpixel;
	EdgeChain tedge;
	bool yxflag = false;				// 如果dx/dy > 1, 则交换x,y后续需要交换回来;				p1(100,10)p2(10,5) --> p1(10,100)p2(5,10)
	bool xflag = false;					// 默认dx/dy<1,如果p1.x > p2.x,交换p1,p2,后续需要交换回来;	p1(10,100)p2(5,10) --> p1(5,10)p2(10,100)
	bool yflag = false;					// 默认dx/dy<1,如果p1.y > p2.y，ty -= 1,后续无需交换回来;   -----------------------------------------
	size_t currentiter = edge.size();	//记录插入散点的起始位置

	if (abs(int(p1.x - p2.x)) > abs(int(p1.y - p2.y))) {
		yxflag = true;
		Pixel tp1(p1.y, p1.x);
		p1 = tp1;
		Pixel tp2(p2.y, p2.x);
		p2 = tp2;
	}

	if (p1.x > p2.x) {
		Pixel tp = p1;
		p1 = p2;
		p2 = tp;
		xflag = true;
	}

	if (p1.y > p2.y) {

		yflag = true;
	}
	//int dx = fabs(p2.x - p1.x);
	//int dy = fabs(p2.y - p1.y);			//不知道为什么这个不行

	int dx = abs(int(p2.x - p1.x)) + 1;
	int dy = abs(int(p2.y - p1.y)) + 1;
	int radio = dy / dx;
	int rest = (dy % dx) / 2;
	int ty = p1.y;
	//前部分补点 + 1(or 0),
	if ((dy % dx) % 2) {
		for (int i = 1; i <= rest + 1; i++) {

			//Pixel tpixel(p1.x, ty);
			tpixel.x = p1.x;
			tpixel.y = ty;
			tedge.push_back(tpixel);
			ty += (1 + (-2)*yflag);
		}
	}
	else {
		for (int i = 1; i <= rest; i++) {

			//Pixel tpixel(p1.x, ty);
			tpixel.x = p1.x;
			tpixel.y = ty;
			tedge.push_back(tpixel);
			ty += (1 + (-2)*yflag);
		}
	}
	//中间的点,均分
	for (int j = p1.x; j <= p2.x; j++) {
		for (int i = 1; i <= radio; i++) {

			//Pixel tpixel(j, ty);
			tpixel.x = j;
			tpixel.y = ty;
			tedge.push_back(tpixel);
			ty += (1 + (-2)*yflag);
		}
	}
	//后部分补点
	for (int i = 0; i < rest; i++) {

		//Pixel tpixel(p2.x, ty);
		tpixel.x = p2.x;
		tpixel.y = ty;
		tedge.push_back(tpixel);
		ty += (1 + (-2)*yflag);
	}

	//写入edge
	if (xflag) {
		vector<Pixel>::reverse_iterator riter;
		for (riter = tedge.rbegin(); riter != tedge.rend(); riter++) {

			edge.push_back(*riter);
		}
	}
	else {
		vector<Pixel>::iterator riter;
		for (riter = tedge.begin(); riter != tedge.end(); riter++) {

			edge.push_back(*riter);
		}
	}
	if (yxflag) {
		for (size_t i = currentiter; i < edge.size(); i++) {

			int temp = 0;
			temp = edge[i].x;
			edge[i].x = edge[i].y;
			edge[i].y = temp;
		}
	}
}

/*
功能：由受约束的点扩展得到用来制作计算flow距离场的种子点，方式是将线段两点及两点之间的所有点判断为为种子点

	这段代码估计是白写了，因为直接用edges，已经是连续边缘了，就不用再连接线段两顶点了
*/
void EdgeChainsToJFASeeds(EdgeChain& Seeds, const EdgeChains& edges) {

	Seeds.clear();
	Pixel p1;
	Pixel p2;
	for (size_t i = 0; i < edges.numOfEdges; i++) {

		for (unsigned int j = edges.sId[i]; j + 1 < edges.sId[i + 1]; j++) {

			p1.x = edges.xCors[j];		p1.y = edges.yCors[j];
			p2.x = edges.xCors[j + 1];	p2.y = edges.yCors[j + 1];
			PixelsLinkToEdgeChain(Seeds, p1, p2);
			Seeds.pop_back();				//删除重复元素
		}
		unsigned int k = edges.sId[i];
		Pixel p(edges.xCors[k], edges.yCors[k]);
		Seeds.push_back(p);					//插入该这线段最后一个点
	}
}
/*----------------------------------------------------------------------*/
int DrawPlots(EdgeChains edges, unsigned int xsize,
	unsigned int ysize)
{

	Mat img(ysize, xsize, CV_8UC3, Scalar(0, 0, 0));
	vector<Point> P;
	size_t num = edges.xCors.size();
	P.resize(num);
	for (int i = 0; i < num; i++) {

		P[i].x = edges.xCors[i];
		P[i].y = edges.yCors[i];
		circle(img, P[i], 0, Scalar(255, 255, 255), -1);

	}
	imshow("Plots of Image", img);
	waitKey();
	return 0;
}

int DrawLines(EdgeChains edges, unsigned int xsize,
	unsigned int ysize)
{
	//cout << rand() * 255 << "\t" << rand() * 255 << "\t" << rand() * 255 << endl;
	Mat img(ysize, xsize, CV_8UC3, Scalar(0, 0, 0));
	vector<Point> P;
	//size_t num = edges.xCors.size();
	P.resize(edges.xCors.size());
	for (int i = 0; i < edges.xCors.size(); i++) {

		P[i].x = edges.xCors[i];
		P[i].y = edges.yCors[i];

	}
	for (uint j = 0; j < edges.numOfEdges; j++) {
		uint red = (rand() % (256 - 0)) + 0;
		uint green = (rand() % (256 - 0)) + 0;
		uint blue = (rand() % (256 - 0)) + 0;
		for (uint i = edges.sId[j]; i < edges.sId[j + 1] - 1; i++) {

			line(img, P[i], P[i + 1], Scalar(red, green, blue), 1, 1);

		}
	}
	imshow("Edges of Image", img);
	waitKey();
	return 0;
}

//int main() {
//
//	ImgSet Set;
//	Set.InitImg = GetInitImg("Mickey2.png");
//	Set.GrayImg = GetGrayImg(Set.InitImg);
//	Set.FilterImg = GetFilterImg(Set.GrayImg);
//	uint xsize = Getxsize(Set.FilterImg);
//	uint ysize = Getysize(Set.FilterImg);
//	uchar *data = GetImgData(Set.FilterImg);
//	EdgeChains edges = GetEdgeChains(data, xsize, ysize);
//
//	VecEdgeChains vedges;
//	EdgeChainsToVecEdgeChains(vedges, edges);
//	//VecEdgeChainsToEdgeChains(edges, vedges);
//
//	EdgeChains temp = edges;
//
//
//	for (int i = 0; i < vedges.size(); i++) {
//		DVector2D *segment = EDToPath(vedges[i]);
//		unsigned int uPointsInCurrentPath = vedges[i].size();
//		unsigned int puPointsInResultPath = vedges[i].size();
//		cout << "before" << uPointsInCurrentPath << endl;
//		compactPath(segment, uPointsInCurrentPath, segment, &uPointsInCurrentPath, 0.00000000001, perpendicularDistanceDeviationMetric);
//		cout << "after" << uPointsInCurrentPath << endl;
//		vedges[i] = PathToED(segment, uPointsInCurrentPath);
//		delete segment;
//	}
//	VecEdgeChainsToEdgeChains(edges, vedges);
//
//	//for (int i = 0; i < edges.xCors.size(); i++) {
//	//	if (edges.sId[i] == temp.sId[i])
//	//		;
//	//	else
//	//		cout << i << endl;
//	//}
//
//	DrawLines(edges, xsize, ysize);
//
//	cout << xsize << "\t" << ysize << endl;
//
//}


//int main() {
//
//	EdgeChains edges;
//
//	edges.xCors = { 1,2,3,4,5,6,7,8 };
//	edges.yCors = { 1,2,3,4,5,6,7,20 };
//	edges.numOfEdges = 8;
//
//	DVector2D *segment = NULL;
//	segment = EDToPath(edges, segment);
//	unsigned int uPointsInCurrentPath = edges.numOfEdges;
//	for (int i = 0; i < edges.numOfEdges; i++)
//		cout << segment[i].dX << " " << segment[i].dY << endl;
//
//	cout << endl;
//
//
//	compactPath(segment, uPointsInCurrentPath, segment, &uPointsInCurrentPath, 5.0, perpendicularDistanceDeviationMetric);
//
//	for (int i = 0; i < uPointsInCurrentPath; i++)
//		cout << segment[i].dX << " " << segment[i].dY << endl;
//
//}