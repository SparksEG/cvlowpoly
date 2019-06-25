#include <assert.h>
#include "JFA.h"
#include <unordered_map>
#define random(x) (rand()%x)
/*
功能：制作在顶点优化迭代过程中位置不变的点，包括经过PathCompact简化过的顶点（即受约束的顶点）和图片的四个顶点

输出：constrainedFlagMap 记录的是否迭代更新的logic（0不迭代/1迭代）图
*/
//void imposeConstrainedPoints(bool* constrainedFlagMap, const EdgeChain constrainedPoints, unsigned int xsize, unsigned int ysize) {
//	
//	for (unsigned int y = 0; y < ysize; ++y) {
//
//		for (unsigned int x = 0; x < xsize; ++x) {
//
//			unsigned int idx = (y * xsize) + x;
//			constrainedFlagMap[idx] = false;
//		}
//	}
//
//	for (size_t i = 0; i < constrainedPoints.size(); i++) {
//
//		int idx = (constrainedPoints[i].y * xsize) + constrainedPoints[i].x;
//		constrainedFlagMap[idx] = true;
//	}
//	constrainedFlagMap[0] = true;
//	constrainedFlagMap[xsize - 1] = true;
//	constrainedFlagMap[(ysize - 1)*xsize] = true;
//	constrainedFlagMap[ysize*xsize - 1] = true;
//}

void imposeConstrainedPoints(bool* constrainedFlagMap, const EdgeChains constrainedPoints, unsigned int xsize, unsigned int ysize) {

	for (unsigned int y = 0; y < ysize; ++y) {

		for (unsigned int x = 0; x < xsize; ++x) {

			unsigned int idx = (y * xsize) + x;
			constrainedFlagMap[idx] = false;
		}
	}

	for (size_t i = 0; i < constrainedPoints.xCors.size(); i++) {

		int idx = (constrainedPoints.yCors[i] * xsize) + constrainedPoints.xCors[i];
		constrainedFlagMap[idx] = true;
	}
	constrainedFlagMap[0] = true;
	constrainedFlagMap[xsize - 1] = true;
	constrainedFlagMap[(ysize - 1)*xsize] = true;
	constrainedFlagMap[ysize*xsize - 1] = true;
}

/*
功能：一次优化迭代顶点的位置

cvtSeeds：所有的顶点，格式已不保留线段的信息，可以看做是一些散点
optTimes：顶点优化的次数
distMap： flow距离谱，记录了图像中的每一个点到所有顶点的距离信息（不简单是欧式距离）
constrainedFlagMap：在顶点优化迭代过程中，位置不变的顶点

输出：cvtSeeds 位置得到更新优化迭代的顶点
*/
//void cvtOptimization(EdgeChain& cvtSeeds,
//	unsigned int xsize, unsigned int ysize, int optTimes, const int* distMap, const bool* constrainedFlagMap ) {
//
//	VoiPixelMap VoiMap;
//	ExecuteJumpFloodingVoi(VoiMap, cvtSeeds, xsize, ysize);
//	for (size_t i = 0; i < cvtSeeds.size(); i++) {
//
//		int seedIdx = cvtSeeds[i].y*xsize + cvtSeeds[i].x;
//		if (constrainedFlagMap[seedIdx] == false) {
//
//			Pixel p = cvtSeeds[i];
//			float iterX = 0;
//			float iterY = 0;
//			float iterW = 0;
//
//			for (VoiPixelMap::iterator pos = VoiMap.lower_bound(seedIdx); pos != VoiMap.upper_bound(seedIdx); ++pos) {
//
//				unsigned int idx = pos->second;   //指向该cvtSeeds种子点周围的最临近点
//				unsigned int y = idx / xsize;
//				unsigned int x = idx - xsize * y;
//				float w = distMap[idx];
//
//				iterX += x * w;
//				iterY += y * w;
//				iterW += w;
//			}
//
//			cvtSeeds[i].x = int(iterX / iterW);
//			cvtSeeds[i].y = int(iterY / iterW);
//		}
//	}
//}

void cvtOptimization(Seeds& cvtSeeds, unsigned int xsize, unsigned int ysize, 
	int optTimes, float* flowDistMap/*, const bool* constrainedFlagMap*/) {


	//记录由edges点的像素位置idx可以索引到edges[i]的i
	unordered_map<int, int> idxOfcvtSeedsId;
	for (size_t i = 0; i < cvtSeeds.size(); i++) {

		int idx = cvtSeeds[i].y * xsize + cvtSeeds[i].x;
		idxOfcvtSeedsId[idx] = i;
	}

	int* VoiMap = (int*)malloc(sizeof(int) * xsize * ysize);
	ExecuteJumpFloodingVoi(VoiMap, cvtSeeds, xsize, ysize);
	
	vector<float> tx(cvtSeeds.size());
	vector<float> ty(cvtSeeds.size());
	vector<float> tw(cvtSeeds.size());
	for (unsigned int y = 0; y < ysize; ++y) {

		for (unsigned int x = 0; x < xsize; ++x) {

			int idx = (y * xsize) + x;
			//if (constrainedFlagMap[idx] == false)
			{

				float w = flowDistMap[idx];
				int i = idxOfcvtSeedsId[VoiMap[idx]];  //由关键字VoiMap[idx](即该像素点所属的种子点的Idx)索引到对应的vector的序号i
				tx[i] += x * w;
				ty[i] += y * w;
				tw[i] += w;

				//tx[i] += x;
				//ty[i] += y;
				//tw[i] ++;

				//unordered_map<int, int>::const_iterator iter = idxOfcvtSeedsId.find(seedsIdx);
				//if (iter != idxOfcvtSeedsId.end()) {

				//	int i = iter->second;//int i = idxOfcvtSeedsId[VoiMap[idx]];
				//	tx[i] += x * w;
				//	ty[i] += y * w;
				//	tw[i] += w;
				//}
			}
		}
	}
	srand((int)time(0));
	for (size_t i = 0; i < cvtSeeds.size(); i++) {

		//int idx = cvtSeeds[i].y * xsize + cvtSeeds[i].x;
		//if (constrainedFlagMap[idx] == false) 
		{
			if (tw[i] == 0)
				cout << "x:" << cvtSeeds[i].x << "y:" << cvtSeeds[i].y << endl;
			cvtSeeds[i].x = int(tx[i] / tw[i]);
			cvtSeeds[i].y = int(ty[i] / tw[i]);

			if (cvtSeeds[i].x < 0 || cvtSeeds[i].x >= xsize) {
				cout << "x:" << cvtSeeds[i].x << "y:" << cvtSeeds[i].y << endl;
				cvtSeeds[i].x = random(xsize);
			}
			if (cvtSeeds[i].y < 0 || cvtSeeds[i].y >= ysize) {
				cout << "x:" << cvtSeeds[i].x << "y:" << cvtSeeds[i].y << endl;
				cvtSeeds[i].y = random(ysize);
			}

		}
	}
}


/*
功能：实际调用函数，顶点优化

constrainedPoints：受约束的点，是经过PathCompact简化过的顶点
cvtSeeds：所有的顶点，格式已不保留线段的信息，可以看做是一些散点
optTimes：顶点优化的次数
distMap： flow距离谱，记录了图像中的每一个点到所有顶点的距离信息（不简单是欧式距离）
constrainedFlagMap：在顶点优化迭代过程中，位置不变的顶点

输出：返回optSeeds 位置得到更新优化迭代的顶点
*/
Seeds VertexOptimization(/*const EdgeChains constrainedPoints,*/ const Seeds cvtSeeds,
	unsigned int xsize, unsigned int ysize, int optTimes, float* flowDistMap) {

	Seeds optSeeds = cvtSeeds;
	//bool* constrainedFlagMap = (bool *)malloc(sizeof(bool)*xsize*ysize);
	//imposeConstrainedPoints(constrainedFlagMap, constrainedPoints, xsize, ysize);
	for (int i = 0; i < optTimes; i++) {
		clock_t begin = clock();

		cvtOptimization(optSeeds, xsize, ysize, optTimes, flowDistMap/*, constrainedFlagMap*/);

		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		cout << elapsed_secs << endl;
	}

	return optSeeds;
};


