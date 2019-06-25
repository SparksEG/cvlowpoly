#include <assert.h>
#include "JFA.h"
#include <unordered_map>
#define random(x) (rand()%x)
/*
���ܣ������ڶ����Ż�����������λ�ò���ĵ㣬��������PathCompact�򻯹��Ķ��㣨����Լ���Ķ��㣩��ͼƬ���ĸ�����

�����constrainedFlagMap ��¼���Ƿ�������µ�logic��0������/1������ͼ
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
���ܣ�һ���Ż����������λ��

cvtSeeds�����еĶ��㣬��ʽ�Ѳ������߶ε���Ϣ�����Կ�����һЩɢ��
optTimes�������Ż��Ĵ���
distMap�� flow�����ף���¼��ͼ���е�ÿһ���㵽���ж���ľ�����Ϣ��������ŷʽ���룩
constrainedFlagMap���ڶ����Ż����������У�λ�ò���Ķ���

�����cvtSeeds λ�õõ������Ż������Ķ���
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
//				unsigned int idx = pos->second;   //ָ���cvtSeeds���ӵ���Χ�����ٽ���
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


	//��¼��edges�������λ��idx����������edges[i]��i
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
				int i = idxOfcvtSeedsId[VoiMap[idx]];  //�ɹؼ���VoiMap[idx](�������ص����������ӵ��Idx)��������Ӧ��vector�����i
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
���ܣ�ʵ�ʵ��ú����������Ż�

constrainedPoints����Լ���ĵ㣬�Ǿ���PathCompact�򻯹��Ķ���
cvtSeeds�����еĶ��㣬��ʽ�Ѳ������߶ε���Ϣ�����Կ�����һЩɢ��
optTimes�������Ż��Ĵ���
distMap�� flow�����ף���¼��ͼ���е�ÿһ���㵽���ж���ľ�����Ϣ��������ŷʽ���룩
constrainedFlagMap���ڶ����Ż����������У�λ�ò���Ķ���

���������optSeeds λ�õõ������Ż������Ķ���
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


