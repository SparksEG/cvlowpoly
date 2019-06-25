#include <map> 
#include <assert.h>
#include "JFA.h"
using namespace std;

/*=================================================================================================
DEFINES
=================================================================================================*/

// Initial window dimensions


#define  MINLength 10 //����  �������е�Li
/*=================================================================================================
FUNCTIONS
=================================================================================================*/

// If the buffers exist, delete them
void ClearBuffers(Pixel* BufferA, Pixel* BufferB) {

	if (BufferA != NULL) {
		free(BufferA);
		BufferA = NULL;
	}

	if (BufferB != NULL) {
		free(BufferB);
		BufferB = NULL;
	}

}


// Jump Flooding Algorithm
void ExecuteJumpFloodingDis(int* disMap, EdgeChain Seeds, unsigned int xsize, unsigned int ysize) {

	// Buffers
	Pixel* BufferA = NULL;
	Pixel* BufferB = NULL;

	// Buffer dimensions
	int BufferWidth = xsize;
	int BufferHeight = ysize;

	// Which buffer are we reading from?
	bool ReadingBufferA = true;


	// No seeds will just give us a black screen :P
	if (Seeds.size() < 1) {
		printf("Please create at least 1 seed.\n");
	}

	printf("Executing the Jump Flooding algorithm...\n");

	// Clear the buffers before we start
	//ClearBuffers();

	// Allocate memory for the two buffers
	int* distmap = NULL;
	distmap = (int*)malloc(sizeof(int) *BufferWidth * BufferHeight);
	BufferA = (Pixel*)malloc(sizeof(Pixel) * BufferWidth * BufferHeight);
	BufferB = (Pixel*)malloc(sizeof(Pixel) * BufferWidth * BufferHeight);

	assert(BufferA != NULL && BufferB != NULL);

	// Initialize BufferA with (-1,-1), indicating an invalid closest seed.
	// We don't need to initialize BufferB because it will be written to in the first round.
	for (int y = 0; y < BufferHeight; ++y) {
		for (int x = 0; x < BufferWidth; ++x) {
			int idx = (y * BufferWidth) + x;
			BufferA[idx].x = -1;
			BufferA[idx].y = -1;
		}
	}

	// Put the seeds into the first buffer
	for (int i = 0; i < Seeds.size(); ++i) {
		Pixel& p = Seeds[i];
		BufferA[(p.y * BufferWidth) + p.x] = p;
	}

	// Initial step length is half the image's size. If the image isn't square,
	// we use the largest dimension.
	int step = BufferWidth > BufferHeight ? BufferWidth / 2 : BufferHeight / 2;

	// We use this boolean to know which buffer we are reading from
	ReadingBufferA = true;

	// We read from the RBuffer and write into the WBuffer
	Pixel* RBuffer;
	Pixel* WBuffer;

	// Carry out the rounds of Jump Flooding
	while (step >= 1) {

		// Set which buffers we'll be using
		if (ReadingBufferA == true) {
			RBuffer = BufferA;
			WBuffer = BufferB;
		}
		else {
			RBuffer = BufferB;
			WBuffer = BufferA;
		}

		// Iterate over each Pixel to find its closest seed
		for (int y = 0; y < BufferHeight; ++y) {					//ѭ��ÿһ����
			for (int x = 0; x < BufferWidth; ++x) {

				// The Pixel's absolute index in the buffer
				int idx = (y * BufferWidth) + x;

				// The Pixel's current closest seed (if any)
				Pixel& p = RBuffer[idx];

				// Go ahead and write our current closest seed, if any. If we don't do this
				// we might lose this information if we don't update our seed this round.
				WBuffer[idx] = p;

				// This is a seed, so skip this Pixel
				if (p.x == x && p.y == y)
					continue;

				// This variable will be used to judge which seed is closest
				float dist;

				if (p.x == -1 || p.y == -1)
					dist = -1; // No closest seed has been found yet   
				else
					dist = (p.x - x)*(p.x - x) + (p.y - y)*(p.y - y); // Current closest seed's distance     
																	  //��ǰ�㵽��ǰ���ѱ�ǵ���ʱ���ڽ���ľ���
																	  // To find each Pixel's closest seed, we look at its 8 neighbors thusly:
																	  //   (x-step,y-step) (x,y-step) (x+step,y-step)
																	  //   (x-step,y     ) (x,y     ) (x+step,y     )
																	  //   (x-step,y+step) (x,y+step) (x+step,y+step)

				for (int ky = -1; ky <= 1; ++ky) {
					for (int kx = -1; kx <= 1; ++kx) {				//��ǰ�㣨�ѱ���ǣ�����������

																	// Calculate neighbor's row and column
						int ny = y + ky * step;
						int nx = x + kx * step;

						// If the neighbor is outside the bounds of the buffer, skip it
						if (nx < 0 || nx >= BufferWidth || ny < 0 || ny >= BufferHeight)
							continue;

						// Calculate neighbor's absolute index
						int nidx = (ny * BufferWidth) + nx;

						// Retrieve the neighbor
						Pixel& pk = RBuffer[nidx];

						// If the neighbor doesn't have a closest seed yet, skip it
						if (pk.x == -1 || pk.y == -1)
							continue;

						//���ֵ�ǰ����ϣ������ң��ĵ��Ѿ����˱��

						// Calculate the distance from us to the neighbor's closest seed
						float newDist = (pk.x - x)*(pk.x - x) + (pk.y - y)*(pk.y - y);		//�µľ�����  ��ǰ�㣨x��y�� �� �ٽ���pk

																							// If dist is -1, it means we have no closest seed, so we might as well take this one

																							//��� ��ǰ��û�б�� �ͽ� ���б�ǵ���ϣ������ң���Ϊ�õ�����ٽ���
																							//��� ��ǰ�����б�� �ҵ�ǰ�㵽���б�ǣ����ٽ������ϣ������ң��ľ��� < ��ǰ�㵽���б�ǵ�ľ��� ����


																							// Otherwise, only adopt this new seed if it's closer than our current closest seed

						if (dist == -1 || newDist < dist) {
							WBuffer[idx] = pk;
							dist = newDist;
						}
					}
				}
				distmap[idx] = int(dist);
			}
		}
		// Halve the step.
		step /= 2;

		// Swap the buffers for the next round
		ReadingBufferA = !ReadingBufferA;
	}

	int minlen = 10;
	for (int i = 0; i < xsize*ysize; i++) {

		distmap[i] = int(sqrt(distmap[i]));
		int tparm = int(distmap[i] / (minlen));
		if (!(tparm % 2)) {

			distmap[i] = (distmap[i] % minlen) * 255 / minlen;
		}
		else {

			distmap[i] = (1 - (distmap[i] % minlen)) * 255 / minlen;
		}
	}

	ClearBuffers(BufferA, BufferB);
}


// Jump Flooding Algorithm
void ExecuteJumpFloodingVoi(VoiPixelMap &VoiMap, const EdgeChain Seeds, unsigned int xsize, unsigned int ysize) {

	// Buffer dimensions
	int BufferWidth = xsize;
	int BufferHeight = ysize;
	// Buffers
	Pixel* BufferA = NULL;
	Pixel* BufferB = NULL;

	// Which buffer are we reading from?
	bool ReadingBufferA = true;

	// No seeds will just give us a black screen :P
	if (Seeds.size() < 1) {
		printf("Please create at least 1 seed.\n");
	}

	printf("Executing the Jump Flooding algorithm for Voi...\n");

	// Clear the buffers before we start
	ClearBuffers(BufferA, BufferB);

	// Allocate memory for the two buffers
	BufferA = (Pixel*)malloc(sizeof(Pixel) * BufferWidth * BufferHeight);
	BufferB = (Pixel*)malloc(sizeof(Pixel) * BufferWidth * BufferHeight);


	//assert(BufferA != NULL && BufferB != NULL);

	// Initialize BufferA with (-1,-1), indicating an invalid closest seed.
	// We don't need to initialize BufferB because it will be written to in the first round.
	for (int y = 0; y < BufferHeight; ++y) {
		for (int x = 0; x < BufferWidth; ++x) {
			int idx = (y * BufferWidth) + x;
			BufferA[idx].x = -1;
			BufferA[idx].y = -1;
		}
	}

	// Put the seeds into the first buffer
	for (int i = 0; i < Seeds.size(); ++i) {
		Pixel p = Seeds[i];
		BufferA[(p.y * BufferWidth) + p.x] = p;
	}

	// Initial step length is half the image's size. If the image isn't square,
	// we use the largest dimension.
	int step = BufferWidth > BufferHeight ? BufferWidth / 2 : BufferHeight / 2;

	// We use this boolean to know which buffer we are reading from
	ReadingBufferA = true;

	// We read from the RBuffer and write into the WBuffer
	Pixel* RBuffer;
	Pixel* WBuffer;

	// Carry out the rounds of Jump Flooding
	while (step >= 1) {

		// Set which buffers we'll be using
		if (ReadingBufferA == true) {
			RBuffer = BufferA;
			WBuffer = BufferB;
		}
		else {
			RBuffer = BufferB;
			WBuffer = BufferA;
		}

		// Iterate over each Pixel to find its closest seed
		for (int y = 0; y < BufferHeight; ++y) {					//ѭ��ÿһ����
			for (int x = 0; x < BufferWidth; ++x) {

				// The Pixel's absolute index in the buffer
				int idx = (y * BufferWidth) + x;

				// The Pixel's current closest seed (if any)
				Pixel& p = RBuffer[idx];

				// Go ahead and write our current closest seed, if any. If we don't do this
				// we might lose this information if we don't update our seed this round.
				WBuffer[idx] = p;

				// This is a seed, so skip this Pixel
				if (p.x == x && p.y == y)
					continue;

				// This variable will be used to judge which seed is closest
				float dist;

				if (p.x == -1 || p.y == -1)
					dist = -1; // No closest seed has been found yet   
				else
					dist = (p.x - x)*(p.x - x) + (p.y - y)*(p.y - y); // Current closest seed's distance     
																	  //��ǰ�㵽��ǰ���ѱ�ǵ���ʱ���ڽ���ľ���
																	  // To find each Pixel's closest seed, we look at its 8 neighbors thusly:
																	  //   (x-step,y-step) (x,y-step) (x+step,y-step)
																	  //   (x-step,y     ) (x,y     ) (x+step,y     )
																	  //   (x-step,y+step) (x,y+step) (x+step,y+step)

				for (int ky = -1; ky <= 1; ++ky) {
					for (int kx = -1; kx <= 1; ++kx) {				//��ǰ�㣨�ѱ���ǣ�����������

																	// Calculate neighbor's row and column
						int ny = y + ky * step;
						int nx = x + kx * step;

						// If the neighbor is outside the bounds of the buffer, skip it
						if (nx < 0 || nx >= BufferWidth || ny < 0 || ny >= BufferHeight)
							continue;

						// Calculate neighbor's absolute index
						int nidx = (ny * BufferWidth) + nx;

						// Retrieve the neighbor
						Pixel& pk = RBuffer[nidx];

						// If the neighbor doesn't have a closest seed yet, skip it
						if (pk.x == -1 || pk.y == -1)
							continue;

						//���ֵ�ǰ����ϣ������ң��ĵ��Ѿ����˱��

						// Calculate the distance from us to the neighbor's closest seed
						float newDist = (pk.x - x)*(pk.x - x) + (pk.y - y)*(pk.y - y);		//�µľ�����  ��ǰ�㣨x��y�� �� �ٽ���pk

																							// If dist is -1, it means we have no closest seed, so we might as well take this one

																							//��� ��ǰ��û�б�� �ͽ� ���б�ǵ���ϣ������ң���Ϊ�õ�����ٽ���
																							//��� ��ǰ�����б�� �ҵ�ǰ�㵽���б�ǣ����ٽ������ϣ������ң��ľ��� < ��ǰ�㵽���б�ǵ�ľ��� ����


																							// Otherwise, only adopt this new seed if it's closer than our current closest seed

						if (dist == -1 || newDist < dist) {
							WBuffer[idx] = pk;
							dist = newDist;
						}

					}
				}
			}
		}
		// Halve the step.
		step /= 2;

		// Swap the buffers for the next round
		ReadingBufferA = !ReadingBufferA;
	}

	clock_t begin = clock();
	Pixel* Buffer = (ReadingBufferA == true) ? BufferA : BufferB;
	for (int y = 0; y < BufferHeight; ++y) {

		for (int x = 0; x < BufferWidth; ++x) {

			int idx = (y * BufferWidth) + x;
			int seedIdx = Buffer[idx].y * BufferWidth + Buffer[idx].x;
			VoiMap.insert(make_pair(seedIdx, idx));
		}
	}
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "make_pair" << elapsed_secs << endl;
	ClearBuffers(BufferA, BufferB);
}