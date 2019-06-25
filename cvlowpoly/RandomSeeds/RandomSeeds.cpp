#include "RandomSeeds.h"

Seeds creatRandomSeeds(int xsize, int ysize, int NumPoints) {
	srand((int)time(0));
	Seeds RandomSeeds;
	for (int k = 0; k < NumPoints; k++) {

		int x = random(xsize);
		int y = random(ysize);
		Pixel point(x, y);
		RandomSeeds.push_back(point);
	}

	return RandomSeeds;
}

