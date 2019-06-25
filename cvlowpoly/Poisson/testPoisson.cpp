#include "Poisson.h"
#include "JFA.h"

//int main(int argc, char** argv) {
//
/////////////////// User selectable parameters ///////////////////////////////
//
//	const int   NumPoints = 200;	// minimal number of points to generate
//	//const int   ImageSize = 1024;	// generate RGB image [ImageSize x ImageSize]
//
//////////////////////////////////////////////////////////////////////////////
//	DefaultPRNG PRNG;
//
//	//const auto Points = GeneratePoissonPoints(NumPoints, PRNG); 
//	const auto Points = GeneratePoissonPoints(NumPoints, PRNG, 10, false, -1.0f);
//
//
//
//	// dump points to a text file
//	std::ofstream File("Poisson.txt", std::ios::out);
//
//	File << "NumPoints = " << Points.size() << std::endl;
//
//	for (const auto& p : Points)
//	{
//		File << p.x << "\t"  << p.y << std::endl;
//	}
//
//	return 0;
//}
