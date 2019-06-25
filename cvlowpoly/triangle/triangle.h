#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include <opencv2/core/utility.hpp>  
#include <opencv2/opencv.hpp>
using namespace cv;

void draw_point(Mat& img, Point2f fp, Scalar color);
void draw_delaunay(Mat& img, Subdiv2D& subdiv, Scalar delaunay_color);


#endif