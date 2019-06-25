
#ifndef __IMAGELIB_H__
#define __IMAGELIB_H__


#include <opencv2/core/utility.hpp>  
#include <opencv2/opencv.hpp>
using namespace cv;


struct ImgSet {
	Mat InitImg;
	Mat GrayImg;
	Mat FilterImg;
};


Mat GetInitImg(String filename);
                                                                                                
Mat GetGrayImg(Mat InitImg);

Mat GetFilterImg(Mat GrayImg);

uchar* GetImgData(Mat Img);

uint Getxsize(Mat Img);

uint Getysize(Mat Img);


#endif