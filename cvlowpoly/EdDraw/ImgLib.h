
#ifndef _IMAGELIB_H_
#define _IMAGELIB_H_


#include <opencv2/core/utility.hpp>  
#include <opencv2/opencv.hpp>
using namespace cv;


struct ImgSet {
	Mat InitImg;
	Mat GrayImg;
	Mat FilterImg;
};


Mat GetInitImg(String filename) {

	Mat InitImg;
	InitImg = imread(filename, 1);
	if (InitImg.empty()) {
		printf("Error:load image failed.");
		exit(1);
	}
	return InitImg;
	
}

Mat GetGrayImg(Mat InitImg) {

	Mat GrayImg;
	cvtColor(InitImg, GrayImg, CV_BGR2GRAY);
	return GrayImg;

}

Mat GetFilterImg(Mat GrayImg) {

	Mat FilterImg;
	GaussianBlur(GrayImg, FilterImg, Size(5, 5), 0, 0);
	return FilterImg;

}

uchar* GetImgData(Mat Img) {
	return Img.data;
}

uint Getxsize(Mat Img) {
	return Img.cols;
}

uint Getysize(Mat Img) {
	return Img.rows;
}


#endif