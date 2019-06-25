
#ifndef __EDGE_DRAWING_LINE_H_
#define __EDGE_DRAWING_LINE_H_

#include <vector>
#include <list>
#include <array>
using namespace std;
//#include "image_defines.h"
//#ifdef __cplusplus
//extern "C" {
//#endif


typedef struct
{
	int x;
	int y;
	int width;
	int height;
}boundingbox_t;


typedef struct
{
	float startx;
	float starty;
	float endx;
	float endy;
}line_float_t;


/*
@function    EdgeDrawingLineDetector
@param       [in]      src:						  image,single channel
@param       [in]      w:                         width of image
@param       [in]      h:                         height of image
@param       [in]      scaleX:                    downscale factor in X-axis
@param       [in]      scaleY:                    downscale factor in Y-axis
@param       [in]      bbox:                      boundingbox to detect
@param       [in/out]  lines:                     result
@return£º										  0:ok; 1:error
@brief£º


/*
struct define
*/
typedef struct
{
	float u;            //col of pixel
	float v;            //row of pixel
}pixel_float_t;

typedef struct image_int8u_s
{
	unsigned char * data;
	unsigned int xsize, ysize;
} *image_int8u_p;

typedef struct image_int16s_s
{
	short * data;
	unsigned int xsize, ysize;
} *image_int16s_p;

typedef struct image_int32s_s
{
	int * data;
	unsigned int xsize, ysize;
} *image_int32s_p;

typedef struct image_float_s
{
	float * data;
	unsigned int xsize, ysize;
} *image_float_p;

struct SingleLineInfo
{
	/*endPoints, the coordinate origin is the top-left corner of the original image.
	*startPointX = sPointInScaledX * (factor)^scaledCount;	*/
	float startPointX;
	float startPointY;
	float endPointX;
	float endPointY;
	//endPoints, the coordinate origin is the top-left corner of the scaled image.
	float sPointInScaledX;
	float sPointInScaledY;
	float ePointInScaledX;
	float ePointInScaledY;
	//direction of a line, the angle between positive line direction (dark side is in the left) and positive X axis.
	float direction;
	//the summation of gradient magnitudes of pixels on scaled lines
	float salience;
	//the length of line
	float lineLength;
	//number of pixels
	unsigned int numOfPixels;
	//the scaled which this line is detected
	unsigned int scaledCount;
	//the decriptor of line
	std::vector<float> descriptor;
};

// Specifies a vector of lines.
typedef std::vector<SingleLineInfo> LineSet;

typedef std::vector<LineSet> ScaleLineSet;//each element in ScaleLineSet is a vector of lines which corresponds the same line detected in different scaled images.


struct ScaledLine {
	unsigned int scaledCount;//the scaled which this line is detected
	unsigned int lineIDInScaled;//the line ID in that scaled image
	unsigned int lineIDInScaleLineVec;//the line ID in Scale line vector
	float lineLength; //the length of line in original image scale
};


struct Pixel {
	unsigned int x;//X coordinate
	unsigned int y;//Y coordinate

	Pixel() {};
	Pixel(int tx, int ty) {
		x = tx;
		y = ty;
	}

	bool operator == (const Pixel &pp) const{

		if (this->x == pp.x && this->y == pp.y)
			return true;
		else
			return false;
	}

	bool operator < (const Pixel &pp) const{

		if (this->x < pp.x && this->y < pp.y)
			return true;
		else
			return false;
	}
};

struct PixelHash {
	
	std::size_t operator () (const Pixel &t) const {
		return  t.x * 100 + t.y;
	}
};

struct EdgeChains {
	std::vector<unsigned int> xCors;//all the x coordinates of edge points
	std::vector<unsigned int> yCors;//all the y coordinates of edge points
	std::vector<unsigned int> sId;  //the start index of each edge in the coordinate arrays
	unsigned int numOfEdges;//the number of edges whose length are larger than minLineLen; numOfEdges < sId.size;
};
struct LineChains {
	std::vector<unsigned int> xCors;//all the x coordinates of line points
	std::vector<unsigned int> yCors;//all the y coordinates of line points
	std::vector<unsigned int> sId;  //the start index of each line in the coordinate arrays
	unsigned int numOfLines;//the number of lines whose length are larger than minLineLen; numOfLines < sId.size;
};

//typedef  std::list<Pixel> PixelChain;//each edge is a pixel chain


struct EDLineParam {

	float gradientThreshold;
	float anchorThreshold;
	int scanIntervals;
	int minLineLen;
	float lineFitErrThreshold;
};

typedef enum _ORIENT_CODE
{
	ORIENT_HORIZONAL = 1,       // horizotal
	ORIENT_VERTICAL = 2,        // vertical


}ORIENT_CODE;

typedef struct ntuple_list_s
{
	unsigned int size;
	unsigned int max_size;
	unsigned int dim;
	float * values;
} *ntuple_list;


/*
class define
*/
class EDLineDetector
{
public:
	EDLineDetector();
	EDLineDetector(EDLineParam param);
	~EDLineDetector();

public:

	/*extract edges from image
	*image:    In, gray image;
	*edges:    Out, store the edges, each edge is a pixel chain
	*smoothed: In, flag to mark whether the image has already been smoothed by Gaussian filter.
	*return 1: error happen
	*/
	int EdgeDrawing(image_int8u_p image, EdgeChains &edgeChains, bool smoothed = false);

	/*extract lines from image
	*image:    In, gray image;
	*lines:    Out, store the extracted lines,
	*smoothed: In, flag to mark whether the image has already been smoothed by Gaussian filter.
	*return 1: error happen
	*/
	int EDline(image_int8u_p image, LineChains &lines, bool smoothed = false);

	/*extract line from image, and store them*/
	int EDline(image_int8u_p image, bool smoothed = false);


public:
	image_int16s_p dxImg_;
	image_int16s_p dyImg_;
	//store the gradient image without threshold;
	image_int16s_p gImgWO_;

	LineChains lines_; //store the detected line chains;
					   //store the line Equation coefficients, vec3=[w1,w2,w3] for line w1*x + w2*y + w3=0;
	std::vector<std::array<float, 3> > lineEquations_;
	//store the line endpoints, [x1,y1,x2,y3]
	std::vector<std::array<float, 4> > lineEndpoints_;
	//store the line direction
	std::vector<float>  lineDirection_;
	//store the line salience, which is the summation of gradients of pixels on line
	std::vector<float>  lineSalience_;
	unsigned int imageWidth;
	unsigned int imageHeight;

private:
	void InitEDLine_();
	/*For an input edge chain, find the best fit line, the default chain length is minLineLen_
	*xCors:  In, pointer to the X coordinates of pixel chain;
	*yCors:  In, pointer to the Y coordinates of pixel chain;
	*offsetS:In, start index of this chain in array;
	*lineEquation: Out, [a,b] which are the coefficient of lines y=ax+b(horizontal) or x=ay+b(vertical);
	*return:  line fit error; 1:error happens;
	*/
	float LeastSquaresLineFit_(unsigned int *xCors, unsigned int *yCors,
		unsigned int offsetS, std::array<float, 2> &lineEquation);
	/*For an input pixel chain, find the best fit line. Only do the update based on new points.
	*For A*x=v,  Least square estimation of x = Inv(A^T * A) * (A^T * v);
	*If some new observations are added, i.e, [A; A'] * x = [v; v'],
	*then x' = Inv(A^T * A + (A')^T * A') * (A^T * v + (A')^T * v');
	*xCors:  In, pointer to the X coordinates of pixel chain;
	*yCors:  In, pointer to the Y coordinates of pixel chain;
	*offsetS:In, start index of this chain in array;
	*newOffsetS: In, start index of extended part;
	*offsetE:In, end index of this chain in array;
	*lineEquation: Out, [a,b] which are the coefficient of lines y=ax+b(horizontal) or x=ay+b(vertical);
	*return:  line fit error; 1:error happens;
	*/
	float LeastSquaresLineFit_(unsigned int *xCors, unsigned int *yCors,
		unsigned int offsetS, unsigned int newOffsetS,
		unsigned int offsetE, std::array<float, 2> &lineEquation);
	/* Validate line based on the Helmholtz principle, which basically states that
	* for a structure to be perceptually meaningful, the expectation of this structure
	* by chance must be very low.
	*/
	bool LineValidation_(unsigned int *xCors, unsigned int *yCors,
		unsigned int offsetS, unsigned int offsetE,
		std::array<float, 3> &lineEquation, float &direction);


private:

	bool bValidate_;//flag to decide whether line will be validated

					/*the threshold of pixel gradient magnitude.
					*Only those pixel whose gradient magnitude are larger than this threshold will be
					*taken as possible edge points. Default value is 36*/
	short gradienThreshold_;
	/*If the pixel's gradient value is bigger than both of its neighbors by a
	*certain threshold (ANCHOR_THRESHOLD), the pixel is marked to be an anchor.
	*Default value is 8*/
	unsigned char anchorThreshold_;
	/*anchor testing can be performed at different scan intervals, i.e.,
	*every row/column, every second row/column etc.
	*Default value is 2*/
	unsigned int scanIntervals_;
	int minLineLen_;//minimal acceptable line length

					/*This type of storage order is because of the order of edge detection process.
					*For each edge, start from one anchor point, first go right, then go left or first go down, then go up*/
	unsigned int *pFirstPartEdgeX_;//store the X coordinates of the first part of the pixels for chains
	unsigned int *pFirstPartEdgeY_;//store the Y coordinates of the first part of the pixels for chains
	unsigned int *pFirstPartEdgeS_;//store the start index of every edge chain in the first part arrays
	unsigned int *pSecondPartEdgeX_;//store the X coordinates of the second part of the pixels for chains
	unsigned int *pSecondPartEdgeY_;//store the Y coordinates of the second part of the pixels for chains
	unsigned int *pSecondPartEdgeS_;//store the start index of every edge chain in the second part arrays
	unsigned int *pAnchorX_;//store the X coordinates of anchors
	unsigned int *pAnchorY_;//store the Y coordinates of anchors

	image_int8u_p edgeImage_;

	float lineFitErrThreshold_;

	//store the gradient image;
	image_int16s_p gImg_;
	//store the direction image
	image_int8u_p dirImg_;

	float logNT_;

	image_float_p ATA;	//the previous matrix of A^T * A;
	image_float_p ATV;	//the previous vector of A^T * V;
	image_float_p fitMatT;	//the matrix used in line fit function;
	image_float_p fitVec;	//the vector used in line fit function;
	image_float_p tempMatLineFit;	//the matrix used in line fit function;
	image_float_p tempVecLineFit;	//the vector used in line fit function;

	int float_equal(float a, float b);
	float log_gamma_lanczos(float x);
	float log_gamma_windschitl(float x);
	float nfa(int n, int k, float p, float  logNT);
};


/* This class is used to generate the line descriptors from multi-scale images  */
class LineDescriptor
{
public:
	LineDescriptor();
	LineDescriptor(unsigned int numOfBand, unsigned int widthOfBand);
	~LineDescriptor();

public:

	/* Interface.*/
	int Run(float scaleX, float scaleY, boundingbox_t bbox,
		image_int8u_p image, ScaleLineSet & keyLines);

private:

	/*This function is used to detect lines from multi-scale images.*/
	int ScaledKeyLines(image_int8u_p image, ScaleLineSet &keyLines);

	/*This function is used to get numbers of pixels in a line from image.*/
	int GetLinePixelsNums(float startX, float startY, float endX, float endY);

	/*This function is used to get information of lines from downsampled image.*/
	void InverseGaussianSamplerLines(pixel_float_t gs_scale, ScaleLineSet &keyLines);



private:

	/*For each scaled of image, we define an EDLineDetector, because we can get gradient images (dxImg, dyImg, gImg)
	*from the EDLineDetector class without extra computation cost. Another reason is that, if we use
	*a single EDLineDetector to detect lines in different scaled of images, then we need to allocate and release
	*memory for gradient images (dxImg, dyImg, gImg) repeatedly for their varying size*/
	std::vector<EDLineDetector*> edLineVec_;

	//int ksize_; //the size of Gaussian kernel: ksize X ksize, default value is 5.

	unsigned int  numOfBand_;//the number of band used to compute line descriptor
	unsigned int  widthOfBand_;//the width of band;
	std::vector<float> gaussCoefL_;//the local gaussian coefficient apply to the orthogonal line direction within each band;
	std::vector<float> gaussCoefG_;//the global gaussian coefficient apply to each Row within line support region

};


static image_int8u_p new_image_int8u_ptr(unsigned int xsize, unsigned int ysize, unsigned char * data);
static image_int8u_p new_image_int8u(unsigned int xsize, unsigned int ysize);
static image_int16s_p new_image_int16s(unsigned int xsize, unsigned int ysize);
static image_float_p new_image_float(unsigned int xsize, unsigned int ysize);
void free_image_float(image_float_p i);
void free_image_int16s(image_int16s_p i);
static void free_image_int32s(image_int32s_p i);
static void free_image_int8u(image_int8u_p i);
static void sobel_edge(ORIENT_CODE oriention, image_int8u_p src, image_int16s_p dst);
static void mcv_sobel(ORIENT_CODE oriention, image_int8u_p src, image_int16s_p dst);
static void array_abs(image_int16s_p src, image_int16s_p dst);
static void mcv_abs(image_int16s_p src, image_int16s_p dst);
static void array_add(image_int16s_p src1, image_int16s_p src2, image_int16s_p dst);
static int array_add(image_float_p src1, image_float_p src2, image_float_p dst);
static void mcv_add(image_int16s_p src1, image_int16s_p src2, image_int16s_p dst);
static void mcv_add(image_float_p src1, image_float_p src2, image_float_p dst);
static image_int16s_p array_threshold(short thresh, image_int16s_p src);
static void array_threshold(short thresh, image_int16s_p src, image_int16s_p dst);
static void mcv_threshold(short thresh, image_int16s_p src, image_int16s_p dst);
static void array_compare_lt(image_int16s_p src1, image_int16s_p src2, image_int8u_p dst);
static void mcv_compare_CMP_LT(image_int16s_p src1, image_int16s_p src2, image_int8u_p dst);
static image_int16s_p array_devide(int n, image_int16s_p src);
static void array_devide(int n, image_int16s_p src, image_int16s_p dst);
static void mcv_mat_divide(int n, image_int16s_p src, image_int16s_p dst);
static int array_multiply_transpose_float(image_float_p src, image_float_p dst);
static int array_multiply2_transpose_float(image_float_p src1, image_float_p src2, image_float_p dst);
static int array_multiply(image_float_p src1, image_float_p src2, image_float_p dst);
static void mcv_multiply_float(image_float_p src1, image_float_p src2, image_float_p dst);
static void mcv_multiply_transpose_float(image_float_p src, image_float_p dst);
static void mcv_multiply2_transpose_float(image_float_p src1, image_float_p src2, image_float_p dst);
static void enlarge_ntuple_list(ntuple_list n_tuple);
static ntuple_list new_ntuple_list(unsigned int dim);
static void free_ntuple_list(ntuple_list in);
static void gaussian_kernel(ntuple_list kernel, float sigma, float mean);
static image_int8u_p gaussian_sampler_byte_bbox(image_int8u_p in, boundingbox_t bbox, pixel_float_t scale, float sigma_scale);
static image_int8u_p gaussian_sampler_byte(image_int8u_p in, pixel_float_t scale, float sigma_scale);
static void InverseBoundingBoxLines(boundingbox_t bbox, ScaleLineSet & keyLines);
int _edge_drawing_line_detector(unsigned char *src, int w, int h,
	float scaleX, float scaleY, boundingbox_t bbox, std::vector<line_float_t> &lines);
int EdgeDrawingLineDetector(unsigned char *src, int w, int h,
	float scaleX, float scaleY, boundingbox_t bbox, std::vector<line_float_t> &lines);



EdgeChains GetEdgeChains(unsigned char* data, unsigned int xsize, unsigned int ysize);

//
//#ifdef __cplusplus
//}
//#endif

#endif