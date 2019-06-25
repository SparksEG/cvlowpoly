
#ifndef _ED_H_
#define _ED_H_

#include "edlines.h"
#include <list>
#include <array>



#ifndef PI
#define PI (3.1415926535897932384626433832795)
#endif

#ifndef ZERO
#define ZERO (0)
#endif

#ifndef ROUND
#define ROUND (0.5F)
#endif

#ifndef MIN
#define MIN(a,b)  ((a) > (b) ? (b) : (a))
#endif

#ifndef MAX
#define MAX(a,b)  ((a) < (b) ? (b) : (a))
#endif

//if |dx|<|dy|;
#ifndef Horizontal
#define Horizontal (255)
#endif
//if |dy|<=|dx|;
#ifndef Vertical
#define Vertical    0
#endif

#ifndef UpDir
#define UpDir       1
#endif

#ifndef RightDir
#define RightDir    2
#endif

#ifndef DownDir
#define DownDir     3
#endif

#ifndef LeftDir
#define LeftDir     4
#endif

#ifndef TryTime
#define TryTime     6
#endif

#ifndef SkipEdgePoint
#define SkipEdgePoint 2
#endif

#ifndef RELATIVE_ERROR_FACTOR
#define RELATIVE_ERROR_FACTOR   100.0f
#endif

#ifndef M_LN10
#define M_LN10   2.302585093f 
#endif

#ifndef log_gamma
#define log_gamma(x)    ((x)>15.f?log_gamma_windschitl(x):log_gamma_lanczos(x))
#endif

#ifndef SalienceScale
#define SalienceScale 0.9F//0.9
#endif

#ifndef ONE
#define ONE (1)
#endif


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

struct EdgeChains {
	std::vector<unsigned int> xCors;//all the x coordinates of edge points
	std::vector<unsigned int> yCors;//all the y coordinates of edge points
	std::vector<unsigned int> sId;  //the start index of each edge in the coordinate arrays
	unsigned int numOfEdges;//the number of edges whose length are larger than minLineLen; numOfEdges < sId.size;
};


struct EDLineParam {

	float gradientThreshold;
	float anchorThreshold;
	int scanIntervals;
	int minLineLen;
	float lineFitErrThreshold;
};


static image_int8u_p new_image_int8u_ptr(unsigned int xsize,
	unsigned int ysize, unsigned char * data)
{
	image_int8u_p image = NULL;

	/* get memory */
	image = new image_int8u_s[1];
	if (!image)exit(1);

	/* set image */
	image->xsize = xsize;
	image->ysize = ysize;
	image->data = data;

	return image;
}


static image_int8u_p new_image_int8u(unsigned int xsize, unsigned int ysize)
{
	image_int8u_p image = NULL;

	/* get memory */
	image = new image_int8u_s[1];
	image->data = new unsigned char[xsize*ysize];

	/* set image size */
	image->xsize = xsize;
	image->ysize = ysize;

	return image;
}

static image_int16s_p new_image_int16s(unsigned int xsize, unsigned int ysize)
{
	image_int16s_p image = NULL;

	/* get memory */
	image = new image_int16s_s[1];
	image->data = new short[xsize*ysize];

	/* set image size */
	image->xsize = xsize;
	image->ysize = ysize;

	return image;
}



static image_float_p new_image_float(unsigned int xsize, unsigned int ysize)
{
	image_float_p image = NULL;

	/* get memory */
	image = new image_float_s[1];
	image->data = new float[xsize*ysize];

	/* set image size */
	image->xsize = xsize;
	image->ysize = ysize;

	return image;
}

//修改
static void free_image_float(image_float_p i)
{
	delete[]i->data;
	delete[]i;
}

//修改
static void free_image_int16s(image_int16s_p i)
{
	delete[]i->data;
	delete[]i;
}


static void free_image_int32s(image_int32s_p i)
{
	delete[]i->data;
	delete[]i;
}

static void free_image_int8u(image_int8u_p i)
{
	delete[]i->data;
	delete[]i;
}


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
	int EdgeDrawing(image_int8u_p image, EdgeChains &edgeChains);

public:
	image_int16s_p dxImg_;
	image_int16s_p dyImg_;
	//store the gradient image without threshold;
	image_int16s_p gImgWO_;

	unsigned int imageWidth;
	unsigned int imageHeight;

private:
	void InitEDLine_();


private:

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


									/** Compare floats by relative error.*/
	static int float_equal(float a, float b)
	{
		float abs_diff, aa, bb, abs_max;
		/* trivial case */
		if (a == b) return true;
		abs_diff = fabs(a - b);
		aa = fabs(a);
		bb = fabs(b);
		abs_max = aa > bb ? aa : bb;

		if (abs_max < FLT_MIN) abs_max = FLT_MIN;
		/* equal if relative error <= factor x eps */
		return (abs_diff / abs_max) <= (RELATIVE_ERROR_FACTOR * DBL_EPSILON);
	}
	/** Computes the natural logarithm of the absolute value of
	the gamma function of x using the Lanczos approximation.*/
	static float log_gamma_lanczos(float x)
	{
		static float q[7] = { 75122.6331530f, 80916.6278952f, 36308.2951477f,
			8687.24529705f, 1168.92649479f, 83.8676043424f,
			2.50662827511f };
		float a = (x + 0.5f) * log(x + 5.5f) - (x + 5.5f);
		float b = 0.f;
		int n;
		for (n = 0; n<7; n++) {
			a -= log(x + (float)n);
			b += q[n] * pow(x, (float)n);
		}
		return a + log(b);
	}
	/** Computes the natural logarithm of the absolute value of
	the gamma function of x using Windschitl method.*/
	static float log_gamma_windschitl(float x)
	{
		return 0.918938533204673f + (x - 0.5f)*log(x) - x
			+ 0.5f*x*log(x*sinh(1 / x) + 1 / (810.f*pow(x, 6.f)));
	}
	/** Computes -log10(NFA).*/
	static float nfa(int n, int 	k, float p, float  logNT)
	{
		float tolerance = 0.1f;       /* an error of 10% in the result is accepted */
		float log1term, term, bin_term, mult_term, bin_tail, err, p_term;
		int i;

		/* check parameters */
		if (n<0 || k<0 || k>n || p <= 0.f || p >= 1.f) {
			printf("nfa: wrong n, k or p values.\n");
			exit(1);
		}
		/* trivial cases */
		if (n == 0 || k == 0) return -logNT;
		if (n == k) return -logNT - (float)n * log10(p);

		/* probability term */
		p_term = p / (1.f - p);

		/* compute the first term of the series */
		log1term = log_gamma((float)n + 1.f) - log_gamma((float)k + 1.f)
			- log_gamma((float)(n - k) + 1.f)
			+ (float)k * log(p) + (float)(n - k) * log(1.f - p);
		term = exp(log1term);

		/* in some cases no more computations are needed */
		if (float_equal(term, 0.f)) {  /* the first term is almost zero */
			if ((float)k > (float)n * p)     /* at begin or end of the tail?  */
				return -log1term / M_LN10 - logNT;  /* end: use just the first term  */
			else
				return -logNT;                      /* begin: the tail is roughly 1  */
		}

		/* compute more terms if needed */
		bin_tail = term;
		for (i = k + 1; i <= n; i++) {

			bin_term = (float)(n - i + 1) / (float)i;
			mult_term = bin_term * p_term;
			term *= mult_term;
			bin_tail += term;
			if (bin_term<1.f) {

				err = term * ((1.f - pow(mult_term, (float)(n - i + 1)))
					/ (1.f - mult_term) - 1.f);

				if (err < tolerance * fabs(-log10(bin_tail) - logNT) * bin_tail) break;
			}
		}
		return -log10(bin_tail) - logNT;
	}
};

//修改
inline EDLineDetector::EDLineDetector()
{
	//set parameters for line segment detection

	//gradienThreshold_ = 80; // ***** ORIGINAL WAS 25
	//anchorThreshold_ = 2;//8
	//scanIntervals_ = 2;//2
	//minLineLen_ = 15;
	//lineFitErrThreshold_ = 1.4f;
	gradienThreshold_ = 36; // ***** ORIGINAL WAS 25
	anchorThreshold_ = 8;//8
	scanIntervals_ = 2;//2
	minLineLen_ = 10;
	lineFitErrThreshold_ = 1.4f;

	InitEDLine_();
}
//修改
inline EDLineDetector::EDLineDetector(EDLineParam param)
{
	//set parameters for line segment detection
	gradienThreshold_ = (short)param.gradientThreshold;
	anchorThreshold_ = (unsigned char)param.anchorThreshold;
	scanIntervals_ = param.scanIntervals;
	minLineLen_ = param.minLineLen;
	lineFitErrThreshold_ = param.lineFitErrThreshold;
	InitEDLine_();
}
//修改
inline void EDLineDetector::InitEDLine_()
{

	ATA = new_image_float(2, 2);
	ATV = new_image_float(1, 2);
	tempMatLineFit = new_image_float(2, 2);
	tempVecLineFit = new_image_float(1, 2);
	fitMatT = new_image_float(minLineLen_, 2);
	fitVec = new_image_float(minLineLen_, 1);


	for (int i = 0; i<minLineLen_; i++) {
		fitMatT->data[1 * fitMatT->xsize + i] = 1;
	}

	pFirstPartEdgeX_ = NULL;
	pFirstPartEdgeY_ = NULL;
	pFirstPartEdgeS_ = NULL;
	pSecondPartEdgeX_ = NULL;
	pSecondPartEdgeY_ = NULL;
	pSecondPartEdgeS_ = NULL;
	pAnchorX_ = NULL;
	pAnchorY_ = NULL;


	dirImg_ = NULL;
	gImgWO_ = NULL;
	gImg_ = NULL;
	dxImg_ = NULL;
	dyImg_ = NULL;
	edgeImage_ = NULL;
}
//修改
inline EDLineDetector::~EDLineDetector() {
	if (pFirstPartEdgeX_ != NULL) {
		delete[] pFirstPartEdgeX_; pFirstPartEdgeX_ = NULL;
		delete[] pFirstPartEdgeY_; pFirstPartEdgeY_ = NULL;
		delete[] pSecondPartEdgeX_; pSecondPartEdgeX_ = NULL;
		delete[] pSecondPartEdgeY_; pSecondPartEdgeY_ = NULL;
		delete[] pAnchorX_; pAnchorX_ = NULL;
		delete[] pAnchorY_; pAnchorY_ = NULL;
	}
	if (pFirstPartEdgeS_ != NULL) {
		delete[] pFirstPartEdgeS_; pFirstPartEdgeS_ = NULL;
		delete[] pSecondPartEdgeS_; pSecondPartEdgeS_ = NULL;
	}

	free_image_int8u(edgeImage_);
	free_image_int8u(dirImg_);
	free_image_int16s(gImgWO_);
	free_image_int16s(gImg_);
	free_image_int16s(dxImg_);
	free_image_int16s(dyImg_);

	free_image_float(ATA);
	free_image_float(ATV);
	free_image_float(fitMatT);
	free_image_float(fitVec);
	free_image_float(tempMatLineFit);
	free_image_float(tempVecLineFit);


}


typedef enum _ORIENT_CODE
{
	ORIENT_HORIZONAL = 1,       // horizotal
	ORIENT_VERTICAL = 2,        // vertical


}ORIENT_CODE;

static void sobel_edge(ORIENT_CODE oriention, image_int8u_p src, image_int16s_p dst)
{
	unsigned char *psrc = NULL;
	unsigned int nsize;
	unsigned int i, j, _center, offset_up, offset_down;
	unsigned int _tp, _td, _t;

	nsize = src->xsize * src->ysize;

	//no edge processing
	//memset(dst->data,ZERO,sizeof(short));

	psrc = src->data;
	switch (oriention)
	{
	case ORIENT_HORIZONAL:
		for (i = 1; i < src->ysize - 1; i++)
		{
			_center = i * src->xsize;
			offset_up = _center - src->xsize;
			offset_down = _center + src->xsize;
			for (j = 1; j < src->xsize - 1; j++)
			{
				_tp = offset_up + j;
				_td = offset_down + j;
				_t = _center + j;
				dst->data[_t] = ((short)psrc[_tp + 1] - (short)psrc[_tp - 1])
					+ ((short)psrc[_td + 1] - (short)psrc[_td - 1])
					+ (((short)psrc[_t + 1] - (short)psrc[_t - 1]) << 1);
			}
		}
		break;

	case  ORIENT_VERTICAL:
		for (i = 1; i < src->ysize - 1; i++)
		{
			_center = i * src->xsize;
			offset_up = _center - src->xsize;
			offset_down = _center + src->xsize;
			for (j = 1; j < src->xsize - 1; j++)
			{
				_tp = offset_up + j;
				_td = offset_down + j;
				_t = _center + j;

				dst->data[_t] = -((short)psrc[_tp - 1] + (((short)psrc[_tp]) << 1) + (short)psrc[_tp + 1])
					+ ((short)psrc[_td - 1] + (((short)psrc[_td]) << 1) + (short)psrc[_td + 1]);
			}
		}
		break;

	default:
		printf("sobel oriention is wrong!");
		break;
	}

	psrc = NULL;

}


static void mcv_sobel(ORIENT_CODE oriention, image_int8u_p src, image_int16s_p dst)
{
	sobel_edge(oriention, src, dst);
}

static void array_abs(image_int16s_p src, image_int16s_p dst)
{
	int nsize;
	int k;
	short *psrc = NULL, *pdst = NULL;

	nsize = src->xsize*src->ysize;

	psrc = src->data;
	pdst = dst->data;
	for (k = 0; k < nsize; k++)
	{
		*pdst = (*psrc >= ZERO) ? (*psrc) : (-*psrc);
		psrc++;
		pdst++;
	}

	psrc = NULL;
	pdst = NULL;

}


static void mcv_abs(image_int16s_p src, image_int16s_p dst)
{
	array_abs(src, dst);
}


static int array_add(image_float_p src1, image_float_p src2, image_float_p dst)
{
	int _code = 0;
	int k, nsize;
	float *psrc1 = NULL, *psrc2 = NULL, *pdst = NULL;

	if ((src1 == NULL) || (src2 == NULL) || (dst == NULL)
		|| (src1->xsize != src2->xsize) || (src1->ysize != src2->ysize)
		|| (src1->xsize != dst->xsize) || (src1->ysize != dst->ysize))
		return 1;

	nsize = src1->xsize*src1->ysize;

	psrc1 = src1->data;
	psrc2 = src2->data;
	pdst = dst->data;
	for (k = 0; k < nsize; k++)
	{
		*pdst = *psrc1 + *psrc2;
		pdst++;
		psrc1++;
		psrc2++;
	}

	psrc1 = NULL;
	psrc2 = NULL;
	pdst = NULL;

	return 0;
}

static void array_add(image_int16s_p src1, image_int16s_p src2, image_int16s_p dst)
{
	short *psrc1 = NULL, *psrc2 = NULL, *pdst = NULL;
	int nsize;
	int k;

	nsize = src1->xsize*src1->ysize;

	psrc1 = src1->data;
	psrc2 = src2->data;
	pdst = dst->data;

	for (k = 0; k < nsize; k++)
	{
		*pdst = *psrc1 + *psrc2;
		psrc1++;
		psrc2++;
		pdst++;
	}

}

static void mcv_add(image_int16s_p src1, image_int16s_p src2, image_int16s_p dst)
{
	array_add(src1, src2, dst);
}

static void mcv_add(image_float_p src1, image_float_p src2, image_float_p dst)
{
	array_add(src1, src2, dst);
}


static image_int16s_p array_threshold(short thresh, image_int16s_p src)
{
	image_int16s_p dst = NULL;
	int nsize;
	int k;
	short *psrc = NULL, *pdst = NULL;

	dst = new image_int16s_s[1];

	dst->xsize = src->xsize;
	dst->ysize = src->ysize;
	nsize = src->xsize*src->ysize;

	dst->data = new short[nsize];

	psrc = src->data;
	pdst = dst->data;
	for (k = 0; k<nsize; k++)
	{
		*pdst = (*psrc > thresh) ? (*psrc) : (ZERO);
		psrc++;
		pdst++;
	}

	psrc = NULL;
	pdst = NULL;

	return dst;

}


static void array_threshold(short thresh, image_int16s_p src, image_int16s_p dst)
{
	int nsize;
	int k;
	short *psrc = NULL, *pdst = NULL;

	nsize = src->xsize*src->ysize;

	psrc = src->data;
	pdst = dst->data;
	for (k = 0; k<nsize; k++)
	{
		*pdst = (*psrc > thresh) ? (*psrc) : (ZERO);
		psrc++;
		pdst++;
	}

	psrc = NULL;
	pdst = NULL;


}


static void mcv_threshold(short thresh, image_int16s_p src, image_int16s_p dst)
{
	array_threshold(thresh, src, dst);

}


static void array_compare_lt(image_int16s_p src1, image_int16s_p src2, image_int8u_p dst)
{
	short *psrc1 = NULL, *psrc2 = NULL;
	unsigned char *pdst = NULL;
	int nsize, k;

	psrc1 = src1->data;
	psrc2 = src2->data;
	pdst = dst->data;

	nsize = src1->xsize*src1->ysize;

	for (k = 0; k < nsize; k++)
	{
		*pdst = (*psrc1 < *psrc2) ? (255) : (ZERO);
		pdst++;
		psrc1++;
		psrc2++;
	}

	psrc1 = psrc2 = NULL;
	pdst = NULL;

}

static void mcv_compare_CMP_LT(image_int16s_p src1, image_int16s_p src2, image_int8u_p dst)
{

	array_compare_lt(src1, src2, dst);

}


static image_int16s_p array_devide(int n, image_int16s_p src)
{
	image_int16s_p dst = NULL;
	int nsize;
	int k;
	float n_inv;
	short *psrc = NULL, *pdst = NULL;

	dst = new image_int16s_s[1];

	dst->xsize = src->xsize;
	dst->ysize = src->ysize;
	nsize = src->xsize*src->ysize;

	dst->data = new short[nsize];

	psrc = src->data;
	pdst = dst->data;

	n_inv = 1.f / n;
	for (k = 0; k < nsize; k++)
	{
		*pdst = (short)(n_inv * (*psrc));
		psrc++;
		pdst++;
	}

	psrc = NULL;
	pdst = NULL;

	return dst;
}

static void array_devide(int n, image_int16s_p src, image_int16s_p dst)
{
	int nsize;
	int k;
	short *psrc = NULL, *pdst = NULL;

	nsize = src->xsize*src->ysize;

	psrc = src->data;
	pdst = dst->data;
	for (k = 0; k < nsize; k++)
	{
		*pdst = *psrc / n;
		psrc++;
		pdst++;
	}

	psrc = NULL;
	pdst = NULL;

}


static void mcv_mat_divide(int n, image_int16s_p src, image_int16s_p dst)
{
	array_devide(n, src, dst);
}


//

inline int EDLineDetector::EdgeDrawing(image_int8u_p image, EdgeChains &edgeChains)
{
	imageWidth = image->xsize;
	imageHeight = image->ysize;

	unsigned int pixelNum = imageWidth * imageHeight;
	unsigned int edgePixelArraySize = pixelNum / 5;
	unsigned int maxNumOfEdge = edgePixelArraySize / 20;
	//compute dx, dy images
	if ((gImg_ == NULL) || ((gImg_->xsize != (int)imageWidth) || (gImg_->ysize != (int)imageHeight))) {
		if (pFirstPartEdgeX_ != NULL) {
			delete[] pFirstPartEdgeX_; pFirstPartEdgeX_ = NULL;
			delete[] pFirstPartEdgeY_; pFirstPartEdgeY_ = NULL;
			delete[] pSecondPartEdgeX_; pSecondPartEdgeX_ = NULL;
			delete[] pSecondPartEdgeY_; pSecondPartEdgeY_ = NULL;
			delete[] pFirstPartEdgeS_; pFirstPartEdgeS_ = NULL;
			delete[] pSecondPartEdgeS_; pSecondPartEdgeS_ = NULL;
			delete[] pAnchorX_; pAnchorX_ = NULL;
			delete[] pAnchorY_; pAnchorY_ = NULL;
		}

		dxImg_ = new_image_int16s(imageWidth, imageHeight);
		dyImg_ = new_image_int16s(imageWidth, imageHeight);
		gImgWO_ = new_image_int16s(imageWidth, imageHeight);
		gImg_ = new_image_int16s(imageWidth, imageHeight);
		dirImg_ = new_image_int8u(imageWidth, imageHeight);
		edgeImage_ = new_image_int8u(imageWidth, imageHeight);
		pFirstPartEdgeX_ = new unsigned int[edgePixelArraySize];
		pFirstPartEdgeY_ = new unsigned int[edgePixelArraySize];
		pSecondPartEdgeX_ = new unsigned int[edgePixelArraySize];
		pSecondPartEdgeY_ = new unsigned int[edgePixelArraySize];
		pFirstPartEdgeS_ = new unsigned int[maxNumOfEdge];
		pSecondPartEdgeS_ = new unsigned int[maxNumOfEdge];
		pAnchorX_ = new unsigned int[edgePixelArraySize];
		pAnchorY_ = new unsigned int[edgePixelArraySize];
	}

	mcv_sobel(ORIENT_HORIZONAL, image, dxImg_);
	mcv_sobel(ORIENT_VERTICAL, image, dyImg_);

	//compute gradient and direction images
	image_int16s_p dxABS_m = NULL, dyABS_m = NULL;
	image_int16s_p sumDxDy = NULL;

	dxABS_m = new_image_int16s(dxImg_->xsize, dxImg_->ysize);
	dyABS_m = new_image_int16s(dyImg_->xsize, dyImg_->ysize);
	sumDxDy = new_image_int16s(dyImg_->xsize, dyImg_->ysize);

	mcv_abs(dxImg_, dxABS_m);
	mcv_abs(dyImg_, dyABS_m);

	mcv_add(dyABS_m, dxABS_m, sumDxDy);

	mcv_threshold(gradienThreshold_ + 1, sumDxDy, gImg_);

	mcv_mat_divide(4, gImg_, gImg_);
	mcv_mat_divide(4, sumDxDy, gImgWO_);

	mcv_compare_CMP_LT(dxABS_m, dyABS_m, dirImg_);

	free_image_int16s(sumDxDy);
	free_image_int16s(dxABS_m);
	free_image_int16s(dyABS_m);

	short *pgImg = gImg_->data;

	unsigned char *pdirImg = dirImg_->data;

	//extract the anchors in the gradient image, store into a vector
	memset(pAnchorX_, 0, edgePixelArraySize * sizeof(unsigned int));//initialization
	memset(pAnchorY_, 0, edgePixelArraySize * sizeof(unsigned int));
	unsigned int anchorsSize = 0;
	int offy;
	int indexInArray;
	unsigned char gValue1, gValue2, gValue3;
	for (unsigned int h = 1; h<imageHeight - 1; h = h + scanIntervals_) {
		offy = h * imageWidth;
		for (unsigned int w = 1; w<imageWidth - 1; w = w + scanIntervals_) {
			indexInArray = offy + w;
			if (pdirImg[indexInArray] == Horizontal) {//if the direction of pixel is horizontal, then compare with up and down
				if (pgImg[indexInArray] >= pgImg[indexInArray - imageWidth] + anchorThreshold_
					&& pgImg[indexInArray] >= pgImg[indexInArray + imageWidth] + anchorThreshold_) {// (w,h) is accepted as an anchor
					pAnchorX_[anchorsSize] = w;
					pAnchorY_[anchorsSize++] = h;
				}
			}
			else {
				if (pgImg[indexInArray] >= pgImg[indexInArray - 1] + anchorThreshold_
					&& pgImg[indexInArray] >= pgImg[indexInArray + 1] + anchorThreshold_) {// (w,h) is accepted as an anchor
					pAnchorX_[anchorsSize] = w;
					pAnchorY_[anchorsSize++] = h;
				}
			}
		}
	}
	if (anchorsSize>edgePixelArraySize) {
		printf("anchor size is larger than its maximal size. anchorsSize = %d, maximal size = %d\n",
			anchorsSize, edgePixelArraySize);
		return 1;
	}

	//link the anchors by smart routing
	memset(edgeImage_->data, ZERO, edgeImage_->xsize*edgeImage_->ysize * sizeof(unsigned char));
	unsigned char *pEdgeImg = edgeImage_->data;
	memset(pFirstPartEdgeX_, 0, edgePixelArraySize * sizeof(unsigned int));//initialization
	memset(pFirstPartEdgeY_, 0, edgePixelArraySize * sizeof(unsigned int));
	memset(pSecondPartEdgeX_, 0, edgePixelArraySize * sizeof(unsigned int));
	memset(pSecondPartEdgeY_, 0, edgePixelArraySize * sizeof(unsigned int));
	memset(pFirstPartEdgeS_, 0, maxNumOfEdge * sizeof(unsigned int));
	memset(pSecondPartEdgeS_, 0, maxNumOfEdge * sizeof(unsigned int));
	unsigned int offsetPFirst = 0, offsetPSecond = 0;
	unsigned int offsetPS = 0;

	unsigned int x, y;
	unsigned int lastX, lastY;
	unsigned char lastDirection;//up = 1, right = 2, down = 3, left = 4;
	unsigned char shouldGoDirection;//up = 1, right = 2, down = 3, left = 4;
	int edgeLenFirst, edgeLenSecond;

	lastX = lastY = 0;

	for (unsigned int i = 0; i<anchorsSize; i++) {
		x = pAnchorX_[i];
		y = pAnchorY_[i];
		indexInArray = y * imageWidth + x;
		if (pEdgeImg[indexInArray]) {//if anchor i is already been an edge pixel.
			continue;
		}

		pFirstPartEdgeS_[offsetPS] = offsetPFirst;
		if (pdirImg[indexInArray] == Horizontal) {//if the direction of this pixel is horizontal, then go left and right.
												  //fist go right, pixel direction may be different during linking.
			lastDirection = RightDir;
			while (pgImg[indexInArray]>0 && !pEdgeImg[indexInArray]) {
				pEdgeImg[indexInArray] = 1;        // Mark this pixel as an edge pixel
				pFirstPartEdgeX_[offsetPFirst] = x;
				pFirstPartEdgeY_[offsetPFirst++] = y;
				shouldGoDirection = 0;//unknown
				if (pdirImg[indexInArray] == Horizontal) {//should go left or right 
					if (lastDirection == UpDir || lastDirection == DownDir) {//change the pixel direction now
						if (x>lastX) {//should go right
							shouldGoDirection = RightDir;
						}
						else {//should go left
							shouldGoDirection = LeftDir;
						}
					}
					lastX = x;
					lastY = y;
					if (lastDirection == RightDir || shouldGoDirection == RightDir) {//go right
						if (x == imageWidth - 1 || y == 0 || y == imageHeight - 1) {//reach the image border
							break;
						}
						// Look at 3 neighbors to the right and pick the one with the max. gradient value
						gValue1 = (unsigned char)pgImg[indexInArray - imageWidth + 1];
						gValue2 = (unsigned char)pgImg[indexInArray + 1];
						gValue3 = (unsigned char)pgImg[indexInArray + imageWidth + 1];
						if (gValue1 >= gValue2 && gValue1 >= gValue3) {//up-right
							x = x + 1;
							y = y - 1;
						}
						else if (gValue3 >= gValue2 && gValue3 >= gValue1) {//down-right
							x = x + 1;
							y = y + 1;
						}
						else {//straight-right
							x = x + 1;
						}
						lastDirection = RightDir;
					}
					else if (lastDirection == LeftDir || shouldGoDirection == LeftDir) {//go left
						if (x == 0 || y == 0 || y == imageHeight - 1) {//reach the image border
							break;
						}
						// Look at 3 neighbors to the left and pick the one with the max. gradient value
						gValue1 = (unsigned char)pgImg[indexInArray - imageWidth - 1];
						gValue2 = (unsigned char)pgImg[indexInArray - 1];
						gValue3 = (unsigned char)pgImg[indexInArray + imageWidth - 1];
						if (gValue1 >= gValue2 && gValue1 >= gValue3) {//up-left
							x = x - 1;
							y = y - 1;
						}
						else if (gValue3 >= gValue2 && gValue3 >= gValue1) {//down-left
							x = x - 1;
							y = y + 1;
						}
						else {//straight-left
							x = x - 1;
						}
						lastDirection = LeftDir;
					}
				}
				else {//should go up or down.
					if (lastDirection == RightDir || lastDirection == LeftDir) {//change the pixel direction now
						if (y>lastY) {//should go down
							shouldGoDirection = DownDir;
						}
						else {//should go up
							shouldGoDirection = UpDir;
						}
					}
					lastX = x;
					lastY = y;
					if (lastDirection == DownDir || shouldGoDirection == DownDir) {//go down
						if (x == 0 || x == imageWidth - 1 || y == imageHeight - 1) {//reach the image border
							break;
						}
						// Look at 3 neighbors to the down and pick the one with the max. gradient value
						gValue1 = (unsigned char)pgImg[indexInArray + imageWidth + 1];
						gValue2 = (unsigned char)pgImg[indexInArray + imageWidth];
						gValue3 = (unsigned char)pgImg[indexInArray + imageWidth - 1];
						if (gValue1 >= gValue2 && gValue1 >= gValue3) {//down-right
							x = x + 1;
							y = y + 1;
						}
						else if (gValue3 >= gValue2 && gValue3 >= gValue1) {//down-left
							x = x - 1;
							y = y + 1;
						}
						else {//straight-down
							y = y + 1;
						}
						lastDirection = DownDir;
					}
					else if (lastDirection == UpDir || shouldGoDirection == UpDir) {//go up
						if (x == 0 || x == imageWidth - 1 || y == 0) {//reach the image border
							break;
						}
						// Look at 3 neighbors to the up and pick the one with the max. gradient value
						gValue1 = (unsigned char)pgImg[indexInArray - imageWidth + 1];
						gValue2 = (unsigned char)pgImg[indexInArray - imageWidth];
						gValue3 = (unsigned char)pgImg[indexInArray - imageWidth - 1];
						if (gValue1 >= gValue2 && gValue1 >= gValue3) {//up-right
							x = x + 1;
							y = y - 1;
						}
						else if (gValue3 >= gValue2 && gValue3 >= gValue1) {//up-left
							x = x - 1;
							y = y - 1;
						}
						else {//straight-up
							y = y - 1;
						}
						lastDirection = UpDir;
					}
				}
				indexInArray = y * imageWidth + x;
			}//end while go right
			 //then go left, pixel direction may be different during linking.
			x = pAnchorX_[i];
			y = pAnchorY_[i];
			indexInArray = y * imageWidth + x;
			pEdgeImg[indexInArray] = 0;//mark the anchor point be a non-edge pixel and
			lastDirection = LeftDir;
			pSecondPartEdgeS_[offsetPS] = offsetPSecond;
			while (pgImg[indexInArray]>0 && !pEdgeImg[indexInArray]) {
				pEdgeImg[indexInArray] = 1;        // Mark this pixel as an edge pixel
				pSecondPartEdgeX_[offsetPSecond] = x;
				pSecondPartEdgeY_[offsetPSecond++] = y;
				shouldGoDirection = 0;//unknown
				if (pdirImg[indexInArray] == Horizontal) {//should go left or right
					if (lastDirection == UpDir || lastDirection == DownDir) {//change the pixel direction now
						if (x>lastX) {//should go right
							shouldGoDirection = RightDir;
						}
						else {//should go left
							shouldGoDirection = LeftDir;
						}
					}
					lastX = x;
					lastY = y;
					if (lastDirection == RightDir || shouldGoDirection == RightDir) {//go right
						if (x == imageWidth - 1 || y == 0 || y == imageHeight - 1) {//reach the image border
							break;
						}
						// Look at 3 neighbors to the right and pick the one with the max. gradient value
						gValue1 = (unsigned char)pgImg[indexInArray - imageWidth + 1];
						gValue2 = (unsigned char)pgImg[indexInArray + 1];
						gValue3 = (unsigned char)pgImg[indexInArray + imageWidth + 1];
						if (gValue1 >= gValue2 && gValue1 >= gValue3) {//up-right
							x = x + 1;
							y = y - 1;
						}
						else if (gValue3 >= gValue2 && gValue3 >= gValue1) {//down-right
							x = x + 1;
							y = y + 1;
						}
						else {//straight-right
							x = x + 1;
						}
						lastDirection = RightDir;
					}
					else	if (lastDirection == LeftDir || shouldGoDirection == LeftDir) {//go left
						if (x == 0 || y == 0 || y == imageHeight - 1) {//reach the image border
							break;
						}
						// Look at 3 neighbors to the left and pick the one with the max. gradient value
						gValue1 = (unsigned char)pgImg[indexInArray - imageWidth - 1];
						gValue2 = (unsigned char)pgImg[indexInArray - 1];
						gValue3 = (unsigned char)pgImg[indexInArray + imageWidth - 1];
						if (gValue1 >= gValue2 && gValue1 >= gValue3) {//up-left
							x = x - 1;
							y = y - 1;
						}
						else if (gValue3 >= gValue2 && gValue3 >= gValue1) {//down-left
							x = x - 1;
							y = y + 1;
						}
						else {//straight-left
							x = x - 1;
						}
						lastDirection = LeftDir;
					}
				}
				else {//should go up or down.
					if (lastDirection == RightDir || lastDirection == LeftDir) {//change the pixel direction now
						if (y>lastY) {//should go down
							shouldGoDirection = DownDir;
						}
						else {//should go up
							shouldGoDirection = UpDir;
						}
					}
					lastX = x;
					lastY = y;
					if (lastDirection == DownDir || shouldGoDirection == DownDir) {//go down
						if (x == 0 || x == imageWidth - 1 || y == imageHeight - 1) {//reach the image border
							break;
						}
						// Look at 3 neighbors to the down and pick the one with the max. gradient value
						gValue1 = (unsigned char)pgImg[indexInArray + imageWidth + 1];
						gValue2 = (unsigned char)pgImg[indexInArray + imageWidth];
						gValue3 = (unsigned char)pgImg[indexInArray + imageWidth - 1];
						if (gValue1 >= gValue2 && gValue1 >= gValue3) {//down-right
							x = x + 1;
							y = y + 1;
						}
						else if (gValue3 >= gValue2 && gValue3 >= gValue1) {//down-left
							x = x - 1;
							y = y + 1;
						}
						else {//straight-down
							y = y + 1;
						}
						lastDirection = DownDir;
					}
					else	if (lastDirection == UpDir || shouldGoDirection == UpDir) {//go up
						if (x == 0 || x == imageWidth - 1 || y == 0) {//reach the image border
							break;
						}
						// Look at 3 neighbors to the up and pick the one with the max. gradient value
						gValue1 = (unsigned char)pgImg[indexInArray - imageWidth + 1];
						gValue2 = (unsigned char)pgImg[indexInArray - imageWidth];
						gValue3 = (unsigned char)pgImg[indexInArray - imageWidth - 1];
						if (gValue1 >= gValue2 && gValue1 >= gValue3) {//up-right
							x = x + 1;
							y = y - 1;
						}
						else if (gValue3 >= gValue2 && gValue3 >= gValue1) {//up-left
							x = x - 1;
							y = y - 1;
						}
						else {//straight-up
							y = y - 1;
						}
						lastDirection = UpDir;
					}
				}
				indexInArray = y * imageWidth + x;
			}//end while go left
			 //end anchor is Horizontal
		}
		else {//the direction of this pixel is vertical, go up and down
			  //fist go down, pixel direction may be different during linking.
			lastDirection = DownDir;
			while (pgImg[indexInArray]>0 && !pEdgeImg[indexInArray]) {
				pEdgeImg[indexInArray] = 1;        // Mark this pixel as an edge pixel
				pFirstPartEdgeX_[offsetPFirst] = x;
				pFirstPartEdgeY_[offsetPFirst++] = y;
				shouldGoDirection = 0;//unknown
				if (pdirImg[indexInArray] == Horizontal) {//should go left or right
					if (lastDirection == UpDir || lastDirection == DownDir) {//change the pixel direction now
						if (x>lastX) {//should go right
							shouldGoDirection = RightDir;
						}
						else {//should go left
							shouldGoDirection = LeftDir;
						}
					}
					lastX = x;
					lastY = y;
					if (lastDirection == RightDir || shouldGoDirection == RightDir) {//go right
						if (x == imageWidth - 1 || y == 0 || y == imageHeight - 1) {//reach the image border
							break;
						}
						// Look at 3 neighbors to the right and pick the one with the max. gradient value
						gValue1 = (unsigned char)pgImg[indexInArray - imageWidth + 1];
						gValue2 = (unsigned char)pgImg[indexInArray + 1];
						gValue3 = (unsigned char)pgImg[indexInArray + imageWidth + 1];
						if (gValue1 >= gValue2 && gValue1 >= gValue3) {//up-right
							x = x + 1;
							y = y - 1;
						}
						else if (gValue3 >= gValue2 && gValue3 >= gValue1) {//down-right
							x = x + 1;
							y = y + 1;
						}
						else {//straight-right
							x = x + 1;
						}
						lastDirection = RightDir;
					}
					else	if (lastDirection == LeftDir || shouldGoDirection == LeftDir) {//go left
						if (x == 0 || y == 0 || y == imageHeight - 1) {//reach the image border
							break;
						}
						// Look at 3 neighbors to the left and pick the one with the max. gradient value
						gValue1 = (unsigned char)pgImg[indexInArray - imageWidth - 1];
						gValue2 = (unsigned char)pgImg[indexInArray - 1];
						gValue3 = (unsigned char)pgImg[indexInArray + imageWidth - 1];
						if (gValue1 >= gValue2 && gValue1 >= gValue3) {//up-left
							x = x - 1;
							y = y - 1;
						}
						else if (gValue3 >= gValue2 && gValue3 >= gValue1) {//down-left
							x = x - 1;
							y = y + 1;
						}
						else {//straight-left
							x = x - 1;
						}
						lastDirection = LeftDir;
					}
				}
				else {//should go up or down.
					if (lastDirection == RightDir || lastDirection == LeftDir) {//change the pixel direction now
						if (y>lastY) {//should go down
							shouldGoDirection = DownDir;
						}
						else {//should go up
							shouldGoDirection = UpDir;
						}
					}
					lastX = x;
					lastY = y;
					if (lastDirection == DownDir || shouldGoDirection == DownDir) {//go down
						if (x == 0 || x == imageWidth - 1 || y == imageHeight - 1) {//reach the image border
							break;
						}
						// Look at 3 neighbors to the down and pick the one with the max. gradient value
						gValue1 = (unsigned char)pgImg[indexInArray + imageWidth + 1];
						gValue2 = (unsigned char)pgImg[indexInArray + imageWidth];
						gValue3 = (unsigned char)pgImg[indexInArray + imageWidth - 1];
						if (gValue1 >= gValue2 && gValue1 >= gValue3) {//down-right
							x = x + 1;
							y = y + 1;
						}
						else if (gValue3 >= gValue2 && gValue3 >= gValue1) {//down-left
							x = x - 1;
							y = y + 1;
						}
						else {//straight-down
							y = y + 1;
						}
						lastDirection = DownDir;
					}
					else	if (lastDirection == UpDir || shouldGoDirection == UpDir) {//go up
						if (x == 0 || x == imageWidth - 1 || y == 0) {//reach the image border
							break;
						}
						// Look at 3 neighbors to the up and pick the one with the max. gradient value
						gValue1 = (unsigned char)pgImg[indexInArray - imageWidth + 1];
						gValue2 = (unsigned char)pgImg[indexInArray - imageWidth];
						gValue3 = (unsigned char)pgImg[indexInArray - imageWidth - 1];
						if (gValue1 >= gValue2 && gValue1 >= gValue3) {//up-right
							x = x + 1;
							y = y - 1;
						}
						else if (gValue3 >= gValue2 && gValue3 >= gValue1) {//up-left
							x = x - 1;
							y = y - 1;
						}
						else {//straight-up
							y = y - 1;
						}
						lastDirection = UpDir;
					}
				}
				indexInArray = y * imageWidth + x;
			}//end while go down
			 //then go up, pixel direction may be different during linking.
			lastDirection = UpDir;
			x = pAnchorX_[i];
			y = pAnchorY_[i];
			indexInArray = y * imageWidth + x;
			pEdgeImg[indexInArray] = 0;//mark the anchor point be a non-edge pixel and
			pSecondPartEdgeS_[offsetPS] = offsetPSecond;
			while (pgImg[indexInArray]>0 && !pEdgeImg[indexInArray]) {
				pEdgeImg[indexInArray] = 1;        // Mark this pixel as an edge pixel
				pSecondPartEdgeX_[offsetPSecond] = x;
				pSecondPartEdgeY_[offsetPSecond++] = y;
				shouldGoDirection = 0;//unknown
				if (pdirImg[indexInArray] == Horizontal) {//should go left or right
					if (lastDirection == UpDir || lastDirection == DownDir) {//change the pixel direction now
						if (x>lastX) {//should go right
							shouldGoDirection = RightDir;
						}
						else {//should go left
							shouldGoDirection = LeftDir;
						}
					}
					lastX = x;
					lastY = y;
					if (lastDirection == RightDir || shouldGoDirection == RightDir) {//go right
						if (x == imageWidth - 1 || y == 0 || y == imageHeight - 1) {//reach the image border
							break;
						}
						// Look at 3 neighbors to the right and pick the one with the max. gradient value
						gValue1 = (unsigned char)pgImg[indexInArray - imageWidth + 1];
						gValue2 = (unsigned char)pgImg[indexInArray + 1];
						gValue3 = (unsigned char)pgImg[indexInArray + imageWidth + 1];
						if (gValue1 >= gValue2 && gValue1 >= gValue3) {//up-right
							x = x + 1;
							y = y - 1;
						}
						else if (gValue3 >= gValue2 && gValue3 >= gValue1) {//down-right
							x = x + 1;
							y = y + 1;
						}
						else {//straight-right
							x = x + 1;
						}
						lastDirection = RightDir;
					}
					else	if (lastDirection == LeftDir || shouldGoDirection == LeftDir) {//go left
						if (x == 0 || y == 0 || y == imageHeight - 1) {//reach the image border
							break;
						}
						// Look at 3 neighbors to the left and pick the one with the max. gradient value
						gValue1 = (unsigned char)pgImg[indexInArray - imageWidth - 1];
						gValue2 = (unsigned char)pgImg[indexInArray - 1];
						gValue3 = (unsigned char)pgImg[indexInArray + imageWidth - 1];
						if (gValue1 >= gValue2 && gValue1 >= gValue3) {//up-left
							x = x - 1;
							y = y - 1;
						}
						else if (gValue3 >= gValue2 && gValue3 >= gValue1) {//down-left
							x = x - 1;
							y = y + 1;
						}
						else {//straight-left
							x = x - 1;
						}
						lastDirection = LeftDir;
					}
				}
				else {//should go up or down.
					if (lastDirection == RightDir || lastDirection == LeftDir) {//change the pixel direction now
						if (y>lastY) {//should go down
							shouldGoDirection = DownDir;
						}
						else {//should go up
							shouldGoDirection = UpDir;
						}
					}
					lastX = x;
					lastY = y;
					if (lastDirection == DownDir || shouldGoDirection == DownDir) {//go down
						if (x == 0 || x == imageWidth - 1 || y == imageHeight - 1) {//reach the image border
							break;
						}
						// Look at 3 neighbors to the down and pick the one with the max. gradient value
						gValue1 = (unsigned char)pgImg[indexInArray + imageWidth + 1];
						gValue2 = (unsigned char)pgImg[indexInArray + imageWidth];
						gValue3 = (unsigned char)pgImg[indexInArray + imageWidth - 1];
						if (gValue1 >= gValue2 && gValue1 >= gValue3) {//down-right
							x = x + 1;
							y = y + 1;
						}
						else if (gValue3 >= gValue2 && gValue3 >= gValue1) {//down-left
							x = x - 1;
							y = y + 1;
						}
						else {//straight-down
							y = y + 1;
						}
						lastDirection = DownDir;
					}
					else	if (lastDirection == UpDir || shouldGoDirection == UpDir) {//go up
						if (x == 0 || x == imageWidth - 1 || y == 0) {//reach the image border
							break;
						}
						// Look at 3 neighbors to the up and pick the one with the max. gradient value
						gValue1 = (unsigned char)pgImg[indexInArray - imageWidth + 1];
						gValue2 = (unsigned char)pgImg[indexInArray - imageWidth];
						gValue3 = (unsigned char)pgImg[indexInArray - imageWidth - 1];
						if (gValue1 >= gValue2 && gValue1 >= gValue3) {//up-right
							x = x + 1;
							y = y - 1;
						}
						else if (gValue3 >= gValue2 && gValue3 >= gValue1) {//up-left
							x = x - 1;
							y = y - 1;
						}
						else {//straight-up
							y = y - 1;
						}
						lastDirection = UpDir;
					}
				}
				indexInArray = y * imageWidth + x;
			}//end while go up
		}//end anchor is Vertical
		 //only keep the edge chains whose length is larger than the minLineLen_;
		edgeLenFirst = offsetPFirst - pFirstPartEdgeS_[offsetPS];
		edgeLenSecond = offsetPSecond - pSecondPartEdgeS_[offsetPS];
		if (edgeLenFirst + edgeLenSecond<minLineLen_ + 1) {//short edge, drop it
			offsetPFirst = pFirstPartEdgeS_[offsetPS];
			offsetPSecond = pSecondPartEdgeS_[offsetPS];
		}
		else {
			offsetPS++;
		}
	}
	//store the last index
	pFirstPartEdgeS_[offsetPS] = offsetPFirst;
	pSecondPartEdgeS_[offsetPS] = offsetPSecond;
	if (offsetPS>maxNumOfEdge) {
		printf("Edge drawing Error: The total number of edges is larger than MaxNumOfEdge, "
			"numofedge = %d, MaxNumOfEdge = %d\n", offsetPS, maxNumOfEdge);
		return 1;
	}
	if (offsetPFirst>edgePixelArraySize || offsetPSecond>edgePixelArraySize) {
		printf("Edge drawing Error: The total number of edge pixels is larger than MaxNumOfEdgePixels, "
			"numofedgePixel1 = &d,  numofedgePixel2 = %d, MaxNumOfEdgePixel = %d\n", offsetPFirst, offsetPSecond, edgePixelArraySize);
		return 1;
	}

	int tempID;
	edgeChains.xCors.resize(offsetPFirst + offsetPSecond);
	edgeChains.yCors.resize(offsetPFirst + offsetPSecond);
	edgeChains.sId.resize(offsetPS + 1);
	unsigned int *pxCors = edgeChains.xCors.data();
	unsigned int *pyCors = edgeChains.yCors.data();
	unsigned int *psId = edgeChains.sId.data();
	offsetPFirst = 0;
	offsetPSecond = 0;
	unsigned int indexInCors = 0;
	unsigned int numOfEdges = 0;
	for (unsigned int edgeId = 0; edgeId<offsetPS; edgeId++) {

		psId[numOfEdges++] = indexInCors;
		indexInArray = pFirstPartEdgeS_[edgeId];
		offsetPFirst = pFirstPartEdgeS_[edgeId + 1];
		for (tempID = offsetPFirst - 1; tempID >= indexInArray; tempID--) {//add first part edge
			pxCors[indexInCors] = pFirstPartEdgeX_[tempID];
			pyCors[indexInCors++] = pFirstPartEdgeY_[tempID];
		}
		indexInArray = pSecondPartEdgeS_[edgeId];
		offsetPSecond = pSecondPartEdgeS_[edgeId + 1];
		for (tempID = indexInArray + 1; tempID<(int)offsetPSecond; tempID++) {//add second part edge
			pxCors[indexInCors] = pSecondPartEdgeX_[tempID];
			pyCors[indexInCors++] = pSecondPartEdgeY_[tempID];
		}
	}
	psId[numOfEdges] = indexInCors;//the end index of the last edge
	edgeChains.numOfEdges = numOfEdges;


	return 0;
}


#endif


//int DrawLines(EdgeChains edges, unsigned int xsize,
//	unsigned int ysize)
//{
//
//	Mat img(ysize, xsize, CV_8UC3, Scalar(0, 0, 0));
//	vector<Point> P;
//	size_t num = edges.xCors.size();
//	P.resize(num);
//	for (int i = 0; i < num; i++) {
//
//		P[i].x = edges.xCors[i];
//		P[i].y = edges.yCors[i];
//		circle(img, P[i], 0, Scalar(255, 255, 255),-1 );
//
//	}
//	imshow("Edge of Image", img);
//
//	waitKey();
//	return 0;
//}



//int main() {
//
//	ImgSet Set;
//	Set.InitImg = GetInitImg("Mickey25.png");
//	Set.FilterImg = GetFilterImg(Set.InitImg);
//
//	uint xsize, ysize;
//	uchar *data;
//	GetFImgPara(xsize, ysize, data, Set);
//
//
//	EdgeChains edges;
//	EDLineDetector* EDL = NULL;
//	EDL = new EDLineDetector[1];
//	
//	image_int8u_p image = NULL;
//	image =	new_image_int8u_ptr(xsize, ysize, data);
//
//	if (EDL->EdgeDrawing(image, edges)) {
//		printf("Line Detection not finished\n");
//		return 1;
//	}
//	else
//		printf("Line Detection is finished\n");
//
//
//	DrawLines(edges, xsize, ysize);
//
//	cout << xsize << "\t" << ysize << endl;
//
//	system("pause");
//}

