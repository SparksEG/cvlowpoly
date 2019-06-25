#ifndef __UTILS__
#define __UTILS__

#include "PathCompacter.h"
#include "edlines.h"
#include "ImgLib.h"
#include "JFA.h"

typedef vector<Pixel> EdgeChain;
typedef vector<vector<Pixel> > VecEdgeChains;

/****************************
���ܣ��� �㷨edlines�����ݽṹ ��EdgeChains�� ת���� ɢ�� �����ݽṹ��EdgeChain��
���룺��ʽvec�洢������ͼ���Ե edges
�������vec�洢������ͼ���Ե
���أ��Ƿ�ɹ�
****************************/
int EdgeChainsToSeeds(EdgeChain &seeds, const EdgeChains edges);

/****************************
���ܣ��� �㷨edlines�����ݽṹ ��EdgeChains�� ת����  vec<vec<Pixel> > �����ݽṹ��VecEdgeChains��
���룺��ʽvec�洢������ͼ���Ե edges
�����˫vec�洢������ͼ���Ե
���أ��Ƿ�ɹ�
****************************/
int EdgeChainsToVecEdgeChains(VecEdgeChains &vedges, const EdgeChains edges);


/****************************
���ܣ��� vec<vec> Pixel �����ݽṹ��VecEdgeChains�� ת����  �㷨edlines�����ݽṹ ��EdgeChains��
���룺˫vec�洢������ͼ���Ե
�������ʽvec�洢������ͼ���Ե edges
���أ��Ƿ�ɹ�
****************************/
int VecEdgeChainsToEdgeChains(EdgeChains &edges, const VecEdgeChains vedges);


/****************************
���ܣ��� vec<Pixel>  �����ݽṹ��EdgeChain��  ת����  �㷨PathCompact�����ݽṹ ��DVector2D*��
���룺vec�洢��һ��ͼ���Ե
���أ�����ָ��洢��һ��ͼ���Ե
****************************/
DVector2D* EDToPath(EdgeChain edge);


/****************************
���ܣ��� �㷨PathCompact�����ݽṹ��DVector2D*�� ת����  vec<Pixel> �����ݽṹ��EdgeChain��
���룺����ָ��洢��һ��ͼ���Եedge���ñ�Ե�����ص����
���أ�vec�洢��һ��ͼ���Ե
****************************/
EdgeChain PathToED(DVector2D* segment, unsigned int uPointsInCurrentPath);


/****************************
���ܣ������㼰�侭�������м��߶��ϵĵ��ҵ�����������ӵ�����������EdgeChains����
���룺����p1��p2
���أ����º��������EdgeChains��
****************************/
void PixelsLinkToEdgeChain(EdgeChain& edge, Pixel p1, Pixel p2);


/****************************
���ܣ��� �㷨edlines�����ݽṹ ��EdgeChains�� ת����  �㷨JFA�����ݽṹ ��EdgeChain��
���룺��ʽvec�洢������ͼ���Ե edges
���أ�vec�洢�� ���б�Ե�ϵ�ɢ��
****************************/
void EdgeChainsToJFASeeds(EdgeChain& Seeds, const EdgeChains& edges);


/****************************
���ܣ������еı�Ե�ϵĵ㻭�ɵ�
���룺����ͼ���Ե edges��ͼ��xsize ��ysize
****************************/
int DrawPlots(EdgeChains edges, unsigned int xsize,
	unsigned int ysize);


/****************************
���ܣ������еı�Ե������
���룺����ͼ���Ե edges��ͼ���xsize ��ysize
****************************/
int DrawLines(EdgeChains edges, unsigned int xsize,
	unsigned int ysize);


#endif
