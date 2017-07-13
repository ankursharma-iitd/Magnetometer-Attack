/**
 * @author Panagiotakis Konstantinos
 * @version     v1.0
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "alloc.h"

#define MAX_ERROR_THRESHOLD 590.0
#define Pi 3.1415927
#define INF 100000000.0
#define ABS(a) ((a) > 0 ? (a) : -(a))
#define SQUARE(a) ((a) * (a))
#define MAX(a,b) ((a) > b ? (a) : (b))
#define MIN(a,b) ((a) < b ? (a) : (b))
#define SIGN(a) ((a) >= 0 ? (1) : -(1))
#define DISTANCE_E(x0,y0,x1,y1) (sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)))
#define CROSS2D(x0,y0,x1,y1) (x0*x1+y0*y1)
#define CROSS3D(x0,y0,z0,x1,y1,z1) (x0*x1+y0*y1+z0*z1)
#define DISTANCE_E2(x0,y0,x1,y1) ((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1))
#define DISTANCE_E4(x0,y0,x1,y1) (((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1))*((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)))
#define DISTANCE_E2_3D(x0,y0,z0,x1,y1,z1) ((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1))
#define DISTANCE_E_3D(x0,y0,z0,x1,y1,z1) (sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1)))


#define SWAP(a,b,t) {t = a; a = b; b = t;}

#define BACKGROUND 0
#define OBJECT 255

//DEFINES GIA to EP

#define MAXIMUMDEGREE 10

#define L1_METRIC 1 //MANHATTAN
#define L2_METRIC 2 //EUCLIDEAN
#define LINF_METRIC 3
#define POLYGONERROR_LISE_METRIC 4
#define POLYGONERROR_ToleranceZone_METRIC 5
#define LABS_METRIC 6
#define LKEYFRAMES_METRIC 7
#define MIN_MAX_KEYFRAMES_METRIC 8
#define L2KEYFRAMES_METRIC 9

#define ADD_TIME_TO_METRIC 0 //EUCLIDEAN :d(x,y)+ |x-y|
#define ADD_EXP_TO_METRIC 1 //EUCLIDEAN :d(x,y)+ |x-y|*exp(-|x-y|)


#define TRUE    1
#define FALSE   0


struct point2DF 
{
	double x;
	double y;
	double d; //value of d(x,y)
	int fatherLine;
	int degree; //degree of point 
	int line[MAXIMUMDEGREE];//8ewroume pws o megistos ari8mos degree einai MAXIMUMDEGREE
};

struct lineSegment
{
	int apo; //pointers to points  //LS: [apo,eos)
	int eos;
	int legal;
	int numChilds;
};


struct Solution
{	
	double **p; //numSol * (N - 1)
	int numSol;
	double *r;//average radious per solution  //numSol
};


struct curveEP
{
	double **D; //lines X lines (lines = M)
	double **C; // lines x dim
	int **MinMaxKF;//lines X lines (for the MINMAX algo for key frames)

	struct point2DF **p;
	struct lineSegment **l;
    int *NumberOfPoints;
	int *NumberOfLines;

	int ***linesId; //[level][x][0..histLines[level][x]] , x \in [0,M-1]
	int **numLines; //[level][x] ===> # eu8eiwn pou brikontai sto x

	struct Solution sol; //lush tou problhmatos

	int N; //# of EP points  
	int M; //# array points
	int n; //curve dimension 

	int level; //level of  Null Plane Curve

	double Threshold;//Threshold for keyframes

	char outDir[500];
	char outName[500];
};

struct Box{
	int xmin,ymin,xmax,ymax;
}; 


struct Sunoro
{
	struct Point *Head;
	struct Point *Tail;
	int size;
};


struct Point
{
        int x,y;
        struct Point *next,*prev;
};


int lines,cols,xmin,xmax,ymin,ymax; 
double gThresh;
struct curveEP *ep;
int metricChoise;

/*****************/
/* MinMaxEP		*/
/***************/
extern void getDistanceMatrixFromCurveUsingMinMax(char *name, int choice);
extern void printfMinMaxSolutions();

/*****************/
/* appr-eqp.c	*/
/***************/
extern void runApprEQP();

/*******
* region.c
********
*/

extern double getL2Norm(double *P0,double *P1,int dim);
extern double getDistancePointLineSegment(double *P,double *P0, double *P1,int dim);
extern double sufraceOfTrig(int x0,int y0,int x1,int y1,int x2,int y2);
extern void cross(double *r,double *x1,double *x2);
extern double getDistFromLineSegment(double Mx,double My,double Nx,double Ny,double x,double y);
extern void normalize(double *r,int len);
extern void getRotationMatrix(double DOFx,double DOFy,double DOFz,double **R);
extern void Midpoint(int x0, int y0, int x1, int y1, int width, int ***image,int r0,int g0,int b0); 
extern int Midpoint2D(int x0, int y0, int x1, int y1, int width, int **image,int r0);
extern struct Sunoro *newSunoro();
extern void deletesunoro(struct Sunoro *A);
extern void addpoint(struct Point *new,struct Sunoro *A);
extern void deletepoint(struct Point *del,struct Sunoro *A);
extern int findRegion2(int **X,double v,double nv,int *xc,int *yc,int lines,int cols);
extern void findRegion(int **X,double v,double nv);
extern void findRegionForHumanTracking(int **X,double v,double nv,int xc,int yc,int x0,int y0,int x1,int y1,int lines,int cols,int height);
extern void changeRegion(int **X,double v,double nv,int x0,int y0,int Thre);
extern void getSufraceOfArea(int **X,int **V,double v,int x0,int y0,int label,int *Res);
extern void getSufraceOfArea1(int **X,int **V,double v,int x0,int y0,int label,int *Res,int Thres);
extern int getSufraceOfArea3(int **X,int **V,int x0,int y0,int thres);
extern void getSufraceOfArea_CenterMass(int **X,int **V,double v,int x0,int y0,int label,int *Res,int *xc,int *yc);
extern int getLastPointOfArea(int **X,int **V, int x0,int y0,int v,int *xp,int *yp);
extern double getNearestPointFormArea(int **X,int x0,int y0,int V,int *xp,int *yp);
extern void Ison(int **A,int xmin,int xmax,int ymin,int ymax,int **C);
extern int getLastPointOfArea(int **X,int **V, int x0,int y0,int v,int *xp,int *yp);
extern double getNearestPointFormArea(int **X,int x0,int y0,int V,int *xp,int *yp);
extern void RotatePointF(double i,double j,double fi,double xc,double yc,double *xn,double *yn);
extern void pollpin(double **A,int lines,int cols1,double **B,int lines2,int cols,double **C);
extern void expansionFilter(int ***In,int ***Out,int **BitIn,int Back,int **BitOut,int lines,int cols);
extern void Diabrwsh(int **In,int **Out,int lines,int cols,int r);
extern void quicksort1D( double *v, int left, int right );
extern void print_box2(int **A,int x1,int y1,int x2,int y2,int r0);


