/**
 * @author Panagiotakis Konstantinos
 * @version     v1.0
 */


/*Functions for input - output and init*/

#include "EP.h"


void runMidpoint2D(int x0, int y0, int x1, int y1, int width, int **image,int r0,int r1,int *p0,int *p1);



int FRAMEV;
int **Gt;
double **TempF;
int **Temp1;

//f = x.y (eswteriko ginomeno)
double getProduct(double *x,double *y,int dim)
{
	int i;
	double f = 0;

	for (i = 0; i < dim; ++i)
	{
		f += x[i]*y[i];		
	}

	return f;
}



//f = |x-y|_2 L2 norm
double getL2Distance(double *x,double *y,int dim)
{
	int i;
	double f = 0,w;

	
	for (i = 0; i < dim; ++i)
	{
		w = x[i] - y[i];
		f += w*w;
	}

	f = (double)(sqrt(f));

	return f;
}


/*
distance( Point P, Segment P0:P1 )
{
      v = P1 - P0
      w = P - P0
      if ( (c1 = w·v) <= 0 )
            return d(P, P0)
      if ( (c2 = v·v) <= c1 )
            return d(P, P1)
      b = c1 / c2
      Pb = P0 + bv
      return d(P, Pb)
}
*/

//Distance between the point ep->C[m] and E. t  ep->C[n1]ep->C[n2]
double getDistancePointLineSegment(double *P,double *P0, double *P1,int dim)
{
	double *v = (double *)malloc(dim*sizeof(double));
	double *w = (double *)malloc(dim*sizeof(double));
	double *Pb = (double *)malloc(dim*sizeof(double));
	double b,c1,c2;
	int i;

	for (i = 0; i < dim; ++i)
	{
      v[i] = P1[i] - P0[i];
      w[i] = P[i] - P0[i];
	}

	
      if ( (c1 = getProduct(w,v,dim)) <= 0 )
	  {
		  free(v);
		  free(w);	
		  free(Pb);

		  return getL2Distance(P,P0,dim);		//d(P, P0);
	  }
      if ( (c2 = getProduct(v,v,dim)) <= c1 )
	  {
		  free(v);
		  free(w);	
		  free(Pb);

          return getL2Distance(P,P1,dim);		//d(P, P1);
	  }

      b = c1 / c2;

	  for (i = 0; i < dim; ++i)
	  {
		Pb[i] = P0[i] +  b*v[i];
	  }

	 

     // printf("c1 = %f c2 = %f D = %f",c1,c2,getL2Distance(P, Pb, dim));

	  b = getL2Distance(P, Pb, dim);
	
	  free(v);
	  free(w);	
	  free(Pb);

      return b;
}
//Distance between the points ep->C[m] - ep->C[n]
double getL2Norm(double *P0,double *P1,int dim)
{
	double f = 0,w;
	int i;

	for (i = 0; i < dim; ++i)
	{
		w = (P0[i] - P1[i]);
		w = w*w;
		f += w;
	}

	f = (double)sqrt(f);

	return f;
}




//Several functions 
//******************************************************************************************


void cross(double *r,double *x1,double *x2)
{
	r[0] = x1[1]*x2[2] - x1[2]*x2[1];
	r[1] = - x1[0]*x2[2] + x1[2]*x2[0];
	r[2] = x1[0]*x2[1] - x1[1]*x2[0];
}


//Area of triangle 012
double sufraceOfTrig(int x0,int y0,int x1,int y1,int x2,int y2)
{
	double suf = 0;

	suf = (double)(0.5*(x0*(y1-y2) - y0*(x1-x2) + x1*y2 - x2*y1));

	if (suf < 0)
		suf = - suf;

	return suf;
}

void normalize(double *r,int len)
{
	double metro = 0;
	int i;
	

	for (i = 0; i < len; ++i)
		metro += r[i]*r[i];

	if (metro != 0)
	{	
		metro = (double)(  sqrt(metro));

		for (i = 0; i < len; ++i)
		{
			r[i] = r[i] / metro;
		}
	}
}

/*	Distance between the points  <x,y> and the segment MN  */

double getDistFromLineSegment(double Mx,double My,double Nx,double Ny,double x,double y)
{
	
	double l,d3,A,B,C,dxM,dyM,dxN,dyN;

	dxM = x - Mx;
	dyM = y - My;

	dxN = Nx - Mx;
	dyN = Ny - My;
	 

    l = (dxM*dxN + dyM*dyN ) / DISTANCE_E2(Nx,Ny,Mx,My);

    if (l < 0)
		d3 = DISTANCE_E(Mx,My,x,y);
    else if (l > 1)
		d3 = DISTANCE_E(Nx,Ny,x,y);
    else
	{
		A = Ny - My;
		B = -( Nx - Mx);
		C = -Mx*(Ny-My) + My*(Nx-Mx);

        d3 = ABS(A*x+B*y+C) / sqrt(A*A+B*B);
    }
	
	return d3;
}

void Diabrwsh(int **In,int **Out,int lines,int cols,int r)
{
	int i,j,i1,i2;

	for (i = 0; i < lines; i++)
	for (j = 0; j < cols; j++)
	{
	
		Out[i][j] = In[i][j]; 
	}


	for (i = r; i < lines-r; i++)
	for (j = r; j < cols-r; j++)
	{
	
		if (In[i][j] == OBJECT)
		{
			for (i1 = -r; i1 <= r; ++i1)
			for (i2 = -r; i2 <= r; ++i2)
			{
				if (In[i+i1][j+i2] == BACKGROUND)
				{
					Out[i][j] = BACKGROUND;
				}
			}
		}
	}


}

/*
	Enlarges by one pixel the object that have background Back
*/


void expansionFilter(int ***In,int ***Out,int **BitIn,int Back,int **BitOut,int lines,int cols)
{
	int i,ib,jb,j,k,i1,i2,w0,w1,w2,ok;

	for (i = 0; i < lines; i++)
	for (j = 0; j < cols; j++)
	{
		for (k = 0; k < 3; ++k)
		{
			Out[i][j][k] = In[i][j][k];
		}
		BitOut[i][j] = BitIn[i][j]; 
	}

	for (i = 2; i < lines-2; i++)
	for (j = 2; j < cols-2; j++)
	{
		w0 = In[i][j][0];
		w1 = In[i][j][1];
		w2 = In[i][j][2];

		if (BitIn[i][j] == Back)
		{
			ok = 0;
			for (i1 = -1; i1 <= 1; ++i1)
			for (i2 = -1; i2 <= 1; ++i2)
			{
				if (BitIn[i1+i][i2+j] != Back)
				{
					ok = 1;
					ib = i1+i;
					jb = i2+j;
					break;
				}
			}

			if (ok == 1)
			{
				BitOut[i][j] = BitIn[ib][jb];
				for (k = 0; k < 3; ++k)
				{
					Out[i][j][k] = In[ib][jb][k];
				}
				
			}

		}
		
	}
}



double makeInterpolation(double x1,double x2,double pos)
{

	double r = x1 + (pos*(x2-x1));

	return r;
}



double getBlockMatching(int ***LAB,int ***Window,int size,int x0,int y0,int dx,int dy,int *bestX,int *bestY)
{
	int i,j,i1,i2;
	double sf,minsf = 100000000;	
	int x,y;

	for (i = -dx; i <= dx; ++i)
	for (j = -dy; j <= dy; ++j)
	{
		sf = 0;
		x = x0+i;
		y = y0+i;
		for (i1 = -size; i1 <= size; ++i1)
		for (i2 = -size; i2 <= size; ++i2)
		{
			sf += ABS(LAB[x+i1][y+i2][0] - Window[size+i1][size+i2][0]) +  ABS(LAB[x+i1][y+i2][1] - Window[i1+size][i2+size][1]) +  ABS(LAB[x+i1][y+i2][2] - Window[i1+size][i2+size][2]);
			if (sf > minsf)
			{
				i1 = i2 = 2*size;
			}
		}
		if (sf < minsf)
		{
			minsf = sf;
			*bestX = x;
			*bestY = y;
		}
		
	}
	return minsf;

}



void getSufraceOfArea1(int **X,int **V,double v,int x0,int y0,int label,int *Res,int Thres)
{
	int i ,j ,ok = 0;
	struct Sunoro *A = newSunoro();
	struct Point *new;


	while (1)
	{
		if (ok != 0)
		{
			x0 = A->Tail->x;
			y0 = A->Tail->y;
			new = A->Tail;
			deletepoint(new,A);
			new = 0;
		}
		
		ok = 1;
		for (i = x0-1; i <= x0+1; ++i) 
		for (j = y0-1; j <= y0+1; ++j) 
		{
			if (i >= 0 && i < lines && j >= 0 && j < cols && X[i][j] == v && V[i][j] == 0)
			{
				V[i][j] = label;
				*Res += 1;
				new = (struct Point *)malloc(sizeof(struct Point ));
				new->x = i;
				new->y = j;
				addpoint(new,A);

				if (*Res > Thres)
				{
					deletesunoro(A);
					return;
				}
			}
		}
		if (A->Tail == 0 || A->size == 0)
			break;
	}
	deletesunoro(A);
}


int getLastPointOfArea(int **X,int **V, int x0,int y0,int v,int *xp,int *yp)
{
	int i ,j ,ok = 0;
	struct Sunoro *A = newSunoro();
	struct Point *new;
	int res = 1;

	V[x0][y0] = v*v+1;

	*xp = x0;
	*yp = y0;

	while (1)
	{
		if (ok != 0)
		{
			x0 = A->Head->x;
			y0 = A->Head->y;
			new = A->Head;
			deletepoint(new,A);
			new = 0;
		}
		
		ok = 1;
		for (i = x0-1; i <= x0+1; ++i) 
		for (j = y0-1; j <= y0+1; ++j) 
		{
			if (i >= 0 && i < lines && j >= 0 && j < cols && X[i][j] == v && V[i][j] == v)
			{
				res++;
				new = (struct Point *)malloc(sizeof(struct Point ));
				new->x = i;
				new->y = j;
				addpoint(new,A);
				V[i][j] = v*v+1;
				*xp = i;
				*yp = j;

			}
		}
		if (A->Tail == 0 || A->size == 0)
			break;
		
		
	}
	deletesunoro(A);
	return res;
}


int getSufraceOfArea3(int **X,int **V, int x0,int y0,int thres)
{
	int i ,j ,ok = 0;
	struct Sunoro *A = newSunoro();
	struct Point *new;
	int v = X[x0][y0];
	int res = 0;

	V[x0][y0] = 1;

	while (1)
	{
		if (ok != 0)
		{
			x0 = A->Head->x;
			y0 = A->Head->y;
			new = A->Head;
			deletepoint(new,A);
			new = 0;
		}
		
		ok = 1;
		for (i = x0-1; i <= x0+1; ++i) 
		for (j = y0-1; j <= y0+1; ++j) 
		{
			if (i >= 0 && i < lines && j >= 0 && j < cols && X[i][j] == v && V[i][j] == 0)
			{
				res++;
				new = (struct Point *)malloc(sizeof(struct Point ));
				new->x = i;
				new->y = j;
				addpoint(new,A);
				V[i][j] = 1;
			}
		}
		if (A->Tail == 0 || A->size == 0)
			break;
		if (res >= thres-1)
		{

			return 1;
		}
	}
	deletesunoro(A);
	return 0;
}




void getSufraceOfArea_CenterMass(int **X,int **V,double v,int x0,int y0,int label,int *Res,int *xc,int *yc)
{
	int i ,j ,ok = 0;
	struct Sunoro *A = newSunoro();
	struct Point *new;

	*xc = 0;
	*yc = 0;
	while (1)
	{
		if (ok != 0)
		{
			x0 = A->Tail->x;
			y0 = A->Tail->y;
			new = A->Tail;
			deletepoint(new,A);
			new = 0;
		}
		
		ok = 1;
		for (i = x0-1; i <= x0+1; ++i) 
		for (j = y0-1; j <= y0+1; ++j) 
		{
			if (i >= 0 && i < lines && j >= 0 && j < cols && X[i][j] == v && V[i][j] == 0)
			{
				V[i][j] = label;
				*Res += 1;
				*xc += i;
				*yc += j;
				new = (struct Point *)malloc(sizeof(struct Point ));
				new->x = i;
				new->y = j;
				addpoint(new,A);
			}
		}
		if (A->Tail == 0 || A->size == 0)
			break;
	}
	deletesunoro(A);
	*xc = (*xc) / (*Res);
	*yc = (*yc) / (*Res);
}



void getSufraceOfArea(int **X,int **V,double v,int x0,int y0,int label,int *Res)
{
	int i ,j ,ok = 0;
	struct Sunoro *A = newSunoro();
	struct Point *new;


	while (1)
	{
		if (ok != 0)
		{
			x0 = A->Tail->x;
			y0 = A->Tail->y;
			new = A->Tail;
			deletepoint(new,A);
			new = 0;
		}
		
		ok = 1;
		for (i = x0-1; i <= x0+1; ++i) 
		for (j = y0-1; j <= y0+1; ++j) 
		{
			if (i >= 0 && i < lines && j >= 0 && j < cols && X[i][j] == v && V[i][j] == 0)
			{
				V[i][j] = label;
				*Res += 1;
				new = (struct Point *)malloc(sizeof(struct Point ));
				new->x = i;
				new->y = j;
				addpoint(new,A);
			}
		}
		if (A->Tail == 0 || A->size == 0)
			break;
	}
	deletesunoro(A);
}

void changeRegion(int **X,double v,double nv,int x0,int y0,int Thre)
{
	int i,j;
	int **V =(int **) Alloc2D(lines,sizeof(int*),cols,sizeof(int));
	int *D = (int *)calloc(100,sizeof(int));

	getSufraceOfArea1(X,V,v,x0,y0,1,&(D[0]),Thre);
 
	if (D[0] >= Thre)
	{
		free(D);
		FreeInt(V,lines);
		return;
	}
			

	for (i = xmin; i < xmax; i++)
	for (j = ymin; j < ymax; j++)
	{
		if (X[i][j] == v && V[i][j] == 1)
		{
			X[i][j] = (int) nv;
		}
	}

	free(D);
	FreeInt(V,lines);
}



int findRegion2(int **X,double v,double nv,int *xc,int *yc,int lines,int cols)
{
	int i,j;
	int **V =(int **) Alloc2D(lines,sizeof(int*),cols,sizeof(int));
	int *D = (int *)calloc(1000,sizeof(int));
	int maxL = 1,c;
	int maxV = 0,maxI = 0;
	*xc = 0;
	*yc = 0;

	for (i = xmin; i < xmax; i++)
	for (j = ymin; j < ymax; j++)
	{
		if (X[i][j] == v && V[i][j] == 0)
		{
			getSufraceOfArea(X,V,v,i,j,maxL,&(D[maxL]));
			maxL++;
		}
	}

	
	//printf("MAXL = %d \n",maxL);
	for (i = 1; i < maxL; ++i)
	{
	//	printf("E[%d] = %d\n",i,D[i]);
		if (D[i] > maxV)
		{
			maxV = D[i];
			maxI = i;
		}
	}
	c= 0;

	for (i = xmin; i < xmax; i++)
	for (j = ymin; j < ymax; j++)
	{
		if (X[i][j] == v && V[i][j] == maxI)
		{
			X[i][j] = (int)nv;
			*xc += i;
			*yc += j;
			c++;

		}
	}

	if (c == 0)
	{
		printf("ERROR detection\n");
		FreeInt(V,lines);
		free(D);
		return 0;
	}
	*xc = (int)((double)((*xc)/(double)c)+0.5);
	*yc = (int)((double)((*yc)/(double)c)+0.5);

	FreeInt(V,lines);
	free(D);

	return maxV;
}



void findRegionForHumanTracking(int **X,double v,double nv,int xc,int yc,int x0,int y0,int x1,int y1,int lines,int cols,int height)
{
	int i,j,i1,i2;
	int **V =(int **) Alloc2D(lines,sizeof(int*),cols,sizeof(int));
	int *D = (int *)calloc(100,sizeof(int));
	int *Dx = (int *)calloc(100,sizeof(int));
	int *Dy = (int *)calloc(100,sizeof(int));
	int maxL = 1,sizeD = 100;
	int w,d1,hu,bc,d2,d3,ok = 0,maxID,Dmax = 0;
	double c1,c2,c3,max,c0;

	xmin = 2;
	ymin = 2;
	xmax = lines-2;
	ymax = cols-2;
	printf("I1 = %d \n",X[80][193]);

	for (i = xmin; i < xmax; i++)
	for (j = ymin; j < ymax; j++)
	{
		if (X[i][j] == v && V[i][j] == 0)
		{
			getSufraceOfArea_CenterMass(X,V,v,i,j,maxL,&(D[maxL]),&(Dx[maxL]),&(Dy[maxL]));
			if (D[maxL] > Dmax)
			{
			   	Dmax = 	D[maxL];	
			}
			maxL++;
			if (maxL == sizeD-2)
			{
				sizeD = 2*sizeD;
				D = (int *)realloc(D,sizeD*sizeof(int));
				Dx = (int *)realloc(Dx,sizeD*sizeof(int));
				Dy = (int *)realloc(Dy,sizeD*sizeof(int));
				for (i1 = maxL; i1 < sizeD; ++i1)
				{
					D[i1] = 0;
					Dx[i1] = 0;
					Dy[i1] = 0;
				}
			}
		}
	}

	printf("I2 = %d \n",X[80][193]);
	
	w = (height/2);	
	for (i = 1; i < maxL; ++i)
	{
		
		d1 = (int)(DISTANCE_E(xc,yc,Dx[i],Dy[i]));
		d2 = (int)(DISTANCE_E(x0,y0,Dx[i],Dy[i]));
		d3 = (int)(DISTANCE_E(x1,y1,Dx[i],Dy[i]));
		printf("I4 = %d \n",X[80][193]);
		if (xc < 4 || yc < 4 || Dx[i] < 4 ||  Dy[i] < 4 || xc > lines-5 || yc > cols-5 || Dx[i] > lines-5 ||  Dy[i] > cols-5)
		{
			bc = 10000;
			hu = 1;
		}
		else
		{
			runMidpoint2D(xc,yc, Dx[i], Dy[i],3, X,OBJECT,BACKGROUND,&hu,&bc);
		}

		printf("I5 = %d \n",X[80][193]);

		c1 = (0.4*w/(d1+0.0001))*((hu+1.0)/(bc+0.0001));

		c0 = (int)(DISTANCE_E(Dx[i],Dy[i],lines/2.0,cols/2.0));

		c2 = ((0.1*w+2*sqrt(D[i]/Pi))/(0.5*d1+d2+0.0001))  *   ((hu+1.0)/(bc+0.0001));

		//c2 = ((0.25*w)/(d2+0.0001))  *   (3*(hu+1.0)/(bc+0.0001));

		c3 = ((0.1*w+2*sqrt(D[i]/Pi))/(0.5*d1+d3+0.0001))  *   ((hu+1.0)/(bc+0.0001));

		//c3 = ((0.25*w)/(d3+0.0001))  *   (3*(hu+1.0)/(bc+0.0001));
			
		//if ((d1 < 0.5*w || (d2 < 0.2*w+2*sqrt(D[i]/Pi) && d2 < w && hu > 0.2*bc) || (d3 < 0.2*w+2*sqrt(D[i]/Pi) && d3 < w&& hu > 0.2*bc)) || (d1 < 1.5*w && hu > 0.5*bc))
		printf("i = %d c01 = %f D = %d Dx =  %d  %d\n",i,c0,D[i],Dx[i],Dy[i]);
		//getchar();

		if (D[i] < 10 && D[i] < Dmax)
		{
			ok = 1;
			for (i1 = xmin; i1 < xmax; i1++)
			for (i2 = ymin; i2 < ymax; i2++)
				if (X[i1][i2] == v && V[i1][i2] == i)
				{
					X[i1][i2] = (int)nv;
				}
		}
		else if (c0 < 0.15*(lines+cols) && D[i] > 20)
		{

		}
		else if (D[i] > 20 && (c1 > 1 || c2 > 1 || c3 > 1) && (d1 < 1.4*w || d2 < 0.3*w || d3 < 0.3*w))			
		{

		}
		else if (D[i] < Dmax)
		{
			ok = 1;
			for (i1 = xmin; i1 < xmax; i1++)
			for (i2 = ymin; i2 < ymax; i2++)
				if (X[i1][i2] == v && V[i1][i2] == i)
				{
					X[i1][i2] = (int)nv;
				}
		}
	}
	printf("ok = %d\n",ok);
	printf("I3 = %d \n",X[80][193]);
	//getchar();
	if (ok == 0)
	{
		max = -2;
		for (i = 1; i < maxL; ++i)
		{
			if (D[i] > max)
			{
				max = D[i];
				maxID = i;
			}
		}
		i = maxID;

		for (i1 = xmin; i1 < xmax; i1++)
		for (i2 = ymin; i2 < ymax; i2++)
			if (X[i1][i2] == v && V[i1][i2] == i)
			{
				X[i1][i2] = (int)v;
			}

	}



	FreeInt(V,lines);
	free(D);
	free(Dx);
	free(Dy);
}






int computeRegionsForHumanTracking(int *maxL,double *Cost,int **V,int **X,double v,double nv,int xc,int yc,int x0,int y0,int x1,int y1,int lines,int cols,int height,int timeStep,int suf)
{
	int i,j,i1;
	int *D;
	int *Dx;
	int *Dy;
	int sizeD = 100;
	int w,d1,hu,bc,d2,d3;
	double c1,c2,c3;

	printf("Starting REGION %d\n",timeStep);

	if (timeStep == -1)
	{
		FreeInt(V,lines);
		free(Cost);
		return (*maxL);
	}

	
	xmin = 2;
	ymin = 2;
	xmax = lines-2;
	ymax = cols-2;
	

	if (timeStep == 0)
	{
		D = (int *)calloc(100,sizeof(int));
		Dx = (int *)calloc(100,sizeof(int));
		Dy = (int *)calloc(100,sizeof(int));
		*maxL = 1;

		for (i = xmin; i < xmax; i++)
		for (j = ymin; j < ymax; j++)
		{
			V[i][j] = 0;
		}
	
		for (i = xmin; i < xmax; i++)
		for (j = ymin; j < ymax; j++)
		{

			if (X[i][j] == v && V[i][j] == 0)
			{
				getSufraceOfArea_CenterMass(X,V,v,i,j,*maxL,&(D[*maxL]),&(Dx[*maxL]),&(Dy[*maxL]));
				(*maxL)++;
				if (*maxL == sizeD-2)
				{
					sizeD = 2*sizeD;
					D = (int *)realloc(D,sizeD*sizeof(int));
					Dx = (int *)realloc(Dx,sizeD*sizeof(int));
					Dy = (int *)realloc(Dy,sizeD*sizeof(int));
					for (i1 = *maxL; i1 < sizeD; ++i1)
					{
						D[i1] = 0;
						Dx[i1] = 0;
						Dy[i1] = 0;
					}
				}
			}
		}
		

		w = (height/2);	
		for (i = 1; i < (*maxL); ++i)
		{
		
			d1 = (int)(DISTANCE_E(xc,yc,Dx[i],Dy[i]));
			d2 = (int)(DISTANCE_E(x0,y0,Dx[i],Dy[i]));
			d3 = (int)(DISTANCE_E(x1,y1,Dx[i],Dy[i]));

			if (xc < 4 || yc < 4 || Dx[i] < 4 ||  Dy[i] < 4 || xc > lines-5 || yc > cols-5 || Dx[i] > lines-5 ||  Dy[i] > cols-5)
			{	
				bc = 10000;
				hu = 1;
			}
			else
			{
				runMidpoint2D(xc,yc, Dx[i], Dy[i],3, X,OBJECT,BACKGROUND,&hu,&bc);
			}

			/*c1 = (0.5*w/(d1+0.1))*((hu+1.0)/(bc+0.1));

			c2 = ((0.3*w)/(0.5*d1+d2+0.1))  *   ((hu+1.0)/(bc+0.1));

			c3 = ((0.3*w)/(0.5*d1+d3+0.1))  *   ((hu+1.0)/(bc+0.1));*/

			c1 = (0.4*w/(d1+0.0001))*(2.5*(hu+1.0)/(bc+0.0001));
		
			c2 = ((0.3*w+2*sqrt(D[i]/Pi))/(0.5*d1+d2+0.0001))  *   (2.5*(hu+1.0)/(bc+0.0001));
	
			c3 = ((0.3*w+2*sqrt(D[i]/Pi))/(0.5*d1+d3+0.0001))  *   (2.5*(hu+1.0)/(bc+0.0001));



			if ((c1 > 1 || c2 > 1 || c3 > 1) && (d1 < 1.4*w || d2 < 0.3*w || d3 < 0.3*w))			
			{
				Cost[i] = c1 + c2 + c3;	
			}
			else
			{
				Cost[i] = (c1 + c2 + c3) / 10.0;	
			}

//			Cost[i] = MAX(c1,MAX(c2,c3));

//			Cost[i] = c1 + c2 + c3;

			//Cost[i] = c1+0.5*(c2+c3);
			//if (D[i] < 0.02*suf)
			//	Cost[i] = Cost[i] * (D[i]/(0.02*suf));

			printf("Cost[%d] = %f\n",i,Cost[i]);
		
		}
		free(D);
		free(Dx);
		free(Dy);
	}
	return (*maxL);
			
}


double getMaxOf4(double c1,double c2,double c3,double c4)
{

	double max = c1;

	if (c2 > max)
		max = c2;

	if (c3 > max)
		max = c3;

	if (c4 > max)
		max = c4;

	return max;

}

int computeRegionsForHumanTracking4Points(int *maxL,double *Cost,int **V,int **X,double v,double nv,int xc,int yc,int x0,int y0,int x1,int y1,int x2,int y2,int lines,int cols,int height,int timeStep,int suf)
{
	int i,j,i1;
	int *D;
	int *Dx;
	int *Dy;
	int sizeD = 100;
	int w,d1,d2,d3,d4,huC,bcC,huH,bcH,huL,bcL,huLR,bcLR;
	double d,dist;

	printf("Starting REGION %d\n",timeStep);

	if (timeStep == -1)
	{
		FreeInt(V,lines);
		free(Cost);
		return (*maxL);
	}

	
	xmin = 2;
	ymin = 2;
	xmax = lines-2;
	ymax = cols-2;
	

	if (timeStep == 0)
	{
		D = (int *)calloc(100,sizeof(int));
		Dx = (int *)calloc(100,sizeof(int));
		Dy = (int *)calloc(100,sizeof(int));
		*maxL = 1;

		for (i = xmin; i < xmax; i++)
		for (j = ymin; j < ymax; j++)
		{
			V[i][j] = 0;
		}
	
		for (i = xmin; i < xmax; i++)
		for (j = ymin; j < ymax; j++)
		{

			if (X[i][j] == v && V[i][j] == 0)
			{
				getSufraceOfArea_CenterMass(X,V,v,i,j,*maxL,&(D[*maxL]),&(Dx[*maxL]),&(Dy[*maxL]));
				(*maxL)++;
				if (*maxL == sizeD-2)
				{
					sizeD = 2*sizeD;
					D = (int *)realloc(D,sizeD*sizeof(int));
					Dx = (int *)realloc(Dx,sizeD*sizeof(int));
					Dy = (int *)realloc(Dy,sizeD*sizeof(int));
					for (i1 = *maxL; i1 < sizeD; ++i1)
					{
						D[i1] = 0;
						Dx[i1] = 0;
						Dy[i1] = 0;
					}
				}
			}
		}
		
		w = (height/2);	
		for (i = 1; i < (*maxL); ++i)
		{
		
			d1 = (int)(DISTANCE_E(xc,yc,Dx[i],Dy[i]));
			d2 = (int)(getDistFromLineSegment(xc,yc,x0,y0,Dx[i],Dy[i]));
			d3 = (int)(getDistFromLineSegment(xc,yc,x1,y1,Dx[i],Dy[i]));
			d4 = (int)(getDistFromLineSegment(xc,yc,x2,y2,Dx[i],Dy[i]));

			d = dist = -getMaxOf4(-d1,-d2,-d3,-d4);
			
			//Cost[i] = kostos gia na pame sto object i, oso pio megalo toso pio duskolo na pame
			Cost[i] = MAX(0,dist - 0.25*sqrt(D[i]));

			if (d > w || xc < 4 || yc < 4 || Dx[i] < 4 ||  Dy[i] < 4 || xc > lines-5 || yc > cols-5 || Dx[i] > lines-5 ||  Dy[i] > cols-5)
			{	
				Cost[i] = 1000000; //aporiptete logw terastiou la8ous se apostash h embadon
			}
			else if (Cost[i] < 0.05*w)
			{
				Cost[i] = 0;
				d = 0;
			}
			else
			{
				//huC,bcC,huH,bcH,	huL,bcL,huLR,bcLR
				runMidpoint2D(xc,yc, Dx[i], Dy[i],3, X,OBJECT,BACKGROUND,&huC,&bcC);
				runMidpoint2D(x0,y0, Dx[i], Dy[i],3, X,OBJECT,BACKGROUND,&huH,&bcH);
				runMidpoint2D(x1,y1, Dx[i], Dy[i],3, X,OBJECT,BACKGROUND,&huL,&bcL);
				runMidpoint2D(x2,y2, Dx[i], Dy[i],3, X,OBJECT,BACKGROUND,&huLR,&bcLR);
				d = getMaxOf4((huC+0.01) /(bcC+0.01) ,(huH+0.01) /(bcH+0.01),(huL+0.01) /(bcL+0.01),(huLR+0.01) /(bcLR+0.01));//getMin
				Cost[i] = Cost[i]/d;
			}
			printf("Cost[%d] = %f d = %f dist = %f suf = %d\n",i,Cost[i],d,dist,D[i]);
		
		}

		free(D);
		free(Dx);
		free(Dy);
	}
	return (*maxL);
			
}
	

double getNextRegion(int maxL,double *Cost,int **V,int **X,double v,double nv,int lines,int cols)
{
	double min;
	int minID,xmin,ymin,xmax,ymax;
	int i,i1,i2;


	xmin = 2;
	ymin = 2;
	xmax = lines-2;
	ymax = cols-2;


	min = INF;
	minID = 1;
	for (i = 1; i < maxL; ++i)
	{
		if (Cost[i] != -1 && Cost[i] < min)
		{
			minID = i;
			min = Cost[i];
		}
	}
	
	Cost[minID] = -1;
	
	for (i1 = xmin; i1 < xmax; i1++)
	for (i2 = ymin; i2 < ymax; i2++)
		if (X[i1][i2] == v && V[i1][i2] == minID)
		{
			X[i1][i2] = (int)nv;
		}

	
	return min;
}

double getNextRegionMax(int maxL,double *Cost,int **V,int **X,double v,double nv,int lines,int cols)
{
	double max;
	int maxID,xmin,ymin,xmax,ymax;
	int i,i1,i2;


	xmin = 2;
	ymin = 2;
	xmax = lines-2;
	ymax = cols-2;


	max = -INF;
	maxID = 1;
	for (i = 1; i < maxL; ++i)
	{
		if (Cost[i] != -1 && Cost[i] > max)
		{
			maxID = i;
			max = Cost[i];
		}
	}
	
	Cost[maxID] = -1;
	
	for (i1 = xmin; i1 < xmax; i1++)
	for (i2 = ymin; i2 < ymax; i2++)
		if (X[i1][i2] == v && V[i1][i2] == maxID)
		{
			X[i1][i2] = (int)nv;
		}

	
	return max;
}

void findRegion(int **X,double v,double nv)
{
	int i,j,i1;
	int **V =(int **) Alloc2D(lines,sizeof(int*),cols,sizeof(int));
	int *D = (int *)calloc(100,sizeof(int));
	int maxL = 1;
	int maxV = 0,maxI = 0;
	int sizeD = 100;

	for (i = 0; i < lines; i++)
	for (j = 0; j < cols; j++)
	{
		if (X[i][j] == v && V[i][j] == 0)
		{
			getSufraceOfArea(X,V,v,i,j,maxL,&(D[maxL]));
			maxL++;
			if (maxL == sizeD-2)
			{
				sizeD = 2*sizeD;
				D = (int *)realloc(D,sizeD*sizeof(int));
				for (i1 = maxL; i1 < sizeD; ++i1)
					D[i1] = 0;
			}
		}
	}
			
	for (i = 1; i < maxL; ++i)
	{
		printf("E[%d] = %d\n",i,D[i]);
		if (D[i] > maxV)
		{
			
			maxV = D[i];
			maxI = i;
		}
	}

	for (i = 0; i < lines; i++)
	for (j = 0; j < cols; j++)
	{
		if (X[i][j] == v && V[i][j] != maxI)
		{
			X[i][j] = (int)nv;
		}
	}
	FreeInt(V,lines);
	free(D);
}



void addpin(double **A,int lines,int cols,double **B,double **C,char pr)
{
	int i,j;

	for (i = 0; i < lines; ++i)
	for (j = 0; j < cols; ++j)
	{
		if (pr == '-')
			C[i][j] = A[i][j]-B[i][j];		
		else
			C[i][j] = A[i][j]+B[i][j];		
			
	}

}


void Ison(int **A,int xmin,int xmax,int ymin,int ymax,int **C)
{
	int i,j;

	for (i = xmin; i <= xmax; ++i)
	for (j = ymin; j <= ymax; ++j)
	{
		C[i][j] = A[i][j];
	}

}

void Anastrofos(double **A,int lines,int cols,double **C)
{
	int i,j;

	for (i = 0; i < lines; ++i)
	for (j = 0; j < cols; ++j)
	{
		C[j][i] = A[i][j];
	}

}



void pollpin(double **A,int lines,int cols1,double **B,int lines2,int cols,double **C)
{
	int i,j,m;
	
	if (cols1 != lines2)
	{
		printf("La8os pollaplasiasmos pinakwn \n");
	}


	for (i = 0; i <lines; ++i)
		for (j = 0; j < cols; ++j)
			C[i][j] = 0;

	for (i = 0; i <lines; ++i)
		for (j = 0; j < cols; ++j)
		{
			for (m = 0; m < cols1; ++m)
			{
				C[i][j] += A[i][m]*B[m][j];
			}
		}


}

int existRoad(int x0, int y0, int x1, int y1, double VIn,double VIn2,int okVin2,int **image,int width) 
{
	int dy, dx, d, pos, neg, x, y, tmp, wd; 
	if (okVin2 == 0)
		VIn2 = (double)8210081.797349;

	if(x0 > x1) {
		SWAP(x0, x1, tmp) 
		SWAP(y0, y1, tmp)
	}  
	if(x0 == x1  && y0 > y1) SWAP(y0, y1, tmp)

	dx = x1 - x0;
	dy = y1 - y0;
  
	x = x0;
	y = y0;

	/** bput_point_c(0, 0, value); **/

	//bput_point_c(x, y, value);
	if (dy > dx || dy < -dx) {
		for (wd=0;wd<width;wd++) { /** left - rigth **/
    		if (wd%2==0) {
				//bput_point_c(x-(wd/2+wd%2),y, value); 
				if (image[x-(wd/2+wd%2)][y] != VIn && image[x-(wd/2+wd%2)][y] != VIn2)
					return 0;
			}
    		else {
				//bput_point_c(x+(wd/2+wd%2),y,value);
				if (image[x+(wd/2+wd%2)][y] != VIn && image[x+(wd/2+wd%2)][y] != VIn2)
					return 0;
			}
        }
	}
	else {
    	for (wd=0;wd<width;wd++) { /** up - down **/
			if (wd%2==0) {
				//bput_point_c(x,y-(wd/2+wd%2), value); 
				if (image[x][y-(wd/2+wd%2)] != VIn && image[x][y-(wd/2+wd%2)] != VIn2)
					return 0;
			}
            else {
				//bput_point_c(x,y+(wd/2+wd%2),value);
				if (image[x][y+(wd/2+wd%2)] != VIn)
					return 0;
			}
		}
	}


	if(dy >= 0)            
		if(dy > dx) {        //slope > 1                          
			pos = 2*(-dx);        //incrN			   
			neg = 2*(dy - dx);    //incrNE
			d = dy - 2*dx;

			while(y < y1) {
				if(d <= 0) {
					d += neg; 
					x++; 
					y++;
				}
				else {
					d += pos;
					y++;
				}
				for (wd=0;wd<width;wd++) { /** left - rigth **/
              		if (wd%2==0) {
						//bput_point_c(x-(wd/2+wd%2),y, value); 
						if (image[x-(wd/2+wd%2)][y] != VIn && image[x-(wd/2+wd%2)][y]   != VIn2)
							return 0;
						
					}
               		else {
						//bput_point_c(x+(wd/2+wd%2),y,value);
						if (image[x+(wd/2+wd%2)][y] != VIn && image[x+(wd/2+wd%2)][y]  != VIn2)
							return 0;	
					}
                }
			} 
		}

		else {               //0 <= slope <= 1 
			pos = 2*(dy - dx);    //incrNE			   
			neg = 2*dy;           //incrE
			d = 2*dy - dx; 

			while(x < x1) {
				if(d <= 0) {
					d += neg;  
					x++;
				}
				else  {
					d += pos;
					x++;
					y++;
				}
            	for (wd=0;wd<width;wd++) { /** up - down **/
					if (wd%2==0) {
						//bput_point_c(x,y-(wd/2+wd%2), value); 
						if (image[x][y-(wd/2+wd%2)] != VIn && image[x][y-(wd/2+wd%2)] != VIn2)
							return 0;
					}
                    else {
						//bput_point_c(x,y+(wd/2+wd%2),value);
						if (image[x][y+(wd/2+wd%2)] != VIn && image[x][y+(wd/2+wd%2)] != VIn2)
							return 0;
						
					}
				}
			}
		}

	else
		if(-dy > dx) {      //slope < -1
			pos = 2*(dy + dx);    //incrSE
			neg = 2*dx;           //incrS
			d = dy + 2*dx;

			while(y > y1) {
				if(d <= 0) {
					d += neg;  
					y--;
				}
				else {
					d += pos;
					x++;
					y--;
				}
				for (wd=0;wd<width;wd++) { /** left - right **/
					if (wd%2==0) {
						//bput_point_c(x+(wd/2+wd%2),y, value);
						if (image[x+(wd/2+wd%2)][y] != VIn && image[x+(wd/2+wd%2)][y] != VIn2)
							return 0;
					}
                    else {
						//bput_point_c(x-(wd/2+wd%2),y,value);
						if (image[x-(wd/2+wd%2)][y] != VIn && image[x-(wd/2+wd%2)][y]   != VIn2)
							return 0;
					}
				}
			} 
		}

		else {               //-1 <= slope < 0
			pos = 2*dy;         //incrE
			neg = 2*(dy + dx);       //incrSE
			d = 2*dy + dx; 

			while(x < x1) {
				if(d <= 0) {
					d += neg;  
					x++;
					y--;
				}
				else{
					d += pos;
					x++;
				}
				for (wd=0;wd<width;wd++) {
                	if (wd%2==0) {
						//bput_point_c(x,y-(wd/2+wd%2), value);
						if (image[x][y-(wd/2+wd%2)] != VIn && image[x][y-(wd/2+wd%2)]  != VIn2  )
							return 0;
					}
                    else {
						//bput_point_c(x,y+(wd/2+wd%2),value);	
						if (image[x][y+(wd/2+wd%2)] != VIn && image[x][y+(wd/2+wd%2)] != VIn2)
							return 0;
							
					}
                }
			} 
		}
	return 1;
}


void MidpointInVin(int x0, int y0, int x1, int y1, int width, int **image,int r0,int **Vin,int vin) 
{
	int dy, dx, d, pos, neg, x, y, tmp, wd; 
  
	if(x0 > x1) {
		SWAP(x0, x1, tmp) 
		SWAP(y0, y1, tmp)
	}  
	if(x0 == x1  && y0 > y1) SWAP(y0, y1, tmp)

	dx = x1 - x0;
	dy = y1 - y0;
  
	x = x0;
	y = y0;

	/** bput_point_c(0, 0, value); **/

	//bput_point_c(x, y, value);
	image[x][y] = r0;
	if (dy > dx || dy < -dx) {
		for (wd=0;wd<width;wd++) { /** left - rigth **/
    		if (wd%2==0) {
				if (Vin[x][y] == vin)
				//bput_point_c(x-(wd/2+wd%2),y, value); 
					image[x-(wd/2+wd%2)][y] = r0;
			}
    		else {
				//bput_point_c(x+(wd/2+wd%2),y,value);
				if (Vin[x+(wd/2+wd%2)][y] == vin)
					image[x+(wd/2+wd%2)][y] = r0;
			}
        }
	}
	else {
    	for (wd=0;wd<width;wd++) { /** up - down **/
			if (wd%2==0) {
				
				//bput_point_c(x,y-(wd/2+wd%2), value); 
				if (Vin[x][y-(wd/2+wd%2)] == vin)
					image[x][y-(wd/2+wd%2)] = r0;
			}
            else {
				//bput_point_c(x,y+(wd/2+wd%2),value);
				if (Vin[x][y+(wd/2+wd%2)] == vin)
					image[x][y+(wd/2+wd%2)] = r0;
			}
		}
	}


	if(dy >= 0)            
		if(dy > dx) {        //slope > 1                          
			pos = 2*(-dx);        //incrN			   
			neg = 2*(dy - dx);    //incrNE
			d = dy - 2*dx;

			while(y < y1) {
				if(d <= 0) {
					d += neg; 
					x++; 
					y++;
				}
				else {
					d += pos;
					y++;
				}
				for (wd=0;wd<width;wd++) { /** left - rigth **/
              		if (wd%2==0) {
						//bput_point_c(x-(wd/2+wd%2),y, value); 
						if (Vin[x-(wd/2+wd%2)][y] == vin)
							image[x-(wd/2+wd%2)][y] = r0;
					}
               		else {
						//bput_point_c(x+(wd/2+wd%2),y,value);
						if (Vin[x+(wd/2+wd%2)][y] == vin)
							image[x+(wd/2+wd%2)][y] = r0;
					}
                }
			} 
		}

		else {               //0 <= slope <= 1 
			pos = 2*(dy - dx);    //incrNE			   
			neg = 2*dy;           //incrE
			d = 2*dy - dx; 

			while(x < x1) {
				if(d <= 0) {
					d += neg;  
					x++;
				}
				else  {
					d += pos;
					x++;
					y++;
				}
            	for (wd=0;wd<width;wd++) { /** up - down **/
					if (wd%2==0) {
						//bput_point_c(x,y-(wd/2+wd%2), value); 
						if (Vin[x][y-(wd/2+wd%2)] == vin)
							image[x][y-(wd/2+wd%2)] = r0;
					}
                    else {
						//bput_point_c(x,y+(wd/2+wd%2),value);
						
						if (Vin[x][y+(wd/2+wd%2)] == vin)
							image[x][y+(wd/2+wd%2)] = r0;
					}
				}
			}
		}

	else
		if(-dy > dx) {      //slope < -1
			pos = 2*(dy + dx);    //incrSE
			neg = 2*dx;           //incrS
			d = dy + 2*dx;

			while(y > y1) {
				if(d <= 0) {
					d += neg;  
					y--;
				}
				else {
					d += pos;
					x++;
					y--;
				}
				for (wd=0;wd<width;wd++) { /** left - right **/
					if (wd%2==0) {
						//bput_point_c(x+(wd/2+wd%2),y, value);
						if (Vin[x+(wd/2+wd%2)][y] == vin)
							image[x+(wd/2+wd%2)][y] = r0;
					}
                    else {
						//bput_point_c(x-(wd/2+wd%2),y,value);
						if (Vin[x-(wd/2+wd%2)][y] == vin)
							image[x-(wd/2+wd%2)][y] = r0;
					}
				}
			} 
		}

		else {               //-1 <= slope < 0
			pos = 2*dy;         //incrE
			neg = 2*(dy + dx);       //incrSE
			d = 2*dy + dx; 

			while(x < x1) {
				if(d <= 0) {
					d += neg;  
					x++;
					y--;
				}
				else{
					d += pos;
					x++;
				}
				for (wd=0;wd<width;wd++) {
                	if (wd%2==0) {
						//bput_point_c(x,y-(wd/2+wd%2), value);
						if (Vin[x][y-(wd/2+wd%2)] == vin)
							image[x][y-(wd/2+wd%2)] = r0;
					}
                    else {
						//bput_point_c(x,y+(wd/2+wd%2),value);	
						if (Vin[x][y+(wd/2+wd%2)] == vin)
							image[x][y+(wd/2+wd%2)] = r0;
					}
                }
			} 
		}
}






int Midpoint2D(int x0, int y0, int x1, int y1, int width, int **image,int r0)
{
	int dy, dx, d, pos, neg, x, y, tmp, wd; 
	int count = 0;
  
	if(x0 > x1) {
		SWAP(x0, x1, tmp) 
		SWAP(y0, y1, tmp)
	}  
	if(x0 == x1  && y0 > y1) SWAP(y0, y1, tmp)

	dx = x1 - x0;
	dy = y1 - y0;
  
	x = x0;
	y = y0;

	/** bput_point_c(0, 0, value); **/

	//bput_point_c(x, y, value);
	image[x][y] = r0;
	if (dy > dx || dy < -dx) {
		for (wd=0;wd<width;wd++) { /** left - rigth **/
    		if (wd%2==0) {
				//bput_point_c(x-(wd/2+wd%2),y, value); 
				image[x-(wd/2+wd%2)][y] = r0;
				++count;
			}
    		else {
				//bput_point_c(x+(wd/2+wd%2),y,value);
				image[x+(wd/2+wd%2)][y] = r0;
				++count;
			}
        }
	}
	else {
    	for (wd=0;wd<width;wd++) { /** up - down **/
			if (wd%2==0) {
				//bput_point_c(x,y-(wd/2+wd%2), value); 
				image[x][y-(wd/2+wd%2)] = r0;
				++count;
			}
            else {
				//bput_point_c(x,y+(wd/2+wd%2),value);
				image[x][y+(wd/2+wd%2)] = r0;
				++count;
			}
		}
	}


	if(dy >= 0)            
		if(dy > dx) {        //slope > 1                          
			pos = 2*(-dx);        //incrN			   
			neg = 2*(dy - dx);    //incrNE
			d = dy - 2*dx;

			while(y < y1) {
				if(d <= 0) {
					d += neg; 
					x++; 
					y++;
				}
				else {
					d += pos;
					y++;
				}
				for (wd=0;wd<width;wd++) { /** left - rigth **/
              		if (wd%2==0) {
						//bput_point_c(x-(wd/2+wd%2),y, value); 
						image[x-(wd/2+wd%2)][y] = r0;
						++count;
					}
               		else {
						//bput_point_c(x+(wd/2+wd%2),y,value);
						image[x+(wd/2+wd%2)][y] = r0;
						++count;
					}
                }
			} 
		}

		else {               //0 <= slope <= 1 
			pos = 2*(dy - dx);    //incrNE			   
			neg = 2*dy;           //incrE
			d = 2*dy - dx; 

			while(x < x1) {
				if(d <= 0) {
					d += neg;  
					x++;
				}
				else  {
					d += pos;
					x++;
					y++;
				}
            	for (wd=0;wd<width;wd++) { /** up - down **/
					if (wd%2==0) {
						//bput_point_c(x,y-(wd/2+wd%2), value); 
						image[x][y-(wd/2+wd%2)] = r0;
						++count;
					}
                    else {
						//bput_point_c(x,y+(wd/2+wd%2),value);
						image[x][y+(wd/2+wd%2)] = r0;
						++count;
					}
				}
			}
		}

	else
		if(-dy > dx) {      //slope < -1
			pos = 2*(dy + dx);    //incrSE
			neg = 2*dx;           //incrS
			d = dy + 2*dx;

			while(y > y1) {
				if(d <= 0) {
					d += neg;  
					y--;
				}
				else {
					d += pos;
					x++;
					y--;
				}
				for (wd=0;wd<width;wd++) { /** left - right **/
					if (wd%2==0) {
						//bput_point_c(x+(wd/2+wd%2),y, value);
						image[x+(wd/2+wd%2)][y] = r0;
						++count;
					}
                    else {
						//bput_point_c(x-(wd/2+wd%2),y,value);
						image[x-(wd/2+wd%2)][y] = r0;
						++count;
					}
				}
			} 
		}

		else {               //-1 <= slope < 0
			pos = 2*dy;         //incrE
			neg = 2*(dy + dx);       //incrSE
			d = 2*dy + dx; 

			while(x < x1) {
				if(d <= 0) {
					d += neg;  
					x++;
					y--;
				}
				else{
					d += pos;
					x++;
				}
				for (wd=0;wd<width;wd++) {
                	if (wd%2==0) {
						//bput_point_c(x,y-(wd/2+wd%2), value);
						image[x][y-(wd/2+wd%2)] = r0;
				++count;
					}
                    else {
						//bput_point_c(x,y+(wd/2+wd%2),value);	
						image[x][y+(wd/2+wd%2)] = r0;
				++count;
					}
                }
			} 
		}
	return count;
}



void runMidpoint2D(int x0, int y0, int x1, int y1, int width, int **image,int r0,int r1,int *p0,int *p1)
{
	int dy, dx, d, pos, neg, x, y, tmp, wd; 
	int val;

	*p0 = 0;
	*p1 = 0;
	if(x0 > x1) {
		SWAP(x0, x1, tmp) 
		SWAP(y0, y1, tmp)
	}  
	if(x0 == x1  && y0 > y1) SWAP(y0, y1, tmp)

	dx = x1 - x0;
	dy = y1 - y0;
  
	x = x0;
	y = y0;



	//bput_point_c(x, y, value);
	//image[x][y] = r0;
	if (dy > dx || dy < -dx) {
		for (wd=0;wd<width;wd++) { // left - rigth 
    		if (wd%2==0) {
				val = image[x-(wd/2+wd%2)][y];
				//bput_point_c(x-(wd/2+wd%2),y, value); 
				if (val == r0)
					++(*p0);
				else if (val == r1)
					++(*p1);
			}
    		else {
				//bput_point_c(x+(wd/2+wd%2),y,value);
				val = image[x+(wd/2+wd%2)][y];
				if (val == r0)
					++(*p0);
				else if (val == r1)
					++(*p1);
			}
        }
	}
	else {
    	for (wd=0;wd<width;wd++) { // up - down 
			if (wd%2==0) {
				//bput_point_c(x,y-(wd/2+wd%2), value); 
				val = image[x][y-(wd/2+wd%2)];

				if (val == r0)
					++(*p0);
				else if (val == r1)
					++(*p1);
				
				
			}
            else {
				//bput_point_c(x,y+(wd/2+wd%2),value);
				val = image[x][y+(wd/2+wd%2)];
				if (val == r0)
					++(*p0);
				else if (val == r1)
					++(*p1);
				
			}
		}
	}


	if(dy >= 0)            
		if(dy > dx) {        //slope > 1                          
			pos = 2*(-dx);        //incrN			   
			neg = 2*(dy - dx);    //incrNE
			d = dy - 2*dx;

			while(y < y1) {
				if(d <= 0) {
					d += neg; 
					x++; 
					y++;
				}
				else {
					d += pos;
					y++;
				}
				for (wd=0;wd<width;wd++) { // left - rigth
              		if (wd%2==0) {
						//bput_point_c(x-(wd/2+wd%2),y, value); 
						val = image[x-(wd/2+wd%2)][y];
						if (val == r0)
							++(*p0);
						else if (val == r1)
							++(*p1);
						
					}
               		else {
						//bput_point_c(x+(wd/2+wd%2),y,value);
						val = image[x+(wd/2+wd%2)][y];
						if (val == r0)
							++(*p0);
						else if (val == r1)
							++(*p1);
					
					}
                }
			} 
		}

		else {               //0 <= slope <= 1 
			pos = 2*(dy - dx);    //incrNE			   
			neg = 2*dy;           //incrE
			d = 2*dy - dx; 

			while(x < x1) {
				if(d <= 0) {
					d += neg;  
					x++;
				}
				else  {
					d += pos;
					x++;
					y++;
				}
            	for (wd=0;wd<width;wd++) { // up - down 
					if (wd%2==0) {
						//bput_point_c(x,y-(wd/2+wd%2), value); 
						val = image[x][y-(wd/2+wd%2)];
						if (val == r0)
							++(*p0);
						else if (val == r1)
							++(*p1);
					
					}
                    else {
						//bput_point_c(x,y+(wd/2+wd%2),value);
						val = image[x][y+(wd/2+wd%2)];
						if (val == r0)
							++(*p0);
						else if (val == r1)
							++(*p1);
					
					}
				}
			}
		}

	else
		if(-dy > dx) {      //slope < -1
			pos = 2*(dy + dx);    //incrSE
			neg = 2*dx;           //incrS
			d = dy + 2*dx;

			while(y > y1) {
				if(d <= 0) {
					d += neg;  
					y--;
				}
				else {
					d += pos;
					x++;
					y--;
				}
				for (wd=0;wd<width;wd++) { // left - right 
					if (wd%2==0) {
						//bput_point_c(x+(wd/2+wd%2),y, value);
						val = image[x+(wd/2+wd%2)][y];
						if (val == r0)
							++(*p0);
						else if (val == r1)
							++(*p1);
					
					}
                    else {
						//bput_point_c(x-(wd/2+wd%2),y,value);
						val = image[x-(wd/2+wd%2)][y];
						if (val == r0)
							++(*p0);
						else if (val == r1)
							++(*p1);
					
					}
				}
			} 
		}

		else {               //-1 <= slope < 0
			pos = 2*dy;         //incrE
			neg = 2*(dy + dx);       //incrSE
			d = 2*dy + dx; 

			while(x < x1) {
				if(d <= 0) {
					d += neg;  
					x++;
					y--;
				}
				else{
					d += pos;
					x++;
				}
				for (wd=0;wd<width;wd++) {
                	if (wd%2==0) {
						//bput_point_c(x,y-(wd/2+wd%2), value);
						val = image[x][y-(wd/2+wd%2)];
						if (val == r0)
							++(*p0);
						else if (val == r1)
							++(*p1);
					}
                    else {
						//bput_point_c(x,y+(wd/2+wd%2),value);	
						val = image[x][y+(wd/2+wd%2)];
						if (val == r0)
							++(*p0);
						else if (val == r1)
							++(*p1);
					}
                }
			} 
		}
}




void Midpoint(int x0, int y0, int x1, int y1, int width, int ***image,int r0,int g0,int b0) 
{
	int dy, dx, d, pos, neg, x, y, tmp, wd; 
  
	if(x0 > x1) {
		SWAP(x0, x1, tmp) 
		SWAP(y0, y1, tmp)
	}  
	if(x0 == x1  && y0 > y1) SWAP(y0, y1, tmp)

	dx = x1 - x0;
	dy = y1 - y0;
  
	x = x0;
	y = y0;

	/** bput_point_c(0, 0, value); **/

	//bput_point_c(x, y, value);
	image[x][y][0] = r0;
	image[x][y][1] = g0;
	image[x][y][2] = b0;
	if (dy > dx || dy < -dx) {
		for (wd=0;wd<width;wd++) { /** left - rigth **/
    		if (wd%2==0) {
				//bput_point_c(x-(wd/2+wd%2),y, value); 
				image[x-(wd/2+wd%2)][y][0] = r0;
				image[x-(wd/2+wd%2)][y][1] = b0;
				image[x-(wd/2+wd%2)][y][2] = g0;
			}
    		else {
				//bput_point_c(x+(wd/2+wd%2),y,value);
				image[x+(wd/2+wd%2)][y][0] = r0;
				image[x+(wd/2+wd%2)][y][1] = g0;
				image[x+(wd/2+wd%2)][y][2] = b0;
			}
        }
	}
	else {
    	for (wd=0;wd<width;wd++) { /** up - down **/
			if (wd%2==0) {
				//bput_point_c(x,y-(wd/2+wd%2), value); 
				image[x][y-(wd/2+wd%2)][0] = r0;
				image[x][y-(wd/2+wd%2)][1] = g0;
				image[x][y-(wd/2+wd%2)][2] = b0;
			}
            else {
				//bput_point_c(x,y+(wd/2+wd%2),value);
				image[x][y+(wd/2+wd%2)][0] = r0;
				image[x][y+(wd/2+wd%2)][1] = g0;
				image[x][y+(wd/2+wd%2)][2] = b0;
			}
		}
	}


	if(dy >= 0)            
		if(dy > dx) {        //slope > 1                          
			pos = 2*(-dx);        //incrN			   
			neg = 2*(dy - dx);    //incrNE
			d = dy - 2*dx;

			while(y < y1) {
				if(d <= 0) {
					d += neg; 
					x++; 
					y++;
				}
				else {
					d += pos;
					y++;
				}
				for (wd=0;wd<width;wd++) { /** left - rigth **/
              		if (wd%2==0) {
						//bput_point_c(x-(wd/2+wd%2),y, value); 
						image[x-(wd/2+wd%2)][y][0] = r0;
						image[x-(wd/2+wd%2)][y][1] = g0;
						image[x-(wd/2+wd%2)][y][2] = b0;
					}
               		else {
						//bput_point_c(x+(wd/2+wd%2),y,value);
						image[x+(wd/2+wd%2)][y][0] = r0;
						image[x+(wd/2+wd%2)][y][1] = g0;
						image[x+(wd/2+wd%2)][y][2] = b0;
					}
                }
			} 
		}

		else {               //0 <= slope <= 1 
			pos = 2*(dy - dx);    //incrNE			   
			neg = 2*dy;           //incrE
			d = 2*dy - dx; 

			while(x < x1) {
				if(d <= 0) {
					d += neg;  
					x++;
				}
				else  {
					d += pos;
					x++;
					y++;
				}
            	for (wd=0;wd<width;wd++) { /** up - down **/
					if (wd%2==0) {
						//bput_point_c(x,y-(wd/2+wd%2), value); 
						image[x][y-(wd/2+wd%2)][0] = r0;
						image[x][y-(wd/2+wd%2)][1] = g0;
						image[x][y-(wd/2+wd%2)][2] = b0;
					}
                    else {
						//bput_point_c(x,y+(wd/2+wd%2),value);
						image[x][y+(wd/2+wd%2)][0] = r0;
						image[x][y+(wd/2+wd%2)][1] = g0;
						image[x][y+(wd/2+wd%2)][2] = b0;
					}
				}
			}
		}

	else
		if(-dy > dx) {      //slope < -1
			pos = 2*(dy + dx);    //incrSE
			neg = 2*dx;           //incrS
			d = dy + 2*dx;

			while(y > y1) {
				if(d <= 0) {
					d += neg;  
					y--;
				}
				else {
					d += pos;
					x++;
					y--;
				}
				for (wd=0;wd<width;wd++) { /** left - right **/
					if (wd%2==0) {
						//bput_point_c(x+(wd/2+wd%2),y, value);
						image[x+(wd/2+wd%2)][y][0] = r0;
						image[x+(wd/2+wd%2)][y][1] = g0;
						image[x+(wd/2+wd%2)][y][2] = b0;
					}
                    else {
						//bput_point_c(x-(wd/2+wd%2),y,value);
						image[x-(wd/2+wd%2)][y][0] = r0;
						image[x-(wd/2+wd%2)][y][1] = g0;
						image[x-(wd/2+wd%2)][y][2] = b0;
					}
				}
			} 
		}

		else {               //-1 <= slope < 0
			pos = 2*dy;         //incrE
			neg = 2*(dy + dx);       //incrSE
			d = 2*dy + dx; 

			while(x < x1) {
				if(d <= 0) {
					d += neg;  
					x++;
					y--;
				}
				else{
					d += pos;
					x++;
				}
				for (wd=0;wd<width;wd++) {
                	if (wd%2==0) {
						//bput_point_c(x,y-(wd/2+wd%2), value);
						image[x][y-(wd/2+wd%2)][0] = r0;
						image[x][y-(wd/2+wd%2)][1] = g0;
						image[x][y-(wd/2+wd%2)][2] = b0;
					}
                    else {
						//bput_point_c(x,y+(wd/2+wd%2),value);	
						image[x][y+(wd/2+wd%2)][0] = r0;
						image[x][y+(wd/2+wd%2)][1] = g0;
						image[x][y+(wd/2+wd%2)][2] = b0;
					}
                }
			} 
		}
}


void print_box3(int ***A,int x1,int y1,int x2,int y2,int r0,int g0,int b0)
{
	
	int i,j;

	if (x1 < 0)
		x1 = 0;
	if (y1 < 0)
		y1 = 0;

	if (x2 >= lines)
		x2 = lines-1;
	if (y2 >= cols)
		y2 = cols-1;

	if (x2 < 0)
		x2 = 0;
	if (y2 < 0)
		y2 = 0;

	if (x1 >= lines)
		x1 = lines-1;
	if (y1 >= cols)
		y1 = cols-1;

	
	for (i = x1,j = y1; j <= y2; ++j)
	{
		A[i][j][0] = r0;			
		A[i][j][1] = g0;			
		A[i][j][2] = b0;			

	}	
	
	for (i = x1,j = y1; i <= x2; ++i)
	{
		A[i][j][0] = r0;			
		A[i][j][1] = g0;			
		A[i][j][2] = b0;			

	}	
	
	for (i = x1,j = y2; i <= x2; ++i)
	{
		A[i][j][0] = r0;			
		A[i][j][1] = g0;			
		A[i][j][2] = b0;			

	}	
	
	for (i = x2,j = y1; j <= y2; ++j)
	{
		A[i][j][0] = r0;			
		A[i][j][1] = g0;			
		A[i][j][2] = b0;			

	}

}


void print_box2(int **A,int x1,int y1,int x2,int y2,int r0)
{
	
	int i,j;

	if (x1 < 0)
		x1 = 0;
	if (y1 < 0)
		y1 = 0;

	if (x2 >= lines)
		x2 = lines-1;
	if (y2 >= cols)
		y2 = cols-1;

	if (x2 < 0)
		x2 = 0;
	if (y2 < 0)
		y2 = 0;

	if (x1 >= lines)
		x1 = lines-1;
	if (y1 >= cols)
		y1 = cols-1;

	
	for (i = x1,j = y1; j <= y2; ++j)
	{
		A[i][j] = r0;			

	}	
	
	for (i = x1,j = y1; i <= x2; ++i)
	{
		A[i][j] = r0;			

	}	
	
	for (i = x1,j = y2; i <= x2; ++i)
	{
		A[i][j] = r0;			

	}	
	
	for (i = x2,j = y1; j <= y2; ++j)
	{
		A[i][j] = r0;			

	}

}


int KentroMazas(int ***X,int r0,int g0,int b0,int *XC,int *YC,int **Y)
{
	int i,j;
	double kx = 0.0,ky = 0.0,num = 0.0;

	*XC = -1;
	*YC = -1;

	if (X == 0)
	{
		for (i = 0; i < lines; ++i)
		for (j = 0; j < cols; ++j)
		{
			if (Y[i][j] == r0)
			{
				//printf("%d %d %d\n",i,j,r0);
				kx += i;
				ky += j;
				++num;
			}
		}
		
	
		if (num != 0)
		{
			*XC = (int)(kx / num);
			*YC = (int)(ky / num);
		}	
		return (int)num;
	}
		
	


	for (i = 0; i < lines; ++i)
	for (j = 0; j < cols; ++j)
	{
		if (X[i][j][0] == r0 && X[i][j][1] == g0 && X[i][j][2] == b0)
		{
			kx += i;
			ky += j;
			++num;
		}
		
	}
	if (num != 0)
	{
		*XC = (int)(kx / num);
		*YC = (int)(ky / num);
	}	
	return (int)num;
}


double OrientationC(int **X,double r0,double *eken,int choice,int *xcc,int *ycc)
{
	int i,j;
	double m20=0,m11=0,m02=0,fi,area = 0;
	double xc = 0,yc = 0,sum = (double)0.0001;
	double XC = 0;

	for (i = 0; i < lines; ++i)
	for (j = 0; j < cols; ++j)
	{
			if (X[i][j] == r0)
			{
				XC += i;
				sum++;
			}
	}

	XC = XC / sum;
	sum = (double)0.0001;	
	for (i = 0; i < lines; ++i)
	for (j = 0; j < cols; ++j)
	{
		if (choice == 1)
		{
			if (i < XC && X[i][j] == r0)
			{
				xc += i;
				yc += j;
				sum++;
			}
		}
		else
		{
			if (i >= XC && X[i][j] == r0)
			{
				xc += i;
				yc += j;
				sum++;
			}
		}
	}
	xc = xc / sum;
	yc = yc / sum;	
	*xcc = (int)(xc+0.5);
	*ycc = (int)(yc+0.5);

	if (choice == 1)
	{	
		for (i = 0; i < lines; ++i)
		for (j = 0; j < cols; ++j)
		if (X[i][j] == r0 && i < XC)
		{
			area += 1;
			m11 += (-i+xc)*(j-yc); 		
			m02 += (-i+xc)*(-i+xc);		
			m20 += (j-yc)*(j-yc);	

		}
	}
	else
	{
		for (i = 0; i < lines; ++i)
		for (j = 0; j < cols; ++j)
		if (X[i][j] == r0 && i >= XC)
		{
			area += 1;
			m11 += (-i+xc)*(j-yc); 		
			m02 += (-i+xc)*(-i+xc);		
			m20 += (j-yc)*(j-yc);	

		}
	}
	
	fi = (double)(0.5*atan2(2*m11,m20-m02));

	*eken = ((m20-m02)*(m20-m02)+4*m11) / area;
	

	//printf("ek = %f\n",*eken);
	return fi;
}





double Orientation(int **X,int r0,double xc,double yc,double *eken)
{
	int i,j;
	double m20=0,m11=0,m02=0,fi,area = 0;
	double Imin = 0, Imax = 0;

	for (i = 0; i < lines; ++i)
	for (j = 0; j < cols; ++j)
		if (X[i][j] == r0 )
		{
			area += 1;
			m11 += (-i+xc)*(j-yc); 		
			m02 += (-i+xc)*(-i+xc);		
			m20 += (j-yc)*(j-yc);	

		}
	
	fi = (double)(0.5*atan2(2*m11,m20-m02));
	//*eken = ((m20-m02)*(m20-m02)+4*m11)/area;

	for (i = 2; i < lines-2; ++i)
	for (j = 2; j < cols-2; ++j)
		if (X[i][j] == r0 && (X[i+1][j] != r0 || X[i][j+1] != r0 || X[i+1][j+1] != r0 || X[i-1][j] != r0 || X[i][j-1] != r0 || X[i-1][j-1] != r0))
		{
			Imin = Imin + SQUARE((-i+xc)*cos(fi) - (j-yc)*sin(fi));
			Imax = Imax + SQUARE((-i+xc)*sin(fi) + (j-yc)*cos(fi));
		}

	
	if (area <= 2)
		*eken = 10;
	else
	{
	//	*eken = (double)((m20+m02+sqrt((m20-m02)*(m20-m02)+4*m11*m11)) / (m20+m02-sqrt((m20-m02)*(m20-m02)+4*m11*m11)));
		//*eken = MAX(m20,m02) / MIN(m20,m02);
		*eken = (double) sqrt(Imax/Imin);	
	}

	printf("ekN = %f M = %f %f %f X = %f %f %d %d a = %f\n",*eken,m02,m20,m11,xc,yc,lines,cols,area);
	//if (*eken < 2)
		//getchar();

	return (double)fi;
}



void projection(int **X,double r0,int hor,double *res)
{
	int i,j;

	if (hor == 1)
	{
		for (i = 0; i < lines; ++i)
			res[i] = 0;
		for (i = 0; i < lines; ++i)
		for (j = 0; j < cols; ++j)
		{
			if (X[i][j] != 0)
				res[i] += 1;
		}	
		return;
	}

	for (i = 0; i < cols; ++i)
		res[i] = 0;
	for (j = 0; j < cols; ++j)
	for (i = 0; i < lines; ++i)
	{
		if (X[i][j] != 0)
			res[j] += 1;
	}	
	return;

}



void RotF(double fi,double x0,double y0,double *x1,double *y1)
{

        *x1 = (double)((cos(fi)*x0)+(-sin(fi)*y0));
        *y1 = (double)((cos(fi)*y0)+(sin(fi)*x0));
}



/* <xn,yn> = R(fi,<x0,y0>) . <i,j> */
void RotatePointF(double i,double j,double fi,double xc,double yc,double *xn,double *yn)
{
	double x,y;

	RotF(fi,i-xc,j-yc,&x,&y);
	x += xc;
	y += yc;
	*xn = x;
	*yn = y;

}



void Rot(double fi,int x0,int y0,double *x1,double *y1)
{

        *x1 = (double)((cos(fi)*x0)+(-sin(fi)*y0));
        *y1 = (double)((cos(fi)*y0)+(sin(fi)*x0));
}

int fround(double f)
{
        int h;

        h = (int)f;

        if (f - h > 0.5)
                h = h+1;

        return h;

}


/* <xn,yn> = R(fi,<x0,y0>) . <i,j> */
void RotatePoint(int i,int j,double fi,int xc,int yc,int *xn,int *yn)
{
	double x,y;

	Rot(fi,i-xc,j-yc,&x,&y);
	x += xc;
	y += yc;

	if (x >= 0 && x < lines && y >= 0 && y < cols)
	{
		*xn = fround(x);
		*yn = fround(y);
	}
	else
	{
		*xn = 0;
		*yn = 0;	
		printf("error Rotation\n");
	}
}




/* Y = R(fi1) . X  an X.x > xc */
/* Y = R(fi2) . X  an X.x  <= xc */



void RotateC(double **X,double **Y,double fi1,double fi2,int xc1,int yc1,int xc2,int yc2,double r0,double sunt,int XC)
{
	int i,j,i1,j1,i2,j2;
	double x,y;
	double di;
	double **baroi;
	double fi;
	int xc,yc;
	
	baroi = TempF;

	for (i = 0; i < lines; ++i)
	for (j = 0; j < cols; ++j)
	{
		baroi[i][j] = 0;
		Y[i][j] = 0;
	}

	for (i = 0; i < lines; ++i)
	for (j = 0; j < cols; ++j)
	{
		if (X[i][j] != r0)
		{
			if (i < XC)
			{
				fi = fi1;
				xc = xc1;
				yc = yc1;
			}
			else
			{
				fi = fi2;
				xc = xc2;
				yc = yc2;
			}

			Rot(fi,i-xc,j-yc,&x,&y);
			x += xc;
			y += yc;

			if (x >= 0 && x < lines && y >= 0 && y < cols)
			{
				i1 = fround(x);
				j1 = fround(y);
				for (i2 = i1-1; i2 <= i1+1; ++i2)
					for (j2 = j1-1; j2 <= j1+1; ++j2)
					{
						if (i2 >= 0 && i2 < lines && j2 >=0 && j2 < cols )	
						{
							di = (i2-x)*(i2-x) + (j2-y)*(j2-y);
							if (di < 2)
							{
								baroi[i2][j2] += 1-(di/2);
								Y[i2][j2] += X[i][j]*(1-(di/2));
							}
						}
					}
			}
		}
	}



	for (i = 0; i < lines; ++i)
	for (j = 0; j < cols; ++j)
	{
		
		if (baroi[i][j] > 0)
		{
			if ( baroi[i][j] >= sunt)
			{
				Y[i][j] = fround(Y[i][j] / baroi[i][j]);	
			}
			else
			{
				Y[i][j] = r0;	
			}
		}
		else
		{
			Y[i][j] = r0;	
		}

	}

}





/* Y = R(fi) . X */


void Rotate(int **X,int **RY,double fi,int xc,int yc,double r0,double sunt)
{
	int i,j,i1,j1,i2,j2;
	double x,y;
	double di;
	double **baroi = TempF;
	double **Y;

	Y = (double **)Alloc2D(lines,sizeof(double *),cols,sizeof(double));	
	for (i = 0; i < lines; ++i)
	for (j = 0; j < cols; ++j)
		baroi[i][j] = 0;

	for (i = 0; i < lines; ++i)
	for (j = 0; j < cols; ++j)
	{
		if (X[i][j] != r0)
		{
			Rot(fi,i-xc,j-yc,&x,&y);
			x += xc;
			y += yc;

			if (x >= 0 && x < lines && y >= 0 && y < cols)
			{
				i1 = fround(x);
				j1 = fround(y);
				for (i2 = i1-1; i2 <= i1+1; ++i2)
					for (j2 = j1-1; j2 <= j1+1; ++j2)
					{
						if (i2 >= 0 && i2 < lines && j2 >=0 && j2 < cols )	
						{
							di = (i2-x)*(i2-x) + (j2-y)*(j2-y);
							if (di < 2)
							{
								baroi[i2][j2] += 1-(di/2);
								Y[i2][j2] += X[i][j]*(1-(di/2));
							}
						}
					}
			}
		}
	}



	for (i = 0; i < lines; ++i)
	for (j = 0; j < cols; ++j)
	{
		if (baroi[i][j] > 0)
		{
			if ( baroi[i][j] >= sunt)
			{
				Y[i][j] = fround((double)Y[i][j] / baroi[i][j]);	
			}
			else
			{
				Y[i][j] = r0;	
			}
		}
		else
		{
			Y[i][j] = r0;	
		}
	}

	for (i = 0; i < lines; ++i)
	for (j = 0; j < cols; ++j)
	{
		RY[i][j] = (int)Y[i][j];
	}	

	for (i = 0; i < lines; ++i)
		free(Y[i]);

	free(Y);
}

void deletesunoro(struct Sunoro *A)
{
	struct Point *temp,*temp1;

	temp = A->Head;

	while (A->Head != 0)
	{
		temp1 = temp->next;
		deletepoint(temp,A);
		temp = temp1;
		
	}
	A->Head = A->Tail = 0;
	free(A);
}


struct Sunoro *newSunoro()
{
	struct Sunoro *A = (struct Sunoro *)malloc(sizeof(struct Sunoro));
	A->Head = 0;
	A->Tail = 0;
	A->size = 0;
	
	return A;
}


void addpoint(struct Point *new,struct Sunoro *A)
{

	if (A->Tail == 0)
	{
		A->size++;
		A->Tail = A->Head = new;
		new->next = new->prev = 0;
		return;	
	}
	new->prev = A->Tail;
	A->Tail->next = new;
	new->next = 0;
	A->Tail = new;
	A->size++;
}

void deletepoint(struct Point *del,struct Sunoro *A)
{
	A->size--;
	if (del == A->Head)
	{
		if (A->Head->next != 0)	
			A->Head = A->Head->next;
		else
		{
			free(del);
			A->Tail = A->Head = 0;
			return;
		}
		
		

		A->Head->prev = 0;
		del->next = 0;
		del->prev = 0;
		free(del);
		return;
	}
	
	if (del == A->Tail)
	{
		if (A->Tail->prev != 0)	
			A->Tail = A->Tail->prev;
		else
		{
			A->Head = A->Tail = 0;
			free(del);
			return;
		}

		A->Tail->next = 0;
		del->next = 0;
		del->prev = 0;
		free(del);
		return;
	}
	if (del->prev != 0)	
		del->prev->next = del->next;
	if (del->next != 0)
		del->next->prev = del->prev;		

	del->next = del->prev = 0;
	free(del);

}


int getMin3Dint(int ***Temp,int xapo,int xeos,int yapo,int yeos,int L)
{
	int i,j;
	int min = 10000000;

	for (i = xapo; i <= xeos; ++i)
	for (j = yapo; j <= yeos; ++j)
	{
		if (Temp[i][j][L] < min)
			min = Temp[i][j][L];
	}
	return min;
}




int getMax3Dint(int ***Temp,int xapo,int xeos,int yapo,int yeos,int L)
{
	int i,j;
	int max = -10000000;

	for (i = xapo; i <= xeos; ++i)
	for (j = yapo; j <= yeos; ++j)
	{
		if (Temp[i][j][L] > max)
			max = Temp[i][j][L];
	}
	return max;
}




double getMax3D(double ***Temp,int xapo,int xeos,int yapo,int yeos,int zapo,int zeos,int *xp,int *yp,int *zp)
{
	int i,j,k;
	double max = -10000000;

	for (i = xapo; i <= xeos; ++i)
	for (j = yapo; j <= yeos; ++j)
	for (k = zapo; k <= zeos; ++k)
	{
		if (i >= 0 && i < lines && j >= 0 && j < cols && Temp[i][j][k] > max)
		{
			*xp = i;
			*yp = j;
			*zp = k;
			max = Temp[i][j][k];
		}	

	}
	return max;
}






double getMin(int **Temp,int xapo,int xeos,int yapo,int yeos,int *xp,int *yp)
{
	int i,j;
	double min = 10000000;

	for (i = xapo; i <= xeos; ++i)
	for (j = yapo; j <= yeos; ++j)
	{
		if (Temp[i][j] < min)
		{
			*xp = i;
			*yp = j;
			min = Temp[i][j];
		}	

	}
	return min;
}

int getNearestArea(int **X,int x0,int y0,int V)
{
	int i,i1,i2;
	double dist,mindist =100000000.0 ;


	for (i = 1; i < lines; ++i)
	{
		dist = 100000000.0;
		for (i1 = x0-i,i2 = y0 - i; i1 < x0 + i; ++i1)
		{
			if (X[i1][i2] >= V)
				dist = (x0-i1)*(x0-i1) + (y0-i2)*(y0-i2);
			if (dist < mindist)
				mindist = dist;

		}
		for (i1 = x0+i,i2 = y0 - i; i2 < y0 + i; ++i2)
		{
			if (X[i1][i2] >= V)
				dist = (x0-i1)*(x0-i1) + (y0-i2)*(y0-i2);
			if (dist < mindist)
				mindist = dist;

		}


		for (i1 = x0+i,i2 = y0 + i; i1 > x0 - i; --i1)
		{
			if (X[i1][i2] >= V)
				dist = (x0-i1)*(x0-i1) + (y0-i2)*(y0-i2);
			if (dist < mindist)
				mindist = dist;

		}

		for (i1 = x0-i,i2 = y0 + i; i2 > y0 - i; --i2)
		{
			if (X[i1][i2] >= V)
				dist = (x0-i1)*(x0-i1) + (y0-i2)*(y0-i2);
			if (dist < mindist)
				mindist = dist;

		}

		if (100000000.0 != mindist)
			if (i*i > mindist)
				return (int)(sqrt(mindist));
	}
	return -1;
}









double getNearestPointFormArea(int **X,int x0,int y0,int V,int *xp,int *yp)
{
	int i,i1,i2;
	double dist,mindist = 100000000.0;
	int x = -1,y = -1;

	
	for (i = 1; i < lines; ++i)
	{
		dist = 100000000.0;
		for (i1 = x0-i,i2 = y0 - i; i1 < x0 + i; ++i1)
		{
			if (i1 < 0 || i2 < 0 || i1 >= lines || i2 >= cols)
				continue;

			if (X[i1][i2] == V)
			{
				dist = (x0-i1)*(x0-i1) + (y0-i2)*(y0-i2);
			
				if (dist < mindist)
				{
					mindist = dist;
					x = i1;
					y = i2;
				}
			}

		}
		for (i1 = x0+i,i2 = y0 - i; i2 < y0 + i; ++i2)
		{
			if (i1 < 0 || i2 < 0 || i1 >= lines || i2 >= cols)
				continue;
			if (X[i1][i2] == V)
			{
				dist = (x0-i1)*(x0-i1) + (y0-i2)*(y0-i2);
				if (dist < mindist)
				{
					mindist = dist;
					x = i1;
					y = i2;
				}
			}


		}


		for (i1 = x0+i,i2 = y0 + i; i1 > x0 - i; --i1)
		{
			if (i1 < 0 || i2 < 0 || i1 >= lines || i2 >= cols)
				continue;
			if (X[i1][i2] == V)
			{
				dist = (x0-i1)*(x0-i1) + (y0-i2)*(y0-i2);
				if (dist < mindist)
				{
					mindist = dist;
					x = i1;
					y = i2;
				}
			}


		}

		for (i1 = x0-i,i2 = y0 + i; i2 > y0 - i; --i2)
		{
			if (i1 < 0 || i2 < 0 || i1 >= lines || i2 >= cols)
				continue;
			if (X[i1][i2] == V)
			{
				dist = (x0-i1)*(x0-i1) + (y0-i2)*(y0-i2);
				if (dist < mindist)
				{
					mindist = dist;
					x = i1;
					y = i2;
				}
			}


		}

		if (100000000.0 != mindist)
		{
			if (i*i > mindist)
			{
				*xp = x;
				*yp = y;
				//printf("OUT = %d %d %d\n",x,y,V);
				return (double)(sqrt(mindist));
			}
		}
	}
	*xp = *yp = -1;
	return -1;
}


void swaptab( double **v, int i, int j)
{
        double temp[3];   

        temp[0] = v[i][0];
        temp[1] = v[i][1];
        temp[2] = v[i][2];

        v[i][0] = v[j][0];
        v[i][1] = v[j][1];
        v[i][2] = v[j][2];

        v[j][0] = temp[0];
        v[j][1] = temp[1];
        v[j][2] = temp[2];
}

void quicksort( double **v, int left, int right )
{
        int i, last;   

        if ( left >= right )
                return;
        swaptab( v, left, (left+right)/2 );
        last = left;
        for( i=left+1; i<=right; i++ )
                if( v[i][0] < v[left][0] )
                        swaptab( v, ++last, i );
        swaptab( v, left, last );
        quicksort( v, left, last-1 );
        quicksort( v, last+1, right );
}


int getPrice(double **R,int apo,int eos,double rate)
{
	int f;

	if (apo > eos)
		return -1;
	quicksort( R,apo,eos);
	
	f =  (int)(apo+(eos-apo)*rate+0.5); 

	if (f > eos)
		return eos;
	return f;
}


void diabrwsh(int **X,int **Y,int rate,double pos)
{
	int i1,i2,v1,i,j;

	for (i = rate; i < lines-rate; ++i)
	for (j = rate; j < cols-rate; ++j)
	{
		Y[i][j] = 0;
		v1 = 0;
		for (i1 = -rate; i1 <= rate; ++i1)
		for (i2 = -rate; i2 <= rate; ++i2)
			if (X[i+i1][j+i2] == 1)
				v1++;
				
		if (v1 > pos*(2*rate+1)*(2*rate+1))
			Y[i][j] = 1;

	}
}

int **zoom(int **A,int lines,int cols,int epix,int epiy)
{
	int newlines,newcols,i1,i2,i,j;
	int **res;

	newlines = lines*epix;
	newcols = cols*epiy;
	res = (int **)Alloc2D(newlines,sizeof(int *),newcols,sizeof(int));

	for (i = 0; i < lines; ++i)	
	for (j = 0; j < cols; ++j)
	{
		for (i1 = 0; i1 < epix; ++i1)
		for (i2 = 0; i2 < epiy; ++i2)
			res[i*epix+i1][j*epiy+i2] = A[i][j];
	}	

	return res;
}



void swaptab1D( double *v, int i, int j)
{
        double temp;   

        temp = v[i];

        v[i] = v[j];
        v[j] = temp;
}

void quicksort1D( double *v, int left, int right )
{
        int i, last;   

        if ( left >= right )
                return;
        swaptab1D( v, left, (left+right)/2 );
        last = left;
        for( i=left+1; i<=right; i++ )
                if( v[i] < v[left] )
                        swaptab1D( v, ++last, i );
        swaptab1D( v, left, last );
        quicksort1D( v, left, last-1 );
        quicksort1D( v, last+1, right );
}


