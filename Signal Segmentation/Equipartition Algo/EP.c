/**
 * @author Panagiotakis Konstantinos
 * @version     v1.0
 */

/* EP method*/

#include "EP.h"

//Returns 1 if i1 and i2 are connected
int areConnected(int level,int i1, int i2)
{
	int i,j;

	for (i = 0; i < ep->p[level][i1].degree; ++i)
	for (j = 0; j < ep->p[level][i2].degree; ++j)
	{
		if (ep->p[level][i1].line[i] == ep->p[level][i2].line[j])
			return 1;
	}

	return 0;
}


/*
	Returns the value of y of e.t at x 
*/
double getValueOfLineY(double x0,double y0, double x1,double y1,double x)
{
	double y = y0;
	double r;

	if (x0 != x1)
	{
		 r = (x-x0) / (x1-x0);
		 y = y0 + r*(y1-y0);
	}
	
	return y;
}

/*
	Computes the value of array  P using linear interpolation 
	getPointValueOnDistanceMatrixUsingTriangles
*/
double getValueOfD(double **D,double x,double y)
{
	int X,Y;
	double dx,dy,p0,p1,p2,f;

	X = (int)x;
	Y = (int)y;

	dx = (double)(x-X);
	dy = (double)(y-Y);
	
	if (X >= lines-1 || Y >= cols-1)
	{
		if (X >= lines-1 )
			X = lines-1;
		if (Y >= cols-1)
			Y = cols-1;

		return D[X][Y];
	}

	p0 = D[X][Y];
	p2 = D[X+1][Y+1];
	
	if (dx > dy)
	{
		p1 = D[X+1][Y];
		f = p0 - dx*(p0-p1) - dy*(p1-p2);
	}
	else
	{
		p1 = D[X][Y+1];
		f = p0 + dx*(p2-p1) + dy*(p1-p0);
	}

	return f;
}

/*

	Returns the point of r segment  [x1 y1] + r [x2-x1 y2-y1] 
	r E [0,1]
				xs = x1 + (x2-x1) * q1 / (q1-q2);
*/
double getPointOfLineSegment(double x1,double y1,double f1,double x2,double y2,double f2,double f,double *xs,double *ys)
{
	double r;
	double q1 = (f1-f);
	double q2 = (f2-f);
	static int id =  0;

	*xs = 0;
	*ys = 0;

	if (q1*q2 > 0)
	{
		*xs = (x1+x2)/2;
		*ys = (y1+y2)/2;
		++id;
		printf("i = %d, Error in getPointOfLineSegment q1 = %f q2 = %f \n",id,q1,q2);
		return -1;
	}
	
	if (q1 == 0 && q2 == 0)
		q1 = (double) 0.0001;

	r = q1 / (q1-q2);
	
	*xs = x1 + r*(x2-x1);
	*ys = y1 + r*(y2-y1);

	return r;
}

/*
	Computes the function Q(s) = d(a(s),b(s)) - d(1,a(s))
*/
double getQ_s(double a,double b)
{
	double f1 = getValueOfD(ep->D,a,b);
	double f2 = getValueOfD(ep->D,(double)ep->M-1,a);
	double r;

	r = f1-f2;

	return r;
}

/*
	Returns the roots 
*/
int getQ_s_Roots(double x1,double y1,double x2,double y2,double *xs,double *ys,double *f)
{
	double i1,j1,i2,j2;
	double q1,q2;
	double rx = (double)(ABS(x2-x1) / DISTANCE_E(x2,y2,x1,y1));
	double ry = (double)(ABS(y2-y1) / DISTANCE_E(x2,y2,x1,y1));
	double f1,f2,t;

	for (t = 0; t < (double)DISTANCE_E(x2,y2,x1,y1) / rx; t += (double)0.1)
	{
		i1 = x1 + t*rx;
		j1 = y1 + t*ry;
		i2 = i1 + rx;
		j2 = j1 + ry;

		if (DISTANCE_E2(i2,j2,x1,y1) > DISTANCE_E2(x2,y2,x1,y1))
			break;

		q1 = getQ_s(i1,j1); 
		q2 = getQ_s(i2,j2);
			
		if (q1*q2  <= 0) //riza
		{
			
			f1 = getValueOfD(ep->D, i1,j1);
			f2 = getValueOfD(ep->D, i2,j2);

			*xs = i1 + (i2-i1) * q1 / (q1-q2);
			*ys = j1 + (j2-j1) * q1 / (q1-q2);
			*f = f1 + (f2-f1) * q1 / (q1-q2);
			*f = getValueOfD(ep->D,*xs,*ys);
			
			printf("Roots Computation was done using high level - grid\n");
			
			return 1;
		}
		
	}
	
	q1 = getQ_s(x1,y1); 
	q2 = getQ_s(x2,y2);
	f1 = getValueOfD(ep->D, x1,y1);
	f2 = getValueOfD(ep->D, x2,y2);

	if (q1*q2 <= 0)
	{
		*xs = x1 + (x2-x1) * q1 / (q1-q2);
		*ys = y1 + (y2-y1) * q1 / (q1-q2);
	
		//*f = f1 + (f2-f1)  * q1 / (q1-q2);
		*f = getValueOfD(ep->D,*xs,*ys);
	}

	if (q1*q2  <= 0)
		return 1;
	else
		return 0;
}
/*
	Computes the threshold
*/
void getRealThreshold()
{
	int i,j;
	double max= 0;
	double **D = ep->D;

	for (i = 0; i < lines; ++i)
	for (j = i+1; j < lines; ++j)
	{
		if (D[i][j] > max)
			max = D[i][j];
	}


	ep->Threshold = ep->Threshold*MAX(max,MAX_ERROR_THRESHOLD);

	printf("real ep->Threshold = %f\n",ep->Threshold);
}

void refreshInitMemory(int N,int M,int dim,int InitCace,double ThresholdPercet)
{
	if (InitCace == 1)
	{
		ep = (struct  curveEP *)malloc(sizeof(struct curveEP));
		ep->C = (double **) Alloc2D(M,sizeof(double*),dim,sizeof(double));
		ep->D = (double **) Alloc2D(M,sizeof(double*),M,sizeof(double));
		ep->MinMaxKF = (int **) Alloc2D(M,sizeof(int*),M,sizeof(int));

		ep->M = M;
		ep->N = N;
		ep->Threshold = ThresholdPercet;

		lines = M;
		cols = M;
		ep->n = dim;
		ep->level = 0;

		ep->p = (struct  point2DF **)malloc(N*sizeof(struct point2DF *));
		ep->p[0] = (struct  point2DF *)malloc(M*sizeof(struct point2DF ));

		ep->NumberOfPoints = (int *)calloc(N, sizeof( int));
		ep->NumberOfPoints[0] = M;

		ep->l = (struct lineSegment **)malloc(N*sizeof(struct lineSegment *));
		ep->l[0] = (struct lineSegment *)malloc((M-1)*sizeof(struct lineSegment ));

		ep->NumberOfLines = (int *)calloc(N, sizeof( int));
		ep->NumberOfLines[0] = M-1;

		ep->sol.numSol = 0;

		ep->numLines = (int **) Alloc2D(N,sizeof(int*),M,sizeof(int));
		ep->linesId  = (int ***) Alloc3D(N,sizeof(int**),M,sizeof(int*),1,sizeof(int));
		
		return;
	}

	
	ep->level++;

	//mem init 

	ep->p[ep->level] = (struct  point2DF *)malloc(10*M*sizeof(struct point2DF )); //arxikopoihsh gia points
	ep->NumberOfPoints[ep->level] = 10*M;

	ep->l[ep->level] = (struct lineSegment *)malloc(10*M*sizeof(struct lineSegment ));
	ep->NumberOfLines[ep->level] = 10*M;
}



/*
	  NullPlaneCurves  init (first level)
*/
void initNullPlaneCurves()
{
	int i;

	for (i = 0; i < ep->NumberOfPoints[0]; ++i)
	{
		ep->p[0][i].degree = 0;
	}

	for (i = 0; i < ep->NumberOfPoints[0]; ++i)
	{
		ep->p[0][i].x = (double)i;
		ep->p[0][i].y = 0;
		ep->p[0][i].d = ep->D[0][i];
		ep->p[0][i].fatherLine = 0;
		
		if (i < lines-1)
		{
			ep->l[0][i].apo = i;
			ep->l[0][i].eos = i+1;
            ep->l[0][i].legal = 1;
			ep->p[0][i].line[ep->p[0][i].degree] = i;
			ep->p[0][i+1].line[ep->p[0][i+1].degree] = i;

			ep->numLines[0][i] = 1;
			ep->linesId[0][i][0] = i;

			ep->p[0][i].degree++;
			ep->p[0][i+1].degree++;
		}
	}
}

int simplifyLines(int w1,int w2)
{
    int level = ep->level,apo,eos,Napo,Neos,i;    
    double X[2],Y[2],x[2],y[2],d[2],g,m,mx,my,ok;
    
    
    Napo = ep->l[level][w1].apo;
    Neos = ep->l[level][w1].eos;
     
    apo = ep->l[level][w2].apo;
    eos = ep->l[level][w2].eos;
             
    X[0] = ep->p[level][Napo].x;
    Y[0] = ep->p[level][Napo].y;
    X[1] = ep->p[level][Neos].x;
    Y[1] = ep->p[level][Neos].y;

    x[0] = ep->p[level][apo].x;
    y[0] = ep->p[level][apo].y;
    x[1] = ep->p[level][eos].x;
    y[1] = ep->p[level][eos].y;
        
    for (i = 0; i < 2; ++i)
    {
      d[i] = DISTANCE_E2(x[i],y[i],X[1-i],Y[1-i]);

       if (d[i] == 0)
        {
           m = DISTANCE_E(x[0],y[0],x[1],y[1])*DISTANCE_E(X[0],Y[0],X[1],Y[1]);
           g = (x[1]-x[0])*(X[1]-X[0])+ (y[1]-y[0])*(Y[1]-Y[0]);     
           g = ABS(g)/m;
           mx = MIN(x[1-i],X[i]);
           my = MIN(y[1-i],Y[i]);
           
           ok = 0;   
           if (x[1-i] >= (int)mx && x[1-i] <= (int)mx+1  && y[1-i] >= (int)my && y[1-i] <= (int)my+1)  
           {
              if (X[i] >= (int)mx && X[i] <= (int)mx+1  && Y[i] >= (int)my && Y[i] <= (int)my+1)  
                ok = 1;              
           }
              
           if (ok == 1 && g > 0.995)
           {
                ep->l[level][w2].legal = 0;
                
                if (i == 0)
                {
                  ep->l[level][w1].apo = Napo;
                  ep->l[level][w1].eos = eos;
                }
                else
                {
                  ep->l[level][w1].apo = apo;
                  ep->l[level][w1].eos = Neos;                    
                }
                 
                
                return 1;   
           }
        }
    }
    return 0;
}


/*
   Merges the small line segments to reduce the computational cost   
*/
void mergeLineSegments()
{
   int level = ep->level;    
   int i,j,Napo,Neos,id = ep->NumberOfLines[level],window = 20;
   int ok = 0;
    
    
    for (i = 1; i < id; ++i)
	{
        if (ep->l[level][i].legal == 0)
           continue;
           
        Napo = ep->l[level][i].apo;
        Neos = ep->l[level][i].eos;

        for (j = i-1; j >= MAX(0,i-window); --j)
        {
            if (ep->l[level][j].legal == 0)
               continue;
               
            if (simplifyLines(i,j) == 1)
            {
               ok++;
               break;    
            }
        }
    }
    
    printf("NumberOflines after simplification = %d\n",id-ok);
    
}

/*Retures the average length of e.t*/

double getMeanLengthOfLinesSegments()
{
   int level = ep->level;    
   int i,Napo,Neos,id = ep->NumberOfLines[level];
   double X0,X1,Y0,Y1,d1,md = 0;
    
    
    for (i = 0; i < id; ++i)
	{
        Napo = ep->l[level][i].apo;
        Neos = ep->l[level][i].eos;
        
        X0 = ep->p[level][Napo].x;
        Y0 = ep->p[level][Napo].y;
        X1 = ep->p[level][Neos].x;
        Y1 = ep->p[level][Neos].y;
        d1 = DISTANCE_E(X0,Y0,X1,Y1); 

        md += d1;   
	}   
	
	md = md / id;
	
	printf("average length (%d) = %f\n",level,md);
	return md;
}

//Return the new =id --- Allagh!!!!!!!!!!!
int addNewPoint(int level, int id, double x, double y, double f, int father)
{
	double **D = ep->D;
	double w = x-((int)(x));
//	int i;

	/*for (i = 0; i < id; ++i)
	{
		if (DISTANCE_E2(x,y,ep->p[level][i].x,ep->p[level][i].y) == 0.0)
			return id;
	}*/

	if (y > x)
		return id;

	if (id == ep->NumberOfPoints[level])
	{
		ep->NumberOfPoints[level] = 2*ep->NumberOfPoints[level];
		ep->p[level] = (struct point2DF *)realloc(ep->p[level],ep->NumberOfPoints[level]*sizeof(struct point2DF));
	}

	ep->p[level][id].x = x;
	ep->p[level][id].y = y;
	ep->p[level][id].d = f;
	
	ep->p[level][id].fatherLine = father;
	ep->p[level][id].degree = 0;

	ep->l[level-1][father].numChilds++;

	id++;

	return id;
}


/*
	add a new line segment 
*/

int addNewLineSegment(int level,int id,int apo, int eos)
{
    int i,Napo,Neos,temp;
    double x0,x1,y0,y1;
    double X0,X1,Y0,Y1,d1,d2;
    
    x0 = ep->p[level][apo].x;
    y0 = ep->p[level][apo].y;
    x1 = ep->p[level][eos].x;
    y1 = ep->p[level][eos].y;
    
    d1 = DISTANCE_E2(x0,y0,x1,y1); 
    
    if (d1 == 0)
		return id;
    
    for (i = 0; i < id; ++i)
	{
        Napo = ep->l[level][i].apo;
        Neos = ep->l[level][i].eos;
        
        X0 = ep->p[level][Napo].x;
        Y0 = ep->p[level][Napo].y;
        X1 = ep->p[level][Neos].x;
        Y1 = ep->p[level][Neos].y;
        d1 = DISTANCE_E2(x0,y0,X0,Y0) + DISTANCE_E2(x1,y1,X1,Y1); 
        d2 = DISTANCE_E2(x1,y1,X0,Y0) + DISTANCE_E2(x0,y0,X1,Y1); 
        d1 = MIN(d1,d2);
        
		if (d1 == 0)
			return id;
	}
	
	if (x0 > x1)
	{
       SWAP(apo,eos,temp);
    }
	
	
	if (id == ep->NumberOfLines[ep->level])
	{
		ep->NumberOfLines[level] = 2*ep->NumberOfLines[level];
		ep->l[level] = (struct lineSegment *)realloc(ep->l[level],ep->NumberOfLines[level]*sizeof(struct lineSegment));
	}

	ep->l[level][id].apo = apo;
	ep->l[level][id].eos = eos;
	ep->l[level][id].legal = 1;
	ep->l[level][id].numChilds = 0;

	ep->p[level][apo].line[ep->p[level][apo].degree % MAXIMUMDEGREE] = id;
	ep->p[level][eos].line[ep->p[level][eos].degree % MAXIMUMDEGREE] = id;

	
	ep->p[level][apo].degree++;
	ep->p[level][eos].degree++;

	if (ep->p[level][apo].degree > MAXIMUMDEGREE)
		ep->p[level][apo].degree = MAXIMUMDEGREE;

	if (ep->p[level][eos].degree > MAXIMUMDEGREE)
		ep->p[level][eos].degree = MAXIMUMDEGREE;

	id++;

	return id;
}



/*
	Computes the points of new NullPlaneCurves, 
	from level-1 ---> level

	using triangles  in the search space 
*/
void getNextLevelOfNullPlaneCurves()
{
	int i,i1,apo,eos,xmin,lineId = 0, id = 0,id0,ok;
	int level = ep->level,x;
	double X[2],F[2],R[2], dc,**D = ep->D,j1,j2,xs,ys;
	static double **V = 0;
	int kUp[2][100];
    int kDown[2][100];

	if (V == 0)
		V = (double **) Alloc2D(lines,sizeof(double*),4,sizeof(double));


	for (i = 0; i < ep->NumberOfLines[level-1]; ++i)
	{
        if (ep->l[level-1][i].legal == 0)
           continue;
           
		apo = ep->l[level-1][i].apo;
		eos = ep->l[level-1][i].eos;
		
		X[0] = ep->p[level-1][apo].x;
		F[0] = ep->p[level-1][apo].d;
		R[0] = X[0] - ((int)X[0]);

		X[1] = ep->p[level-1][eos].x;
		F[1] = ep->p[level-1][eos].d;
		R[1] = X[1] - ((int)X[1]);
		
		if (R[1] == 0 && X[1]-X[0] > 0)
		{
           R[1] = 1;        
        }
		   
//		printf("X[0] = %f  X[1] = %f\n",X[0],X[1]);

		xmin = (int)MIN(X[0],X[1]);
		

		for (x = (int)xmin; x < lines-1; ++x)
		{
			for (i1 = 0; i1 < 2; ++i1) //basikoi a3ones
			{
				j1 = x;
		  		j2 = X[i1];
				dc = getValueOfD(D,j1,j2);
				V[x][i1] = F[i1] - dc; 
			}

			for (i1 = 0; i1 < 2; ++i1) //diagwnio
			{
				j1 = x + R[i1];
		  		j2 = X[i1];
				dc = getValueOfD(D,j1,j2);
				V[x][i1+2] = F[i1] - dc; 
			}
		}
		ok = id;
		for (x = (int)xmin; x < lines-2; ++x)
		{
			kUp[0][0] = 0;
			kUp[1][0] = 0;
   		    kDown[0][0] = 0;
			kDown[1][0] = 0;

            //panw trigwno			
			for (i1 = 0; i1 < 2; ++i1)
			{
			   if (V[x][i1] == 0 && V[x][i1+2] == 0)
			   {
                    id0 = id;                            
                    xs = x;
                    ys = X[i1];
                    id = addNewPoint(level, id, xs, ys, getValueOfD(D,xs,ys),i);   
                    if (id > id0)
 				    {
						kUp[i1][kUp[i1][0]+1] = id-1;
					    kUp[i1][0] += id-id0;                           
				    }
               }
               else if (V[x][i1]*V[x][i1+2] <= 0)
               {
                    id0 = id;
                    xs = x + R[i1]*V[x][i1] / (V[x][i1] - V[x][i1+2]);
                    ys = X[i1];
                    id = addNewPoint(level, id, xs, ys, getValueOfD(D,xs,ys),i);   
                    if (id > id0)
 				    {
						kUp[i1][kUp[i1][0]+1] = id-1;
                        kUp[i1][0] += id-id0;						
				    }
                 }
            }
            
		   if (V[x][0]*V[x][1] < 0)
		   { 
                 id0 = id;
                 xs = x;
                 ys = X[0] + (X[1]-X[0])*V[x][0] / (V[x][0] - V[x][1]);
                 id = addNewPoint(level, id, xs, ys, getValueOfD(D,xs,ys),i); 
//                 printf("Diag: x = %d X[0] = %f xs = %f ys = %f\n",x,X[0],xs,ys);
  //               getchar();     
                 for (i1 = 0; i1 < 2; ++i1)
                 {     
                    if (id > id0)
 				    {
                        if (kUp[i1][0] == 0 && kUp[1-i1][0] > 0)
                        {
						   kUp[i1][kUp[i1][0]+1] = id-1;
         				   kUp[i1][0] += id-id0;
                        }			     					     	
                        else if (kUp[i1][0] == 0 && kUp[1-i1][0] == 0)
                        {
						   kUp[i1][kUp[i1][0]+1] = id-1;
         				   kUp[i1][0] += id-id0;
         				   break;
                        }			     					                            
				    }
                 }
           }            
         	ok = id-ok;
         	if (ok > 0)
         	{
              // printf("Points on up Triangle = %d (%d) xs = %f ys = %f\n",ok,x,xs,ys);      
             }
             ok = id;
            //katw trigwno
            for (i1 = 0; i1 < 2; ++i1)
			{
			   if (V[x+1][i1] == 0 && V[x][i1+2] == 0)
			   {
                    id0 = id;                              
                    xs = x+1;
                    ys = X[i1];
                    id = addNewPoint(level, id, xs, ys, getValueOfD(D,xs,ys),i);   
                    if (id > id0)
 				    {
						kDown[i1][kDown[i1][0]+1] = id-1;
						kDown[i1][0] += id-id0;                           
				    }
               }
               else if (V[x+1][i1]*V[x][i1+2] <= 0)
               {
                    id0 = id;
                    xs = x + R[i1] + (1-R[i1])*V[x][i1+2] / (V[x][i1+2] - V[x+1][i1]);
                    ys = X[i1];
                    id = addNewPoint(level, id, xs, ys, getValueOfD(D,xs,ys),i);   
                    if (id > id0)
 				    {
						kDown[i1][kDown[i1][0]+1] = id-1;
						kDown[i1][0] += id-id0;
				    }
                 }
            }

		   if (V[x+1][0]*V[x+1][1] < 0)
		   { 
                 id0 = id;
                 xs = x+1;
                 ys = X[0] + (X[1]-X[0])*V[x+1][0] / (V[x+1][0] - V[x+1][1]);
                 id = addNewPoint(level, id, xs, ys, getValueOfD(D,xs,ys),i); 
//                 printf("Diag: x = %d X[0] = %f xs = %f ys = %f\n",x,X[0],xs,ys);
  //               getchar();     
                 for (i1 = 0; i1 < 2; ++i1)
                 {     
                    if (id > id0)
 				    {
                        if (kDown[i1][0] == 0 && kDown[1-i1][0] > 0)
                        {   
                        	kDown[i1][kDown[i1][0]+1] = id-1;
                 			kDown[i1][0] += id-id0;			     					 
                        }
                        else if (kDown[i1][0] == 0 && kDown[1-i1][0] == 0)    	
                        {
                            kDown[i1][kDown[i1][0]+1] = id-1;
                 			kDown[i1][0] += id-id0;			 
                            break; 
                        }
				    }
                 }
           }                        
            ok = id-ok;
         	if (ok > 0)
         	{
            //   printf("Points on down Triangle = %d (%d) xs = %f ys = %f(%d - %d) (%d - %d)\n",ok,x,xs,ys,kUp[0][0],kUp[1][0],kDown[0][0],kDown[1][0]);      
             } 
             ok = id;            
//diagwnios
		   if (V[x][2]*V[x][3] < 0)
		   { 
                 id0 = id;
                 xs = x + R[0] + (X[1]-X[0])*V[x][2] / (V[x][2] - V[x][3]);
                 ys = X[0] + (X[1]-X[0])*V[x][2] / (V[x][2] - V[x][3]);
                 id = addNewPoint(level, id, xs, ys, getValueOfD(D,xs,ys),i); 
//                 printf("Diag: x = %d X[0] = %f xs = %f ys = %f\n",x,X[0],xs,ys);
  //               getchar();     
        
                 if (id > id0)
                 {
                      for (i1 = 0; i1 < 2; ++i1)
                      {     
                        if (kDown[i1][0] == 0 && kDown[1-i1][0] > 0)
                        {   
                        	kDown[i1][kDown[i1][0]+1] = id-1;
                 			kDown[i1][0] += id-id0;			     					 
                        }    	
                        if (kUp[i1][0] == 0 && kUp[1-i1][0] > 0)
                        {
						   kUp[i1][kUp[i1][0]+1] = id-1;
         				   kUp[i1][0] += id-id0;
                        }			     					     	
				    }
                 }
           }
                  
                     ok = id-ok;
         	if (ok > 0)
         	{
              // printf("Points DT = %d (%d) xs = %f ys = %f (%d - %d) (%d - %d)\n",ok,x,xs,ys,kUp[0][0],kUp[1][0],kDown[0][0],kDown[1][0]);      
            //   getchar();
             } 
             ok = id;


       	    for (i1 = 0; i1 < MIN(kUp[0][0],kUp[1][0]); ++i1)
            {
				lineId = addNewLineSegment(level,lineId,kUp[0][1+i1],kUp[1][1+i1]);
			}                        
       	    for (i1 = 0; i1 < MIN(kDown[0][0],kDown[1][0]); ++i1)
			{
				lineId = addNewLineSegment(level,lineId,kDown[0][1+i1],kDown[1][1+i1]);
			}                        

		}
	
	}

	ep->NumberOfPoints[level] = id;
	ep->NumberOfLines[level] = lineId;
	printf("***Level %d := %d points %d lines\n",level,id,lineId);
	mergeLineSegments();
	getMeanLengthOfLinesSegments();	
//	getchar();
	
}

/*
	Computes Thresholds kai corrects the Distance Matrix
*/
double getGlobalThreshold_CorrectDistanceMatrix()
{
	double f = 0.0,c,t,t1,t2;
	int i,j,w = 0;
	int i1,i2,k,added = 0;

	for (i = 1; i < ep->M-1; ++i)
	for (j = i+2; j < ep->M-2; ++j)
	{
		for (i1 = -1; i1 <= 1; ++i1)
		for (i2 = -1; i2 <= 1; ++i2)
		{
			f += ABS(ep->D[i][j] - ep->D[i+i1][j+i2]);
			w++;
		}
	}

	c = ((double)f/w);
	   return c;
/*

added little noise in the case of same matrix points
*/

	if (ADD_TIME_TO_METRIC == 1)
	{
		printf("ADD_TIME_TO_METRIC \n");
		for (i = 0; i < lines; ++i)
		for (j = i; j < lines; ++j)
		{
			ep->D[i][j] += (double)(0.1*c*ABS(i-j))/(double)lines;
			ep->D[j][i] += (double)(0.1*c*ABS(i-j))/(double)lines;
		}
	}

	
	for (i = 1; i < ep->M; ++i)
	for (j = i+1; j < ep->M; ++j)
	{
		for (k = j+1; k < ep->M; ++k)
		{
			if (ep->D[i][j] == ep->D[i][k])
			{
				t1 = ((double)rand());
				t2 = ((double)rand()/RAND_MAX);
				t = t1+t2;
				t = ((double)t/RAND_MAX);

				t = 0.00001*(0.5*ep->D[i][j]+0.5*c)*t;
				ep->D[i][j] += c;
				ep->D[j][i] += c;
				added++;
			}
		}

		for (k = 1; k < i; ++k)
		{
			if (ep->D[i][j] == ep->D[k][j])
			{
				t1 = ((double)rand());
				t2 = ((double)rand()/RAND_MAX);
				t = t1+t2;
				t = ((double)t/RAND_MAX);

				t = 0.00001*(0.5*ep->D[i][j]+0.5*c)*t;
				ep->D[i][j] += c;
				ep->D[j][i] += c;
				added++;
			}
		}
	}

	printf("Number of points with the same D = %d \n",added);

	if (ep->N <= 0)
	{
		getRealThreshold();
	}

	return c;
}



/*
	Computes Thresholds kai corrects the Distance Matrix
*/
double getGlobalThreshold_CorrectDistanceMatrix2()
{
	double f = 0.0,c = INF,t,t1,t2;
	int i,j,w = 0;
	int i1,i2,k,added = 0;


	for (i = 1; i < lines; ++i)
	for (j = i+1; j < lines; ++j)
	{
        if (ep->D[i][j] > 0 && ep->D[i][j] < c)
           c = ep->D[i][j];
	}


/*
added little noise in the case of same matrix points
*/

	if (ADD_TIME_TO_METRIC == 1)
	{
		printf("ADD_TIME_TO_METRIC \n");
		for (i = 0; i < lines; ++i)
		for (j = i; j < lines; ++j)
		{
			ep->D[i][j] += (double)(0.1*c*ABS(i-j))/(double)lines;
			ep->D[j][i] += (double)(0.1*c*ABS(i-j))/(double)lines;
		}
	}

	
	for (i = 1; i < ep->M; ++i)
	for (j = i+1; j < ep->M; ++j)
	{
		for (k = j+1; k < ep->M; ++k)
		{
			if (ep->D[i][j] == ep->D[i][k])
			{
				t1 = ((double)rand());
				t2 = ((double)rand()/RAND_MAX);
				t = t1+t2;
				t = ((double)t/RAND_MAX);

				t = abs(t*0.01*c);
				ep->D[i][j] += t;
				ep->D[j][i] += t;
				added++;
			}
		}
	}

	printf("Ari8mos shmeiwn me idio D = %d \n",added);

	if (ep->N <= 0)
	{
		getRealThreshold();
	}

	return c;
}



double getMaxCostThreshold(double thr, double x1,double y1,double x2,double y2)
{
	double **D = ep->D;
	double f;

	f = gThresh + thr*(ABS(getValueOfD(D,x1,y1) -  getValueOfD(D,x1,y1+1))+ ABS(getValueOfD(D,x1,y1) -  getValueOfD(D,x1+1,y1)) + ABS(getValueOfD(D,x2,y2) -  getValueOfD(D,x2,y2+1))+ABS(getValueOfD(D,x2,y2) -  getValueOfD(D,x2+1,y2)));

	return f;
}





/*
free mem
*/
void freeSolution()
{
	int i;
	
	for (i = 0; i < ep->sol.numSol; ++i)
		free(ep->sol.p[i]);

	free(ep->sol.p);
	free(ep->sol.r);
	ep->sol.numSol = 0;
}

/*
	computes the solutions recursively  
*/


double getSolutions(int functionChoice)
{
	int i,j,father1,father2,father;
	int pointA,pointB;
	int N = ep->N,seg;
    double q1,q2,f,x1,y1,x2,y2,f1,f2,xs,ys,res = INF;
	int ok = 0,riza;


	for (i = 0; i < ep->NumberOfLines[N-2]; ++i)
	{
        if (ep->l[N-2][i].legal == 0)
           continue;
           
		pointA = ep->l[N-2][i].apo;
		pointB = ep->l[N-2][i].eos;
		
        x1 = ep->p[N-2][pointA].x;
		y1 = ep->p[N-2][pointA].y;

		x2 = ep->p[N-2][pointB].x;
		y2 = ep->p[N-2][pointB].y;

		f1 = ep->p[N-2][pointA].d;
		f2 = ep->p[N-2][pointB].d;
			

		q1 = getQ_s(x1,y1); 
		q2 = getQ_s(x2,y2);
		riza = 0;

		if (functionChoice == 1 && q1*q2 <= 0)
		{
			riza = 1;
			xs = x1 + (x2-x1) * q1 / (q1-q2);
			ys = y1 + (y2-y1) * q1 / (q1-q2);
			f = f1 + (f2-f1) * q1 / (q1-q2);   //Error = 0.000219
			//printf("!!! f = %f f1 = %f f2 = %f FF = %f x1 = %f x2 = %f\n",f,f1,f2,ep->p[0][10].d,x1,x2);
		}
		else if (functionChoice == 2 && getQ_s_Roots( x1, y1, x2, y2,&xs,&ys,&f) == 1)
		{
			riza = 1;
		}
		
		if (riza == 1)
		{
			ok = 1;
			seg = i;
		
			ep->sol.numSol++;

			if (ep->sol.numSol == 1)
			{
				ep->sol.r  = (double *)malloc(sizeof(double));

				ep->sol.p = (double **)malloc(sizeof(double *));
				ep->sol.p[0] = (double *)malloc((ep->N+1)*sizeof(double ));

			}
			else
			{
				ep->sol.r  = (double *)realloc(ep->sol.r,ep->sol.numSol*sizeof(double));

				ep->sol.p = (double **)realloc(ep->sol.p,ep->sol.numSol*sizeof(double *));
				ep->sol.p[ep->sol.numSol-1] = (double *)malloc((ep->N+1)*sizeof(double ));
			}

			ep->sol.r[ep->sol.numSol-1] = f;
			ep->sol.p[ep->sol.numSol-1][N] = (double)(ep->M - 1);
			ep->sol.p[ep->sol.numSol-1][N-1] = xs;
			ep->sol.p[ep->sol.numSol-1][N-2] = ys;
			ep->sol.p[ep->sol.numSol-1][0] = (double)0;
			father = i;
	
			for (j = N-3; j > 0; --j)
			{

				father1 = ep->p[j+1][ep->l[j+1][father].apo].fatherLine;
				father2 = ep->p[j+1][ep->l[j+1][father].eos].fatherLine;

				father = father1;
				
				pointA = ep->l[j][father].apo;
				pointB = ep->l[j][father].eos;

				f1 = ep->p[j][pointA].d;
				f2 = ep->p[j][pointB].d;

				if ((f - f1)*(f-f2) > 0)
				{
					father = father2;
				
					pointA = ep->l[j][father].apo;
					pointB = ep->l[j][father].eos;

					f1 = ep->p[j][pointA].d;
					f2 = ep->p[j][pointB].d;
				}
				
				x1 = ep->p[j][pointA].x;
				y1 = ep->p[j][pointA].y;
				x2 = ep->p[j][pointB].x;
				y2 = ep->p[j][pointB].y;


				getPointOfLineSegment( x1, y1, f1, x2, y2,  f2, f, &xs, &ys);

				ep->sol.p[ep->sol.numSol-1][j] = ys;
			}
		}
	}

	

	if (functionChoice == 1 && ok == 0)
	{
		printf("Running\n");
		getSolutions(2);
		
	}

	for (i = 0; i < ep->sol.numSol; ++i)
	{
		res = MIN(res,ep->sol.r[i]);
	}

	printf("# Lusewn = %d res = %f\n\n\n", ep->sol.numSol,res);
	return res;
}


/*
These files will be read by the matlab - plotNullPlaneCurves
*/
void printDataForMatlab()
{
	int i,j;
	FILE *fp1;
	double x0,y0,x1,y1,max;
	int N = ep->N;
	char command[1000];
	int *index;

	
	fp1 = fopen("curve.txt", "w");
	for (i = 0; i < lines; ++i)
	{
		for (j = 0; j < ep->n; ++j)
		{
			fprintf(fp1,"%f ",ep->C[i][j]);
		}
		fprintf(fp1,"\n");
	}
	fclose(fp1);
	
	fp1 = fopen("matrix.txt", "w");
	for (i = 0; i < lines; ++i)
	{
		for (j = 0; j < cols; ++j)
		{
			fprintf(fp1,"%f ",ep->D[i][j]);
		}
		fprintf(fp1,"\n");
	}
	fclose(fp1);
//LINES SEGMENTS 
	fp1 = fopen("segments.txt", "w");

	for (i = 1; i <= ep->level; ++i)
	{
		fprintf(fp1,"0 0 0 0\n");
		for (j = 0; j < ep->NumberOfLines[i]; ++j)
		{
            if (ep->l[i][j].legal == 0)
               continue;
			x0 = ep->p[i][ep->l[i][j].apo].x;
			y0 = ep->p[i][ep->l[i][j].apo].y;
			x1 = ep->p[i][ep->l[i][j].eos].x;
			y1 = ep->p[i][ep->l[i][j].eos].y;
			fprintf(fp1,"%f %f %f %f\n",x0,y0,x1,y1);
		}
	}

	fclose(fp1);

	

	//POINTS

	fp1 = fopen("points.txt", "w");

	for (i = 0; i <= ep->level; ++i)
	{
		fprintf(fp1,"0 0 0 0\n");
		for (j = 0; j < ep->NumberOfPoints[i]; ++j)
		{
			x0 = ep->p[i][j].x;
			y0 = ep->p[i][j].y;
			fprintf(fp1,"%f %f %d %f\n",x0,y0,ep->p[i][j].degree,getValueOfD(ep->D,x0,y0));
		}
	}
	fclose(fp1);

	//Sorting solutions

	index = (int *)malloc((1+ep->sol.numSol)*sizeof(int));
	for (i = 0; i < ep->sol.numSol; ++i)
		index[i] = i;

    for (i = 0; i < ep->sol.numSol; ++i)
	{
		max = ep->sol.r[index[i]];
		for (j = i+1; j < ep->sol.numSol; ++j)
		{
			if (max < ep->sol.r[j])
			{
				index[i] = j;
				index[j] = i;
				max = ep->sol.r[j];
				ep->sol.r[j] = ep->sol.r[i];
				ep->sol.r[i] = max;
			}
		}
	}

	fp1 = fopen("solutions.txt", "w");

		fprintf(fp1,"%f %f \n",(double)ep->level, (double)ep->sol.numSol);

	for (i = 0; i < ep->sol.numSol; ++i)
	{
		for (j = 0; j <= ep->level+2; ++j)
		{
			fprintf(fp1,"%f ",ep->sol.p[index[i]][j]);
		}
		fprintf(fp1,"\n");
	}
	fclose(fp1);

	for (i = 0; i < ep->sol.numSol; ++i)
	{
		printf("r[%d] = %f \n",i,ep->sol.r[i]);
	}

	if (metricChoise == MIN_MAX_KEYFRAMES_METRIC)
	{
		printfMinMaxSolutions(index);
	}
	
	free(index);

//backup
	printf("\n\nResults are copied to %s \n",ep->outDir);

	sprintf(command,"copy curve.txt %s%s_matrix.txt\n",ep->outDir,ep->outName);
	printf(">>%s \n",command);
	system(command);
	
	sprintf(command,"copy matrix.txt %s%s_matrix.txt\n",ep->outDir,ep->outName);
	printf(">>%s \n",command);
	system(command);
	
	sprintf(command,"copy solutions.txt %s%s_solutions.txt\n",ep->outDir,ep->outName);
	printf(">>%s \n",command);
	system(command);

	sprintf(command,"copy points.txt %s%s_points.txt\n",ep->outDir,ep->outName);
	printf(">>%s \n",command);
	system(command);

//	sprintf(command,"cd ..%c..%cMatlabFiles\n",92,92);
//	system(command);
//	printf(">>%s \n",command);

}

/*
	from the curve and its dimensions we compute  M
*/
int getNumberMFromCurve(char *name,int dim)
{
	int M = 0;
	FILE *fp1;
	float f;
	
	fp1 = fopen(name, "r");
	while (fscanf(fp1,"%f",&f) != EOF)
	{
		M++;	
	}
	fclose(fp1);

	return (int)(M / dim);
}


/*
	we compute  M based on D, 
*/
int getNumberM(char *name)
{
	int M = 0;
	FILE *fp1;
	float f;
	
	fp1 = fopen(name, "r");

	while (fscanf(fp1,"%f",&f) != EOF)
	{
		M++;	
	}
	fclose(fp1);

	return (int)(sqrt(M));
}




/* we compute DistanceMatrix D  from curve */

void getDistanceMatrixFromCurve(char *name, int choice)
{
	int i,j,k,k1;
	FILE *fp1;
	double max,d1,d2,d3,fd,fd1,fd2,fd3;
	float f;
	int dim = ep->n;
		
	fp1 = fopen(name, "r");

	for (i = 0; i < lines; ++i)
	for (k = 0; k < dim; ++k)
	{
		fscanf(fp1,"%f ",&f);
		ep->C[i][k] = f;
	}
	fclose(fp1);

	if (choice == L1_METRIC)
	{
		printf("L1_METRIC is used\n");
		for (i = 0; i < lines; ++i)
		for (j = i; j < lines; ++j)
		{
			fd = 0;
			for (k = 0; k < dim; ++k)
			{
				fd += ABS(ep->C[i][k] - ep->C[j][k]); 
			}
			ep->D[i][j] = fd;
			ep->D[j][i] = fd;
		}
	}
	else if (choice == L2_METRIC)
	{
		printf("L2_METRIC is used\n");
		for (i = 0; i < lines; ++i)
		for (j = i; j < lines; ++j)
		{
			fd = 0;
			for (k = 0; k < dim; ++k)
			{
				fd += (ep->C[i][k]-ep->C[j][k])*(ep->C[i][k]-ep->C[j][k]);
			}
			fd = (double)(sqrt(fd));
			ep->D[i][j] = fd;
			ep->D[j][i] = fd;
		}
	}
	else if (choice == LINF_METRIC)
	{
		printf("LINF_METRIC is used\n");
		for (i = 0; i < lines; ++i)
		for (j = i; j < lines; ++j)
		{
			max = 0;
			for (k = 0; k < dim; ++k)
			{
				fd = ABS(ep->C[i][k]-ep->C[j][k]); 
				if (fd > max)
					max = fd;
			}
			fd = max;
			ep->D[i][j] = fd;
			ep->D[j][i] = fd;
		}
	}
	else if (choice == POLYGONERROR_LISE_METRIC)
	{
		printf("POLYGONERROR_LISE_METRIC is used\n");
		for (i = 0; i < lines; ++i)
		for (j = i; j < lines; ++j)
		{
			d2 = 0;
			for (k = i+1; k <= j-1; ++k)
			{
				d3 = getDistancePointLineSegment(ep->C[k],ep->C[i],ep->C[j],ep->n);
				d2 += d3*d3;
			}
			ep->D[i][j] = (double)((((d2))));
			ep->D[j][i] = ep->D[i][j];
		}
	}
	else if (choice == POLYGONERROR_ToleranceZone_METRIC)
	{
		printf("POLYGONERROR_ToleranceZone_METRIC is used\n");
		for (i = 0; i < lines; ++i)
		for (j = i; j < lines; ++j)
		{
			d2 = 0;
			for (k = i+1; k <= j-1; ++k)
			{
				d3 = getDistancePointLineSegment(ep->C[k],ep->C[i],ep->C[j],ep->n);
				if (d3 > d2)
					d2 = d3;
			}
			ep->D[i][j] = (double)((d2));
			ep->D[j][i] = ep->D[i][j];
		}
	}
	else if (choice == LABS_METRIC)
	{
		printf("LABS_METRIC is used\n");
		for (i = 0; i < lines; ++i)
		for (j = i; j < lines; ++j)
		{
			d2 = 0;
			d3 = 0;
			for (k = 0; k < dim; ++k)
			{
				d1 = ep->C[i][k]; 
				d1 = d1*d1;
				d2 += d1; 
				d1 = ep->C[j][k]; 
				d1 = d1*d1;
				d3 += d1; 
			}
			fd = (double)ABS(sqrt(d2)-sqrt(d3));
			ep->D[i][j] = fd;
			ep->D[j][i] = fd;
		}
	}
	else if (choice == LKEYFRAMES_METRIC)
	{
		printf("LKEYFRAMES_METRIC is used\n");
		for (i = 0; i < lines; ++i)
		for (j = i; j < lines; ++j)
		{
			fd = 0;
			for (k = i+1; k <= j-1; ++k)
			{
				d2 = 0;
				d3 = 0;
				for (k1 = 0; k1 < dim; ++k1)
				{
					d2 += (ep->C[i][k1]-ep->C[k][k1])*(ep->C[i][k1]-ep->C[k][k1]);
					d3 += (ep->C[j][k1]-ep->C[k][k1])*(ep->C[j][k1]-ep->C[k][k1]);
				}
				d2 = sqrt(d2);
				d3 = sqrt(d3);
				fd += MAX(d2,d3);
			}
			 
			ep->D[i][j] = fd;
			ep->D[j][i] = fd;
		}
	}
	else if (choice == L2KEYFRAMES_METRIC)
	{
	    printf("L2_METRIC/Key Frames is used\n");
		for (i = 0; i < lines; ++i)
		for (j = i; j < lines; ++j)
		{
			fd1 = 0;
			for (k = 0; k < dim/3; ++k)
			{
				fd1 += (ep->C[i][k]-ep->C[j][k])*(ep->C[i][k]-ep->C[j][k]);
			}
			fd2 = 0;
			for (k = dim/3; k < 2*dim/3; ++k)
			{
				fd2 += (ep->C[i][k]-ep->C[j][k])*(ep->C[i][k]-ep->C[j][k]);
			}

			fd3 = 0;
			for (k = 2*dim/3; k < dim; ++k)
			{
				fd3 += (ep->C[i][k]-ep->C[j][k])*(ep->C[i][k]-ep->C[j][k]);
			}

			fd = (double)(sqrt(fd1))+(double)(sqrt(fd2))+(double)(sqrt(fd3));
			ep->D[i][j] = fd;
			ep->D[j][i] = fd;
		}
    }

	gThresh = getGlobalThreshold_CorrectDistanceMatrix2();

	if (ADD_EXP_TO_METRIC == 1)
	{
		d3 = (double)(0.001*gThresh);
		printf("ADD_EXP_TO_METRIC \n");
		for (i = 0; i < lines; ++i)
		for (j = i; j < lines; ++j)
		{
			d1 = (ABS(i-j))/(double)(0.1*lines);
			d2 = (double)(d3*d1*exp(1-d1));
			ep->D[i][j] += d2;
			ep->D[j][i] += d2;
		}
	}


	printf("READING Curve M = %d N = %d dim = %d\n", ep->M,ep->N,ep->n);
	printf("lines = %d cols = %d dist(0,1) = %f \n",lines,cols,ep->D[0][1]);
}


/* readDistanceMatrix D from a file */

void readDistanceMatrix(char *name)
{	
	int i,j;
	FILE *fp1;
	float f;
		
	fp1 = fopen(name, "r");

	for (i = 0; i < lines; ++i)
	for (j = 0; j < cols; ++j)
	{
		fscanf(fp1,"%f",&f);
		ep->D[i][j] = f;
	}

	printf("READING D M = %d N = %d dim = %d\n", ep->M,ep->N,ep->n);
	printf("lines = %d cols = %d\n",lines,cols);
	fclose(fp1);
/*
		Elegxos gia swsta dedomena
*/

	for (i = 0; i < lines; ++i)
	{
			
		for (j = 0; j < cols; ++j)
		{
		//	printf("%3.3f ",ep->D[i][j]);
			if (ep->D[i][j] != ep->D[j][i])
				printf("Error input (symmetry): D[%d][%d] = %f and D[%d][%d] = %f \n",i,j,ep->D[i][j],j,i,ep->D[j][i]);
			if (i != j && ep->D[i][j] == 0)
				printf("Error input : D[%d][%d] = %f \n",i,j,ep->D[i][j]);

		}
//			printf("\n");
	}
}


void mainAlgoOfEP()
{
	int i;

	gThresh = getGlobalThreshold_CorrectDistanceMatrix2();
	printf("GT = %f\n",gThresh);
	
	initNullPlaneCurves();

	printf("Init ok\n");
	for (i = 1; i < ep->N-1; ++i)
	{
		refreshInitMemory(ep->N,ep->M,ep->n,0,0);
		printf("Refresh %d ok\n",i);
		getNextLevelOfNullPlaneCurves();
	}
	getSolutions(1);
	printDataForMatlab();
	
	printf("Running runApprEQP algo\n");
	runApprEQP();
}

double getConstantParameter()
{
	int i;
	double s =0;

	for (i= 1; i <= ep->M-1; ++i)
	{
		s += ep->D[i-1][i];
	}

	return s;

}

void mainAlgoOfEPThresh()
{
	int i,N = ep->N,maxi;
	double res;
	double W[1000],c = getConstantParameter();

	gThresh = getGlobalThreshold_CorrectDistanceMatrix2();
	printf("GT = %f\n",gThresh);
	
	initNullPlaneCurves();

	printf("Init ok ep->Threshold = %f\n",ep->Threshold);
	for (i = 1; i < ep->N-1; ++i)
	{
		refreshInitMemory(ep->N,ep->M,ep->n,0,0);
		printf("Refresh %d ok\n",i);
		getNextLevelOfNullPlaneCurves();

		ep->N = i+1;
		res = getSolutions(1);
		ep->N = N;

		W[i] = res;
		printf("res(%d) = %f \n",i,res);
		maxi = i;

		if (res < ep->Threshold)
		{
			printf("**** Minimun Error = %f, N = %d\n",res,i+1);
			printDataForMatlab();
			break;
		}
		freeSolution();
	}

	W[0] =  ep->D[0][ep->M-1];

	printf("C = [");
	for (i = 0; i <= maxi; ++i)
		printf(" %3.3f ",c - (i+1)*W[i]);
	printf(" ]\n");


	printf("\n\nCd = [");
	for (i = 1; i <= maxi-1; ++i)
		printf(" %3.3f ", -(i+2)*W[i+1]-i*W[i-1] + 2*(i+1)*W[i]);
	printf(" ]\n");


	printf("\n\n (polygonal case LISE - KeyFramesDistortion\n\n\nC = [");
	for (i = 0; i <= maxi; ++i)
		printf(" %3.3f ", (i+1)*W[i]);
	printf(" ]\n");


	printf("\n\nCd = [");
	for (i = 1; i <= maxi-1; ++i)
		printf(" %3.3f ", ((i+2)*W[i+1]+i*W[i-1] - 2*(i+1)*W[i]) / 1);
	printf(" ]\n");

	printf("\n\n (polygonal case Tol.Zone\n\n\nC = [");
	for (i = 0; i <= maxi; ++i)
		printf(" %3.3f ", W[i]);
	printf(" ]\n");


	printf("\n\nCd = [");
	for (i = 1; i <= maxi-1; ++i)
		printf(" %3.3f ", (W[i+1]+W[i-1] - 2*W[i]) / 1);
	printf(" ]\n");

	//printDataForMatlab();
}


int main(int argc, char *argv[])
{
	int M;

	if (argc !=  5 && argc !=  6 && argc != 7 && argc !=  8) 
	{
		fprintf(stdout,"READDISTANCEMATRIX : %s txtFileDistanceMatrix N outDir outname\n",argv[0]);
		fprintf(stdout,"READCURVE : %s txtFileCurve metricChoise CurveDimension N outDir outname\n\n",argv[0]);

		fprintf(stdout,"READDISTANCEMATRIX : %s THRESHOLD txtFileDistanceMatrix ThresholdPercentage  outDir outname\n",argv[0]);
		fprintf(stdout,"READCURVE : %s THRESHOLD txtFileCurve metricChoise CurveDimension ThresholdPercentage  outDir outname\n\n",argv[0]);

		fprintf(stdout,"metricChoise 1: L_1 norm\n");
		fprintf(stdout,"metricChoise 2: L_2 norm\n");
		fprintf(stdout,"metricChoise 3: L_INF norm\n");
		fprintf(stdout,"metricChoise 4: P_LISE norm\n");
		fprintf(stdout,"metricChoise 5: P_Tolerance zone norm\n");
		fprintf(stdout,"metricChoise 6: LABS_METRIC zone norm\n");
		fprintf(stdout,"metricChoise 7: KeyFrames_METRIC zone norm\n");
		fprintf(stdout,"metricChoise 8: MIN-MAX KeyFrames_METRIC \n");
		fprintf(stdout,"metricChoise 9: L_2 KeyFrames_METRIC \n");

		exit(1);
    }

	if (argc == 5)//READDISTANCEMATRIX - N
	{
		printf("EP!!!Using data from  Distance Matrix!!!\n");
		M = getNumberM(argv[1]);
		refreshInitMemory(atoi(argv[2]),M,1,1,0.0);
		readDistanceMatrix(argv[1]);
		sprintf(ep->outDir,"%s",argv[3]);
		sprintf(ep->outName,"%s",argv[4]);
		mainAlgoOfEP();
	}

    if (argc == 7 && atoi(argv[2]) == MIN_MAX_KEYFRAMES_METRIC)//READ CURVE - Using MIN - MAX Algorithm
    {
 //%s txtFileCurve metricChoise Curvedimension N  outDir outname
//0      1			2				3			4  5      6
   		printf("EP-MINMAX!!!Using data from  Curve !!!\n");  
		metricChoise = MIN_MAX_KEYFRAMES_METRIC;                
  		M = getNumberMFromCurve(argv[1],atoi(argv[3]));
		refreshInitMemory(atoi(argv[4])+1,M,atoi(argv[3]),1,0);
		getDistanceMatrixFromCurveUsingMinMax(argv[1],atoi(argv[2]));

		sprintf(ep->outDir,"%s",argv[5]);
		sprintf(ep->outName,"%s",argv[6]);
		mainAlgoOfEP();

    }
	else if (argc == 7)//READ CURVE - N
	{
//%s txtFileCurve metricChoise Curvedimension N  outDir outname
//0      1			2				3			4  5      6
		printf("EP!!!Using data from  Curve!!!\n");
		M = getNumberMFromCurve(argv[1],atoi(argv[3]));
		refreshInitMemory(atoi(argv[4]),M,atoi(argv[3]),1,0);
		getDistanceMatrixFromCurve(argv[1],atoi(argv[2]));

		sprintf(ep->outDir,"%s",argv[5]);
		sprintf(ep->outName,"%s",argv[6]);
		mainAlgoOfEP();
	}

	if (argc == 6)//READDISTANCEMATRIX - Thresholds
	{//READDISTANCEMATRIX : %s THRESHOLD txtFileDistanceMatrix ThresholdPercentage  outDir outname
		printf("READDISTANCEMATRIX EP!!!Using data from Distance Matrix and PreDefined Threshold!!!\n");
		M = getNumberM(argv[2]);
		refreshInitMemory(100,M,1,1,0.0);
		ep->Threshold = atof(argv[1]);
		
		readDistanceMatrix(argv[2]);
		sprintf(ep->outDir,"%s",argv[4]);
		sprintf(ep->outName,"%s",argv[5]);
		
		mainAlgoOfEPThresh();
	}

	if (argc == 8 && atoi(argv[3]) == MIN_MAX_KEYFRAMES_METRIC)//READ CURVE - Using MIN - MAX Algorithm
	{
//: %s THRESHOLD txtFileCurve 8 CurveDimension ThresholdPercentage  outDir outname
		metricChoise = MIN_MAX_KEYFRAMES_METRIC;
		printf("EP - MINMAX!!!Using data from Distance Curve and PreDefined Threshold!!!\n");
		M = getNumberMFromCurve(argv[2],atoi(argv[4]));
		refreshInitMemory(100,M,atoi(argv[4]),1,0);//Shmantikh allagh trexoume gia to +1
		ep->Threshold = atof(argv[1]);
		getDistanceMatrixFromCurveUsingMinMax(argv[2],atoi(argv[3]));

		sprintf(ep->outDir,"%s",argv[6]);
		sprintf(ep->outName,"%s",argv[7]);
		mainAlgoOfEPThresh();
	}
	else if (argc == 8)//READ CURVE - Thresholds
	{
//: %s THRESHOLD txtFileCurve metricChoise CurveDimension ThresholdPercentage  outDir outname
		metricChoise = atoi(argv[4]);

		printf("EP!!!Using data from Distance Curve and PreDefined Threshold!!!\n");
		M = getNumberMFromCurve(argv[2],atoi(argv[4]));
		refreshInitMemory(100,M,atoi(argv[4]),1,0);
		ep->Threshold = atof(argv[1]);
		getDistanceMatrixFromCurve(argv[2],atoi(argv[3]));

		sprintf(ep->outDir,"%s",argv[6]);
		sprintf(ep->outName,"%s",argv[7]);
		mainAlgoOfEPThresh();
	}



	return 0;
}

