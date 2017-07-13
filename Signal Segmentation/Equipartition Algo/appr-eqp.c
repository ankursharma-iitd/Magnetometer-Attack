/**
 * @author Panagiotakis Konstantinos
 * @version     v1.0
 */



#include "EP.h"


void getPointsFromP(int *points,int **p)
{
     int M = ep->M;     
     int N = ep->N+1; 
     int i,f,k;
     
     points[0] = M-1;
     
     f = M-1;
     k = 0;
     for (i = N-1; i >= 1; --i)
     {
         k++;
         points[k] = p[f][i];
         f = p[f][i];
     }
}

double getC(double **s)
{
  int M = ep->M;     
  double bestDist = 0,temp = 0;
  int i,j,i1,i2;
  
  
  for (i = 1; i < M-1; ++i)
  for (j = 1; j < M-1; ++j)
  {
      for (i1=-1; i1 <= 1; ++i1)    
      for (i2=-1; i2 <= 1; ++i2)
      {
          if (i1 != 0 || i2 != 0)
          {
             temp = ABS(s[i][j] - s[i+i1][j+i2]) / (sqrt(i1*i1+i2*i2));
             if (temp > bestDist)
               bestDist = temp;  
          }
      }
  }
  
  printf("bestDist = %f %f\n",bestDist,ABS(s[3][3]-s[3][4]));
  return bestDist;
}

void runApprEQP()
{
	int i,j,k,u,v,l;
	int N = ep->N+1;
	int M = ep->M,sum;
	char command[1000];
	FILE *fp1 = fopen("solutionsEQP.txt", "w");
    int **p = (int **)(int **) Alloc2D(M,sizeof(int*),N+1,sizeof(int));
    double **s = ep->D,bestDist = 0,lowDist = INF,c,e,dx,dy;
    int *points = (int *)calloc(N+1,sizeof(int));
    int *bestPoints = (int *)calloc(N+1,sizeof(int));
    int *bestMinPoints = (int *)calloc(N+1,sizeof(int));    
 
 /*   for (i = 0; i < M; ++i)
    for (j = 0; j < M; ++j)
    {
        dx = exp((float)i/M)*cos(2*Pi*1.5*((float)i/M) ) - exp((float)j/M)*cos(2*Pi*1.5*((float)j/M) );
        dy = exp((float)i/M)*sin(2*Pi*1.5*((float)i/M) ) - exp((float)j/M)*sin(2*Pi*1.5*((float)j/M) );
        s[i][j] = sqrt(dx*dx+dy*dy);
    }
    */
	fprintf(fp1,"%f \n",(double)ep->N);
    printf("s[1][2] = %f \n",s[1][2]);

   c =  getC(s);
   e = 10.0*c/M;
printf("M = %d N = %d c = %f e = %f\n",M,N,c,e);
    k = -1;
    while (1)
    {
        k++;
        
        for (u = 0; u < M; ++u)
        for (l = 0; l < N+1; ++l)
        {
            p[u][l] = -1;
        }
        
        p[0][0] = 1;
   
        sum = 0;
        for (u = 0; u < M; ++u)
        for (v = u; v < M; ++v)
        {
            if (ABS(s[u][v] -s[0][k]) < e/2)
            {
            
              for (l = 0; l < N; ++l)
                if (p[u][l] != -1)
                {
                    p[v][l+1] = u;
                    sum++;
                }
            }
        }                
        
        
        if (p[M-1][N-1] != -1)
        {
                  
           getPointsFromP(points,p);
           if (s[points[N-2]][points[N-1]] >  bestDist)
           {
              bestDist = s[points[N-2]][points[N-1]];
              for (i = 0; i < N; ++i)
                  bestPoints[i] = points[i];                           
           }
           if (s[points[N-2]][points[N-1]] < lowDist)
           {
              lowDist = s[points[N-2]][points[N-1]];
              for (i = 0; i < N; ++i)
                  bestMinPoints[i] = points[i];                           
           } 
           printf("!!!%d: %f %f\n",k,((float)sum)/(M*N),bestDist);          
        }
        
        //printf("%d: %f %f\n",k,((float)sum)/(M*N),bestDist);
   
       if (k == M-1)
       {
             
          if (bestDist == 0)
          {
           k = -1;
           e = e*2;
          }
          else
           break;
       }
    }

	for (i = N-1; i >= 0; --i)
	{
	     fprintf(fp1,"%d ",bestPoints[i]);
	}
	fprintf(fp1,"\n");
	fprintf(fp1,"\n");
	for (i = N-1; i >= 0; --i)
	{
	     fprintf(fp1,"%d ",bestMinPoints[i]);
	}
	fprintf(fp1,"\n");
	fclose(fp1);

	sprintf(command,"copy solutionsEQP.txt %s%s_solutionsEQP.txt\n",ep->outDir,ep->outName);
	printf(">>%s \n",command);
	system(command);
	
   free(points);
   free(bestPoints);
   free(bestMinPoints);
   for (i = 0; i < M; ++i)
       free(p[i]);
   
   free(p);
}


