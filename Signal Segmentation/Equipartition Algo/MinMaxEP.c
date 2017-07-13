/**
 * @author Panagiotakis Konstantinos
 * @version     v1.0
 */


#include "EP.h"

extern double getGlobalThreshold_CorrectDistanceMatrix();


void getDistanceMatrixFromCurveUsingMinMax(char *name, int choice)
{
	int i,j,k,k1,k2,k1best;
	FILE *fp1;
	float f;
	double max,min,d1,d2,d3,fd1,fd2,fd3,fd,maxmax = 0,minmin = INF;
	int dim = ep->n;
	double **Dist = (double **) Alloc2D(lines,sizeof(double*),lines,sizeof(double));
		
	fp1 = fopen(name, "r");

	for (i = 0; i < lines; ++i)
	for (k = 0; k < dim; ++k)
	{
		fscanf(fp1,"%f ",&f);
		ep->C[i][k] = f;
	}
	fclose(fp1);
	
	printf("lines = %d dim = %d\n",lines,dim);


	for (i = 0; i < lines; ++i)
	for (j = i; j < lines; ++j)
	{
			fd1 = 0;
			for (k = 0; k < dim; ++k)
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
			Dist[i][j] = fd;
			Dist[j][i] = fd;
	}
	printf("MIN_MAX_KEYFRAMES_METRIC is used\n");

		for (i = 0; i < lines; ++i)
		for (j = i; j < lines; ++j)
		{
			min = INF;
			k1best = i;
			for (k1 = i; k1 <= j; ++k1)
			{
				max = 0;
				for (k2 = i; k2 <= j; ++k2)
				{
					if (Dist[k1][k2] > max)
					{
						max = Dist[k1][k2];
						if (max > min)
						{
							break;
						}
					}
				}

				if (min > max)
				{
					min = max;
					k1best = k1;
				}
            }
 
   		ep->MinMaxKF[i][j] = k1best;
		ep->MinMaxKF[j][i] = k1best;
		ep->D[i][j] = min;
		ep->D[j][i] = min;
     
     /*if (i == 80)
        {
              printf("D[%d] = %2.2f ",j,min);
        }*/
       
		if (min > maxmax)
			maxmax = min;

		if (min > 0 && min < minmin)
			minmin = min;

   }		
   

	d3 = (double)(minmin/10.0);
	printf("ADD_EXP_TO_METRIC d3 = %f Max = %f\n",d3,maxmax);

	for (i = 0; i < lines; ++i)
	for (j = i+1; j < lines; ++j)
	{
		d3 = d3*(ABS(i-j))/(double)(lines);
		//d2 = (double)(d3*d1*exp(1-d1));
		ep->D[i][j] += d3;
		ep->D[j][i] += d3;
	}
	
	//gThresh = getGlobalThreshold_CorrectDistanceMatrix();

	FreeDouble(Dist,lines);

	printf("READING Curve M = %d N = %d dim = %d\n", ep->M,ep->N,ep->n);
	printf("lines = %d cols = %d dist(0,1) = %f \n",lines,cols,ep->D[0][1]);
}

/*


*/

void printfMinMaxSolutions(int *index)
{
	int i,j,x,y;
	int N = ep->N;
	char command[1000];
	FILE *fp1 = fopen("solutionsMinMax.txt", "w");

	fprintf(fp1,"%f %f \n",(double)ep->N, (double)ep->sol.numSol);

	for (i = 0; i < ep->sol.numSol; ++i)
	{
		for (j = 0; j <= ep->level+1; ++j)
		{
            x = (int)ep->sol.p[index[i]][j];
            y = (int)ep->sol.p[index[i]][j+1];
//            printf("%d - %d %d\n",x , y, ep->MinMaxKF[x][y]);
			fprintf(fp1,"%d ",ep->MinMaxKF[x][y]);
		}
		fprintf(fp1,"\n");
	}
	fclose(fp1);

    fp1 = fopen("MatrixMinMax.txt", "w");
	for (i = 0; i < lines; ++i)
	{
		for (j = 0; j < cols; ++j)
		{
			fprintf(fp1,"%d ",ep->MinMaxKF[i][j]);
		}
		fprintf(fp1,"\n");
	}
	fclose(fp1);

	sprintf(command,"copy solutionsMinMax.txt %s%s_solutionsMinMax.txt\n",ep->outDir,ep->outName);
	printf(">>%s \n",command);
	system(command);
	
    sprintf(command,"copy MatrixMinMax.txt %s%s_MatrixMinMax.txt\n",ep->outDir,ep->outName);
	printf(">>%s \n",command);
	system(command);
}


