#include <stdio.h>
#include <stdlib.h>
#include "alloc.h"


/*______________________________________________________________
 |    Generic function for allocating initialized memory.       |
 |    Failure causes exiting the program.                       |
 |______________________________________________________________|
*/
void *Calloc(int nobjs , int nbytes)
{
    void *sp;

    sp = (void * )calloc(nobjs , nbytes);
    if (sp == NULL)
    {
        fprintf(stderr, "calloc : Error in allocating memory %d.",nbytes);
        exit(1);
    }
    return sp;
}

/*_________________________________________________________
 | name:    Alloc2D
 | args:    1) dimx: # objects in x dimension
 |          2) dimx_bytes: size of each object in
 |                         x dimension(in bytes)
 |          3) dimy: # objects in y dimension
 |          4) dimy_bytes: size of each object in
 |                         y dimension(in bytes)
 | return:  A generic pointer to a 2D space
 | usage:   Allocates a zero initialized 2D space for any type.
 |_________________________________________________________
*/
void **Alloc2D(int dimx, int dimx_bytes, int dimy, int dimy_bytes)
{
    int i;
    void **Image;

    Image = (void **)Calloc(dimx, dimx_bytes);
    for (i=0; i<dimx; i++)
    {
        Image[i] = (void *)Calloc(dimy, dimy_bytes);
    }
    return Image;
}

/*_________________________________________________________
 | name:    Alloc3D
 | args:    1) dimx: # objects in x dimension
 |          2) dimx_bytes: size of each object in
 |                         x dimension(in bytes)
 |          3) dimy: # objects in y dimension
 |          4) dimy_bytes: size of each object in
 |                         y dimension(in bytes)
 | return:  A generic pointer to a 3D space
 | usage:   Allocates a zero initialized 3D space for any type.
 |_________________________________________________________
*/
void ***Alloc3D(int dimx, int dimx_bytes, int dimy, int dimy_bytes, int dimz, 
               int dimz_bytes)
{
    int i,j;
    void ***Image;

    Image = (void ***)Calloc(dimx, dimx_bytes);
    for (i=0; i<dimx; i++)
    {
        Image[i] = (void **)Calloc(dimy, dimy_bytes);
        for (j=0; j<dimy; j++)
            Image[i][j] = (void *)Calloc(dimz, dimz_bytes);
    }
    return Image;
}



/*_________________________________________________
 |    Prints a 2D float Image into a file          |
 |_________________________________________________|
*/
void Print2D(float **I, int dimx, int dimy, char *file)
{
   int i;
   int j;
   FILE *fp;

   fp = fopen(file, "w");
   for (i=0; i<dimx; i++)
   {
      for (j=0; j<dimy; j++)
      {
	 fprintf(fp, "%f ", I[i][j]);
	 if (j %10 == 0) fprintf(fp, "\n");
      }
   }
   fclose(fp);
}

/*________________________________________________________
 |    Destroys a 2D matrix consisting of dimx rows.       |
 |________________________________________________________|
*/
int **FreeInt(int **I, int dimx)
{
   int i;

   for (i=0; i<dimx; i++)
   {
      free((void *)I[i]);
   }
   free((void *)I);
   return NULL;
}

/*________________________________________________________
 |    Destroys a 2D matrix consisting of dimx rows.       |
 |________________________________________________________|
*/
float **Free(float **I, int dimx)
{
   int i;

   for (i=0; i<dimx; i++)
   {
      free((void *)I[i]);
   }
   free((void *)I);
   return NULL;
}


/*________________________________________________________
 |    Destroys a 2D matrix consisting of dimx rows.       |
 |________________________________________________________|
*/
double **FreeDouble(double **I, int dimx)
{
   int i;

   for (i=0; i<dimx; i++)
   {
      free((void *)I[i]);
   }
   free((void *)I);
   return NULL;
}


/*_________________________________________________________
 | Destroys a three dimension matrix. 
 | dimx, dimy = the first 2 dimensions of the matrix.
 | I = the matrix to be destroyed.
 |_________________________________________________________
*/
int ***FreeInt3D(int ***I, int dimx, int dimy)
{
   int i,j;

   for (i=0; i<dimx; i++)
   {
      for (j=0; j<dimy; j++)
         free((void *)I[i][j]);

      free((void *)I[i]);
   }
   free((void *)I);
   return NULL;
}



/*_________________________________________________________
 | Destroys a three dimension matrix. 
 | dimx, dimy = the first 2 dimensions of the matrix.
 | I = the matrix to be destroyed.
 |_________________________________________________________
*/
float ***Free3D(float ***I, int dimx, int dimy)
{
   int i,j;

   for (i=0; i<dimx; i++)
   {
      for (j=0; j<dimy; j++)
         free((void *)I[i][j]);

      free((void *)I[i]);
   }
   free((void *)I);
   return NULL;
}

/*____________________________________________________________
 | Allocate memory for a tree. The levels of the tree = depth |
 | and each level has 4^level image elements with dimensions  |
 | l/(2^level) x c/(2^level)                                  |
 |____________________________________________________________|
*/
Level *AllocTree(int depth, int l, int c)
{
   int i,j;
   List *item, *tmp1;
   Level *tmp, *m, *m1;

   m=(Level *) malloc(sizeof(Level));
   m->head=(List *) malloc(sizeof(List));
   m->tail=(List *) malloc(sizeof(List));
   m->level=0;
   m->head->l=l;
   m->head->c=c;
   m->tail->l=l;
   m->tail->c=c;
   m->head->image=(float **)Alloc2D(l,sizeof(float *),c,sizeof(float));
   m->tail->image=(float **)Alloc2D(l,sizeof(float *),c,sizeof(float));
   m->next=NULL;
   m->head->next=NULL;
   m->tail->next=NULL;

   m1=m;
   for (i=1; i<=depth; i++)
   {
      tmp=(Level *) malloc(sizeof(Level));
      tmp->next=NULL;
      tmp->head=NULL;
      tmp->tail=NULL;

      for (item=m1->head; item!=NULL; item=item->next)
      {
         for (j=0; j<4; j++)
         {
            tmp1=(List *)malloc(sizeof(List));
            tmp1->l=(item->l)/2;
            tmp1->c=(item->c)/2;
            tmp1->image=(float **)Alloc2D(item->l/2,sizeof(float *),
                                          item->c/2,sizeof(float));
            tmp1->next=NULL;
            if ((tmp->head)==NULL)
                  tmp->head=tmp1;
            else
                  tmp->tail->next=tmp1;
            tmp->tail=tmp1;
         }
      }
      tmp->level=i;
      m1->next=tmp;
      m1=tmp;
   }

   return m;
}

/*_______________________________________
 | Destroys a tree pointed by root       |
 |_______________________________________|
*/
Level *FreeTree(Level *root)
{
   Level *r, *s;
   List  *p, *q;

   for (r=root; r!=NULL; r=s)
   {
      s=r->next;
      for (p=r->head; p!=NULL; p=q)
      {
         q=p->next;
         p->image=Free(p->image,p->l);
         free(p);
      }
      free(r);
   }
   return NULL;
}
