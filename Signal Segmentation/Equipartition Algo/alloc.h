/*______________________________________________________
 |The two structures below implement a tree as a list of|
 |lists. The list of type Level represents the depth of |
 |the tree while the list of type List contains a list  |
 |with the members of each level of the tree            |
 |______________________________________________________|
*/
typedef struct Level {
        int          level;
        struct List  *head;
        struct List  *tail;
        struct Level *next;
} Level;

typedef struct List {
        float       **image;
        int         l,c;
        struct List *next;
} List;


void *Calloc(int nobjs , int nbytes);
void **Alloc2D(int dimx, int dimx_bytes, int dimy, int dimy_bytes);
void ***Alloc3D(int dimx, int dimx_bytes, int dimy, int dimy_bytes, int dimz, 
               int dimz_bytes);
void Print2D(float **I, int dimx, int dimy, char *file);
int **FreeInt(int **I, int dimx);
double **FreeDouble(double **I, int dimx);
int ***FreeInt3D(int ***I, int dimx,int dimy);
float **Free(float **I, int dimx);
float ***Free3D(float ***I, int dimx, int dimy);

Level *AllocTree(int depth, int l, int c);
Level *FreeTree(Level *root);
