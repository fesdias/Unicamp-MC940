#include "ift.h"
#include <math.h>

/* 
   Project 01, course MO445, Prof. Alexandre Falcao. 

   This code crops a region of interest containing a fingerprint with
   as minimum as possible surrounding noise.

   The code contains three functions iftAsfCOBin(), iftCloseBasins(),
   iftErodeBin() that are implemented by the Image Foresting
   Transform. Your task is to substitute them by the corresponding
   functions implemented by you. In order to do that, you should fill
   the code of the functions below. You should use iftAddFrame for
   padding zeroes and iftRemFrame to return the original image size
   whenever needed. Object pixels in these functions are pixels with
   value different from zero (e.g., usually 1 or 255) and background
   pixels are those with value equal to zero.  Object pixels are
   internal border pixels when they have a four-neighbor outside the
   object and background pixels are external border pixels when they
   have a four-neighbor inside an object.
 */


/* it returns pixels at the border of the image */

iftSet *MyImageBorder(iftImage *bin)
{
   iftSet *border = NULL;
   iftAdjRel *v4 = iftCreateAdjRel(1);

   iftImage *imgFrame = iftCopyImage(bin);
   imgFrame = iftAddFrame(bin, 1, 256);

   for (int p = 0; p < imgFrame->n; p++) {
      iftVoxel u = iftGetVoxelCoord(imgFrame, p);

      for (int i = 0; i < v4->n; i++){
         iftVoxel v = iftGetAdjacentVoxel(v4, u, i);

         if (imgFrame->val[p] != 0 && imgFrame->val[iftImgVoxelVal(imgFrame, v)] == 256) {
            iftInsertSet(border, v);
            break;
         }
      }
   }

   return border;
}

/* it returns a set with internal border pixels */

iftSet *MyObjectBorder(iftImage *bin)
{
   iftSet *border = NULL;
   iftAdjRel *v4 = iftCreateAdjRel(1);

   iftImage *imgFrame = iftCopyImage(bin);
   imgFrame = iftAddFrame(bin, 1, 0);

   for (int p = 0; p < imgFrame->n; p++) {
      iftVoxel u = iftGetVoxelCoord(imgFrame, p);

      for (int i = 0; i < v4->n; i++){
         iftVoxel v = iftGetAdjacentVoxel(v4, u, i);

         if (imgFrame->val[p] != 0 && imgFrame->val[iftImgVoxelVal(imgFrame, v)] == 0) {
            iftInsertSet(border, v);
            break;
         }
      }
   }

   return border;
}

/* it returns a set with external border pixels */

iftSet *MyBackgroundBorder(iftImage *bin)
{
   iftSet *border = NULL;
   iftAdjRel *v4 = iftCreateAdjRel(1);

   iftImage *imgFrame = iftCopyImage(bin);
   imgFrame = iftAddFrame(bin, 1, 255);

   for (int p = 0; p < imgFrame->n; p++) {
      iftVoxel u = iftGetVoxelCoord(imgFrame, p);

      for (int i = 0; i < v4->n; i++){
         iftVoxel v = iftGetAdjacentVoxel(v4, u, i);

         if (imgFrame->val[p] == 0 && imgFrame->val[iftImgVoxelVal(imgFrame, v)] != 0) {
            iftInsertSet(border, v);
            break;
         }
      }
   }

   return border;
}

/* it returns cost map */

iftImage *CostMap(iftImage *C, int value)
{
   for (int p = 0; p < C->n; p++)
      C->val[p] = value;

   return C;
}

/* it returns root map */

iftVoxel *VoxelMap(iftVoxel **R, int x, int y)
{
   for (int i = 0; i < x; i++) {
      for (int j = 0; i < y; y++)
         R[x][y] = NULL;
   }

   return R;
}

// Function to calculate distance
int distance(int x1, int y1, int x2, int y2)
{
    // Calculating distance
    return (pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

/* it dilates objects */

iftImage *MyDilateBin(iftImage *bin, iftSet **S, float radius)
{
   // create adjacent relation viz-8
   iftAdjRel *v8 = iftCreateAdjRel(sqrt(2));

   // create priority queue 
   iftGQueue *Q = iftCreateGQueue(1, 0, 0);

   // create cost matrix initialize in infinity 
   iftImage *C = iftCopyImage(bin);
   CostMap(C, 10000);

   // create root map 
   iftVoxel **R;
   VoxelMap(&R, bin->xsize, bin->ysize);

   // create dilatation mask 
   iftImage *D = iftCopyImage(&bin);

   // if doesn't have objects seed, create through object border 
   if (S == NULL)
      S = MyObjectBorder(bin);

   // use the seed to initialize cost map, root map, and priority queue 
   while (S != NULL) {
      iftVoxel p = iftGetVoxelCoord(bin, S[0]);

      // remove pixel 
      iftRemoveSetElem(S, p);

      // set cost to zero 
      C->val[iftGetVoxelIndex(bin, v)] = 0;

      // set root 
      R[iftGetVoxelIndex(bin, p)] = p;

      // insert voxel in the queue 
      iftInsertGQueue(Q, C[iftGetVoxelIndex(bin, p)]);
   }

   while (iftEmptyGQueue(Q) != 1) {

      // find cheapest element of queue 
      iftVoxel p = Q->C.first[1];

      if (C->val[iftGetVoxelIndex(bin, p)] <= pow(radius, 2)) {
         D->val[iftGetVoxelIndex(bin, p)] = bin->val[R->val[iftGetVoxelIndex(bin, p)]];

         // visit every viz-8
         for (int i = 0; i < v8->n; i++) {
            iftVoxel q = iftGetAdjacentVoxel(v8, p, i);

            if (C->val[iftGetVoxelIndex(bin, q)] > C->val[iftGetVoxelIndex(bin, p)] && bin->val[iftGetVoxelIndex(bin, q)] == 0) {
               
               // Calculate euclidean distance
               int q_x = iftGetXCoord(bin, iftGetVoxelIndex(bin, q));
               int q_y = iftGetYCoord(bin, iftGetVoxelIndex(bin, q));

               int Rp_x = iftGetXCoord(bin, iftGetVoxelIndex(bin, R[iftGetVoxelIndex(bin, p)]));
               int Rp_y = iftGetYCoord(bin, iftGetVoxelIndex(bin, R[iftGetVoxelIndex(bin, p)]));

               int tmp = distance(q_x, q_y, Rp_x, Rp_y);

               // set new costs 
               if (tmp < C->val[iftGetVoxelIndex(bin, q)]) {

                  iftRemoveGQueueElem(Q, C[iftGetVoxelIndex(bin, q)]);

                  C->val[iftGetVoxelIndex(bin, q)] = tmp;
                  R[iftGetVoxelIndex(bin, q)] = R[iftGetVoxelIndex(bin, p)];

                  iftInsertGQueue(Q, C[iftGetVoxelIndex(bin, q)]);
               }
            }
         }
      }

      // Coloca p em S
      else 
         iftUnionSetElem(S, p);
   }

   return D;
}

/* it erodes objects */

iftImage *MyErodeBin(iftImage *bin, iftSet **S, float radius)
{
   // create adjacent relation viz-8 
   iftAdjRel *v8 = iftCreateAdjRel(sqrt(2));

   // create priority queue 
   iftGQueue *Q = iftCreateGQueue(1, 0, 0);

   // create cost matrix initialize in infinity 
   iftImage *C = iftCopyImage(bin);
   CostMap(C, 10000);

   // create root map 
   iftVoxel **R;
   VoxelMap(&R, bin->xsize, bin->ysize);

   // create dilatation mask 
   iftImage *E = iftCopyImage(&bin);

   // if doesn't have objects seed, create through object border
   if (S == NULL)
      S = MyBackgroundBorder(bin);

   // use the seed to initialize cost map, root map, and priority queue 
   while (S != NULL) {
      iftVoxel p = iftGetVoxelCoord(bin, S[0]);

      // remove pixel 
      iftRemoveSetElem(S, p);

      // set cost to zero 
      C->val[iftGetVoxelIndex(bin, v)] = 0;

      // set root 
      R[iftGetVoxelIndex(bin, p)] = p;

      // insert voxel in the queue 
      iftInsertGQueue(Q, C[iftGetVoxelIndex(bin, p)]);
   }

   while (iftEmptyGQueue(Q) != 1) {

      // find cheapest element of queue 
      iftVoxel p = Q->C.first[1];

      if (C->val[iftGetVoxelIndex(bin, p)] <= pow(radius, 2)) {
         E->val[iftGetVoxelIndex(bin, p)] = bin->val[R->val[iftGetVoxelIndex(bin, p)]];

         // visit every viz-8 
         for (int i = 0; i < v8->n; i++) {
            iftVoxel q = iftGetAdjacentVoxel(v8, p, i);

            if (C->val[iftGetVoxelIndex(bin, q)] > C->val[iftGetVoxelIndex(bin, p)] && bin->val[iftGetVoxelIndex(bin, q)] != 0) {
               
               // Calculate euclidean distance 
               int q_x = iftGetXCoord(bin, iftGetVoxelIndex(bin, q));
               int q_y = iftGetYCoord(bin, iftGetVoxelIndex(bin, q));

               int Rp_x = iftGetXCoord(bin, iftGetVoxelIndex(bin, R[iftGetVoxelIndex(bin, p)]));
               int Rp_y = iftGetYCoord(bin, iftGetVoxelIndex(bin, R[iftGetVoxelIndex(bin, p)]));

               int tmp = distance(q_x, q_y, Rp_x, Rp_y);

               // set new costs 
               if (tmp < C->val[iftGetVoxelIndex(bin, q)]) {

                  iftRemoveGQueueElem(Q, C[iftGetVoxelIndex(bin, q)]);

                  C->val[iftGetVoxelIndex(bin, q)] = tmp;
                  R[iftGetVoxelIndex(bin, q)] = R[iftGetVoxelIndex(bin, p)];

                  iftInsertGQueue(Q, C[iftGetVoxelIndex(bin, q)]);
               }
            }
         }
      }

      // Coloca p em S
      else 
         iftUnionSetElem(S, p);
   }

   return E;
}

/* it executes dilation followed by erosion */

iftImage *MyCloseBin(iftImage *bin, float radius)
{
   iftImage *dil = NULL;
   iftImage *ero = NULL; 
   iftSet *seed = NULL;

   dil = MyDilateBin(bin, &seed, radius);
   ero = MyErodeBin(dil, &seed, radius);

   return ero;
}

/* it executes erosion followed by dilation */

iftImage *MyOpenBin(iftImage *bin, float radius)
{
   iftImage *ero = NULL; 
   iftImage *dil = NULL;
   iftSet *seed = NULL;

   ero = MyErodeBin(bin, &seed, radius);
   dil = MyDilateBin(ero, &seed, radius);
   
   return dil;
}

/* it executes closing followed by opening */

iftImage *MyAsfCOBin(iftImage *bin, float radius)
{
   iftImage *close = NULL; 
   iftImage *dil = NULL;
   iftSet *seed = NULL;

   iftImage *aux = iftAddFrame(bin, iftRound(radius) + 5, 0);
   dil = MyDilateBin(aux, &seed, radius);
   close = MyErodeBin(dil, &seed, 2.0f*radius);

   iftDestroyImage(&dil);
   iftDestroyImage(&aux);

   aux = MyDilateBin(close, &seed, radius);
   dil = iftRemFrame(aux, iftRound(radius) + 5);

   iftDestroyImage(&close);
   iftDestroyImage(&seed);
   iftDestroyImage(&aux);

   return dil;
}

/* it closes holes in objects */

iftImage *MyCloseBasins(iftImage *bin)
{
   // create adjacent relation viz-4
   iftAdjRel *v4 = iftCreateAdjRel(1);

   // create priority queue 
   iftGQueue *Q = iftCreateGQueue(1, 0, 0);

   // create cost matrix initialize in infinity
   iftImage *C = iftCopyImage(bin);
   CostMap(C, 10000);

   for (int i = 0; i < bin->n; i++) {
      iftVoxel p = iftGetVoxelCoord(bin, i);

      for (int i = 0; i < v4->n; i++) {
         iftVoxel q = iftGetAdjacentVoxel(v4, p, i);

         if (iftValidVoxel(bin, q) != 1) {
            C->val[iftGetVoxelIndex(bin, p)] = bin->val[iftGetVoxelIndex(bin, p)];
            iftInsertGQueue(Q, C[iftGetVoxelIndex(bin, q)]);
            break;
         }
      }
   }

   while (iftEmptyGQueue(Q) != 1) {

      // find cheapest element of queue 
      iftVoxel p = Q->C.first[1];
      iftRemoveGQueueElem(Q, C[iftGetVoxelIndex(bin, p)]);

      for (int i = 0; i < v4->n; i++) {
         iftVoxel q = iftGetAdjacentVoxel(v4, p, i);

         if (iftValidVoxel(bin, q) && C->val[iftGetVoxelIndex(bin, q)] > C->val[iftGetVoxelIndex(bin, p)]) {

            int tmp;

            if (C->val[iftGetVoxelIndex(bin, p)] > bin->val[iftGetVoxelIndex(bin, p)])
               tmp = C->val[iftGetVoxelIndex(bin, p)];

            else
               tmp = bin->val[iftGetVoxelIndex(bin, p)];

            // set new costs 
            if (tmp < C->val[iftGetVoxelIndex(bin, q)]) {
               C->val[iftGetVoxelIndex(bin, q)] = tmp;
               iftInsertGQueue(Q, C[iftGetVoxelIndex(bin, q)]);
            }
         }
      }
   }

   return C;
}

int main(int argc, char *argv[])
{
  timer *tstart=NULL;
  char   filename[200];
  
  /*--------------------------------------------------------*/

  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;

  /*--------------------------------------------------------*/

  
  if (argc != 3) {
    printf("project01 <P1> <P2>\n");
    printf("P1: folder with original images\n");
    printf("P2: folder with cropped images\n");
    exit(0);
  }

  tstart = iftTic();

  iftFileSet *fs   = iftLoadFileSetFromDirBySuffix(argv[1],".png", 1);  
  int nimages      = fs->n;
  char *out_dir    = argv[2];
  iftMakeDir(out_dir);
  iftAdjRel *A     = iftCircular(3.5), *B = iftCircular(1.5);
  
  for (int i=0; i < nimages; i++) {
    char *basename = iftFilename(fs->files[i]->path,".png");
    iftImage *orig = iftReadImageByExt(fs->files[i]->path);
    /* normalize  image */
    iftImage *norm = iftNormalize(orig,0,255);
    /* binarize image */
    iftImage *aux1 = iftBelowAdaptiveThreshold(norm, NULL, A, 0.98, 2, 255);
    /* remove noise components from the background */
    iftImage *aux2 = iftSelectCompAboveArea(aux1,B,100);
    iftDestroyImage(&aux1);
    /* apply morphological filtering to make the fingerprint the
       largest component: this operation must add frame and remove it
       afterwards. */
    aux1           = iftAsfCOBin(aux2,15.0);//MyAsfCOBin(aux2,15.0);
    iftDestroyImage(&aux2);
    /* close holes inside the components to allow subsequent erosion
       from the external borders only */
    aux2           = iftCloseBasins(aux1,NULL,NULL);//MyCloseBasins(aux1);
    iftDestroyImage(&aux1);
    /* erode components and select the largest one to estimate its
       center as close as possible to the center of the fingerprint */    
    iftSet *S = NULL;
    aux1           = iftErodeBin(aux2,&S,30.0);// MyErodeBin(aux2,&S,30.0);
    iftDestroySet(&S);
    iftDestroyImage(&aux2);
    aux2           = iftSelectLargestComp(aux1,B);

    /* crop the normalized image by the minimum bounding box of the
       resulting mask (largest component) */ 

    iftDestroyImage(&aux1);
    iftVoxel pos;
    iftBoundingBox bb = iftMinBoundingBox(aux2, &pos);    
    aux1              = iftExtractROI(norm,bb); 
      
    sprintf(filename,"%s/%s.png",out_dir,basename);
    iftWriteImageByExt(aux1,filename);
    iftDestroyImage(&aux1);
    iftDestroyImage(&aux2);
    iftDestroyImage(&orig);
    iftDestroyImage(&norm);
    iftFree(basename);
  }

  iftDestroyFileSet(&fs);
  iftDestroyAdjRel(&A);
  iftDestroyAdjRel(&B);
  
  puts("\nDone...");
  puts(iftFormattedTime(iftCompTime(tstart, iftToc())));
  
  /* ---------------------------------------------------------- */

  info = mallinfo();
  MemDinFinal = info.uordblks;
  if (MemDinInicial!=MemDinFinal)
    printf("\n\nDinamic memory was not completely deallocated (%d, %d)\n",
	   MemDinInicial,MemDinFinal);   

  return 0;
}
