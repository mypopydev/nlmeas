// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright 2013, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
// All rights reserved.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "utils.inc.c"

#include "fimage2.h"


/*
typedef struct fimage {
   int sz[3];    //[col,row,ch]
   int ncol;
   int nrow;
   int nch;
   int N;
   int planar;
   float *gray;  // must mantain both 
   float *v;     // same as gray
//   float v[];    // flexible array is another possibility
} Fimage;

// determines how the data should be interpreted
#define FIMAGE_PLANAR 1
#define FIMAGE_VEC 0
*/

// multi channel 2d image
Fimage *malloc_fimage4(int nc, int nr, int nch, int planar)
{
   Fimage *f = xcalloc(1, sizeof(*f) + nc * nr * nch * sizeof(float));
   f->v = (float *) (f + 1);
   f->ncol = nc;
   f->nrow = nr;
   f->nch = nch;
   f->planar = planar;
   f->gray = f->v;
   f->sz[0] = f->ncol;
   f->sz[1] = f->nrow;
   f->sz[2] = f->nch;
   f->N = nc * nr * nch;
   return f;
}

// 1d vector
Fimage *malloc_fimage1(int numel)
{
   return malloc_fimage4(numel, 1, 1, FIMAGE_PLANAR);
}

// single channel 2d image
Fimage *malloc_fimage2(int nc, int nr)
{
   return malloc_fimage4(nc, nr, 1, FIMAGE_PLANAR);
}

// multi channel 2d image planar   rrrr ggggg bbbbb
Fimage *malloc_fimage3pla(int nc, int nr, int nch)
{
   return malloc_fimage4(nc, nr, nch, FIMAGE_PLANAR);
}

// multi channel 2d image vectorial   rgb rgb rgb rgb
Fimage *malloc_fimage3vec(int nc, int nr, int nch)
{
   return malloc_fimage4(nc, nr, nch, FIMAGE_VEC);
}

// this is very strange, the structure is allocate only 
Fimage *malloc_fimageWrap(float *v, int nc, int nr, int nch, int planar)
{
   Fimage *f = xcalloc(1, sizeof(*f));
   f->v = v;
   f->ncol = nc;
   f->nrow = nr;
   f->nch = nch;
   f->planar = planar;
   f->gray = f->v;
   f->sz[0] = f->ncol;
   f->sz[1] = f->nrow;
   f->sz[2] = f->nch;
   f->N = nc * nr * nch;
   return f;
}




// in place transform from vec to planar form of the array
void fimage_vec_to_planar(Fimage * u)
{
   int nc = u->ncol, nr = u->nrow, nch = u->nch;
   int x, y, c;
   // copy to a temporary place
   float *tmp = xcalloc(1, nc * nr * nch * sizeof *tmp);
   for (x = 0; x < nc * nr * nch; x++)
      tmp[x] = u->v[x];
   // rearrange
   for (x = 0; x < nc; x++)
      for (y = 0; y < nr; y++)
	 for (c = 0; c < nch; c++)
	    u->v[(x + y * nc) + nc * nr * c] = tmp[(x + y * nc) * nch + c];

   u->planar = 1;
   free(tmp);
}

// in place transform from planar to vec form of the array
void fimage_planar_to_vec(Fimage * u)
{
   int nc = u->ncol, nr = u->nrow, nch = u->nch;
   int x, y, c;
   // copy to a temporary place
   float *tmp = xcalloc(1, nc * nr * nch * sizeof *tmp);
   for (x = 0; x < nc * nr * nch; x++)
      tmp[x] = u->v[x];
   // rearrange
   for (x = 0; x < nc; x++)
      for (y = 0; y < nr; y++)
	 for (c = 0; c < nch; c++)
	    u->v[(x + y * nc) * nch + c] = tmp[(x + y * nc) + nc * nr * c];

   u->planar = 0;
   free(tmp);
}

void fimage_to_vec(Fimage * u)
{
   if (u->planar == 1)
      fimage_planar_to_vec(u);
}

void fimage_to_planar(Fimage * u)
{
   if (u->planar == 0)
      fimage_vec_to_planar(u);
}


/* MACROS FOR ACCESSING FIMAGE PIXELS                          */
/* In this file all the multi-channel Fimages are "planar",    */
/* that means channels are stored as: rrrr ggggg bbbbb         */

/* _p* compute the array position of pixel in Fimage given its coordinates    */
/* _v* access the pixel value of Fimage at a given coordinate                 */
/* _ppla  index of pixel at position (i,j,c) ASSUMING that u IS PLANAR        */
#define _ppla(u,i,j,c) (  (i) + (j)*(u)->ncol + (u)->ncol*(u)->nrow*(c) )
/* _vpla  access value at position (i,j,c) ASSUMING that u IS PLANAR          */
#define _vpla(u,i,j,c) (  (u)->v[ _ppla(u,i,j,c) ] )
/* _vv    access value at position [i] of the array                           */
#define _vv(u,i)    (  (u)->v[ i ] )


///* MACROS FOR SYMMETRIZED ACCESS             */
///* dct type II symmetry at the boundary      */
///* _psym*  position of the symmetrized pixel */
///* _vsym*  value of the symmetrized pixel    */
#define _psympla(u,x,y,c)    (   _ppla(  u,                                    \
         ( (x)<0 ? -(x)-1: ( (x) >= (u)->ncol ? -(x)+2*(u)->ncol-1 : (x)  ) ), \
         ( (y)<0 ? -(y)-1: ( (y) >= (u)->nrow ? -(y)+2*(u)->nrow-1 : (y)  ) ), \
         (c) < (u)->nch ? c: -1  )   )
#define _vsympla(u,x,y,c)    (  (u)->v[ _psympla(u,x,y,c) ]  )

///* MACROS FOR TESTING INTERIOR             */
///* return the position or -1 if not inside */
#define _pinpla(u,i,j,c) ( (i)<(u)->ncol && (j)<(u)->nrow &&                   \
      (c)<(u)->nch && (i)>=0 && (j)>=0 && (c)>=0 ? _ppla(u,i,j,c) : -1 )
