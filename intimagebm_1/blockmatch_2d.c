// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright 2013, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
// All rights reserved.

#include <math.h>
#include "smartparameter.h"
#include "fimage2.h"
#include "utils.inc.c"

#include "blockmatch.h"

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

/************ BM UTILS ***************************************/

/** Shifts IN by (dx,dy) pixels stores it in OUT
 * the image is extended by symmetry at the boundaries*/

static Fimage *intshift(Fimage * in, Fimage * out, int dx, int dy)
{
   for (int c = 0; c < out->nch; c++)
      for (int y = 0; y < out->nrow; y++)
	 for (int x = 0; x < out->ncol; x++)
	    _vpla(out, x, y, c) = _vsympla(in, x + dx, y + dy, c);
   return out;
}



void bmIntSlice(Fimage * in1, Fimage * in2,
		Fimage * disp, Fimage * cost, Fimage * dispr,
		Fimage * costr, int patch_sz[], int hor_range[],
		int ver_range[], char *method)
{
   int nc1 = in1->ncol, nr1 = in1->nrow, nch = in1->nch;
   int nc2 = in2->ncol, nr2 = in2->nrow;
   int h = patch_sz[1];
   int halfh = h / 2;

   // NO NEED FOR 2d block matching : pass directly to stereoIntSlice
   if (ver_range[0] == 0 && ver_range[1] == 0) {
      stereoIntSlice(in1, in2, disp, cost, dispr, costr,
		     patch_sz, hor_range, ver_range, method);
      return;
   }

   /* Initialize the images to accumulate the results */
   for (int i = 0; i < cost->N; i++)
      _vv(cost, i) = INFINITY;
   for (int i = 0; i < costr->N; i++)
      _vv(costr, i) = INFINITY;

   Fimage *tmpin;

   Fimage *tmpdisp = malloc_fimage3pla(nc1, nr1, 2);
   Fimage *tmpcost = malloc_fimage3pla(nc1, nr1, 1);
   Fimage *tmpdispR = malloc_fimage3pla(nc2, nr2, 2);
   Fimage *tmpcostR = malloc_fimage3pla(nc2, nr2, 1);

   tmpin = malloc_fimage3pla(nc2, nr2, nch);
   for (int off_y = ver_range[0]; off_y <= ver_range[1]; off_y++) {
      /*shift right image (in2) vertically by off_y */
      intshift(in2, tmpin, 0, off_y);

      /*compute the horizontal disparities */
      int zero_ver_range[] = { 0, 0 };
      stereoIntSliceSubpix(in1, tmpin, tmpdisp, tmpcost, tmpdispR,
			   tmpcostR, patch_sz, hor_range, zero_ver_range,
			   method);

      /*accumulate the results */
      for (int x = 0; x < nc1; x++)
	 for (int y = 0; y < nr1; y++) {
	    if (_vpla(cost, x, y, 0) > _vpla(tmpcost, x, y, 0)) {
	       _vpla(cost, x, y, 0) = _vpla(tmpcost, x, y, 0);
	       _vpla(disp, x, y, 0) = _vpla(tmpdisp, x, y, 0);
	       _vpla(disp, x, y, 1) = _vpla(tmpdisp, x, y, 1) + off_y;
	    }
	 }
      for (int x = 0; x < nc2; x++)
	 for (int y = 0; y < nr2; y++) {
	    // is the patch center within the bounds of in2?
	    if ((y + off_y) >= halfh && (y + off_y) <= nr2 - h + halfh) {
	       if (_vpla(costr, x, y + off_y, 0) >
		   _vpla(tmpcostR, x, y, 0)) {
		  _vpla(costr, x, y + off_y, 0) = _vpla(tmpcostR, x, y, 0);
		  _vpla(dispr, x, y + off_y, 0) = _vpla(tmpdispR, x, y, 0);
		  _vpla(dispr, x, y + off_y, 1) =
		      _vpla(tmpdispR, x, y, 1) - off_y;
	       }
	    }
	 }
   }
   free(tmpin);

   free(tmpdisp);
   free(tmpcost);
   free(tmpdispR);
   free(tmpcostR);
}



/************ SUBPIXEL UTILS ***************************************/

#include "shear.h"

SMART_PARAMETER(SUBPIXEL, 1)

/* shift in1 vertically by q pixels
 * (could a noninteger factor) save the translated image in out */
static void vshift(Fimage * in1, Fimage * out, float q)
{
   int nc = in1->ncol, nr = in1->nrow, nch = in1->nch;

   float *ti = xmalloc(nc * nr * sizeof(float));
   float *to = xmalloc(nc * nr * sizeof(float));

   for (int c = 0; c < nch; c++) {
      for (int y = 0; y < nr; y++)
	 for (int x = 0; x < nc; x++)
	    ti[y + x * nr] = _vpla(in1, x, y, c);	//transpose while copying

#ifdef _OPENMP
#pragma omp critical		// some issue with fftw3 and thread safety
#endif
      image_shear(ti, to, nr, nc, 0., -q);

      for (int y = 0; y < nr; y++)
	 for (int x = 0; x < nc; x++)
	    _vpla(out, x, y, c) = to[y + x * nr];	//transpose while copying
   }
   free(ti);
   free(to);
}


void bmIntSliceSubpix(Fimage * in1, Fimage * in2,
		      Fimage * disp, Fimage * cost, Fimage * dispr,
		      Fimage * costr, int patch_sz[], int hor_range[],
		      int ver_range[], char *method)
{
   int nc1 = in1->ncol, nr1 = in1->nrow, nch = in1->nch;
   int nc2 = in2->ncol, nr2 = in2->nrow;

   float sub = SUBPIXEL();

   // NO NEED FOR SUBPIXEL: pass directly to bmIntSlice
   if (sub == 1) {
      bmIntSlice(in1, in2, disp, cost, dispr, costr,
		 patch_sz, hor_range, ver_range, method);
      return;
   }

   /* Initialize the images to accumulate the results */
   for (int i = 0; i < cost->N; i++)
      _vv(cost, i) = INFINITY;
   for (int i = 0; i < costr->N; i++)
      _vv(costr, i) = INFINITY;

   Fimage *tmpin;

   Fimage *tmpdisp = malloc_fimage3pla(nc1, nr1, 2);
   Fimage *tmpcost = malloc_fimage3pla(nc1, nr1, 1);
   Fimage *tmpdispR = malloc_fimage3pla(nc2, nr2, 2);
   Fimage *tmpcostR = malloc_fimage3pla(nc2, nr2, 1);

   /*shift right image j/sub  pixels to the up */
   tmpin = malloc_fimage3pla(nc2, nr2, nch);
   for (int j = 0; j < sub; j++) {
      float sh = 1. / sub * j;
      vshift(in2, tmpin, sh);
      bmIntSlice(in1, tmpin, tmpdisp, tmpcost, tmpdispR, tmpcostR,
		 patch_sz, hor_range, ver_range, method);

      /*accumulate */
      for (int y = 0; y < nr1; y++)
	 for (int x = 0; x < nc1; x++) {
	    if (_vpla(cost, x, y, 0) > _vpla(tmpcost, x, y, 0)) {
	       _vpla(cost, x, y, 0) = _vpla(tmpcost, x, y, 0);
	       _vpla(disp, x, y, 0) = _vpla(tmpdisp, x, y, 0);
	       _vpla(disp, x, y, 1) = _vpla(tmpdisp, x, y, 1) + sh;
	    }
	    if (j == 0 && x < nc2) {
	       // spares 1 of the iterations in the other loop
	       if (_vpla(costr, x, y, 0) > _vpla(tmpcostR, x, y, 0)) {
		  _vpla(costr, x, y, 0) = _vpla(tmpcostR, x, y, 0);
		  _vpla(dispr, x, y, 0) = _vpla(tmpdispR, x, y, 0);
		  _vpla(dispr, x, y, 1) = _vpla(tmpdispR, x, y, 1);
	       }
	    }
	 }
   }
   free(tmpin);

   /*shift left image i/sub pixels to the up */
   if (dispr) {
      tmpin = malloc_fimage3pla(nc1, nr1, nch);
      // the first value (i=0) was already done in the previous loop
      for (int i = 1; i < sub; i++) {
	 float sh = 1. / sub * i;
	 vshift(in1, tmpin, sh);
	 bmIntSlice(tmpin, in2, tmpdisp, tmpcost, tmpdispR, tmpcostR,
		    patch_sz, hor_range, ver_range, method);

	 /*accumulate */
	 for (int y = 0; y < nr2; y++)
	    for (int x = 0; x < nc2; x++) {
	       if (_vpla(tmpcostR, x, y, 0) < _vpla(costr, x, y, 0)) {
		  _vpla(costr, x, y, 0) = _vpla(tmpcostR, x, y, 0);
		  _vpla(dispr, x, y, 0) = _vpla(tmpdispR, x, y, 0);
		  _vpla(dispr, x, y, 1) = _vpla(tmpdispR, x, y, 1) + sh;
	       }
	    }
      }
      free(tmpin);
   }

   free(tmpdisp);
   free(tmpcost);
   free(tmpdispR);
   free(tmpcostR);
}
