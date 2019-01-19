// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright 2013, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
// All rights reserved.

#include <stdio.h>
#include <stdlib.h>
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


SMART_PARAMETER(SLICED, 32)	// SLICED=0 means no parallel handling,
    // and SLICED>1 forces a certain slice size
/************ SLICE AND DICE UTILS ***************************************/
// crop a rectangular region of size(w,h) from in
// (x0,y0) is the upper left corner of the crop region
// crop must be pre-allocated and w,h must be compatible with its size
void cropFimage(Fimage * in, Fimage * crop, int x0, int y0, int w, int h)
{
   if (in->ncol - x0 < w || in->nrow - y0 < h ||
       crop->ncol != w || crop->nrow != h || crop->nch != in->nch) {
      printf("crop error\n");
      exit(0);
   }
   for (int c = 0; c < in->nch; c++)
      for (int j = 0; j < h; j++)
	 for (int i = 0; i < w; i++)
	    _vpla(crop, i, j, c) = _vpla(in, i + x0, j + y0, c);
}

// paste a rectangular region of source into target.
// x0,y0 provides the alignement of target over the source, it is the position
// over target of the 0,0 pixel of source
// relative_x/ystart is the starting point of the copy,
// while (w,h) is the size of the area to be copied
void pasteFimage(Fimage * source, Fimage * target,
		 int x0, int y0, int relative_xstart, int relative_ystart,
		 int w, int h)
{
   if (relative_xstart + w > source->ncol ||
       relative_ystart + h > source->nrow ||
       x0 + relative_xstart < 0 ||
       y0 + relative_ystart < 0 ||
       x0 + relative_xstart + w > target->ncol ||
       y0 + relative_ystart + h > target->nrow) {
      printf("paste error\n");
      exit(0);
   }
   // printf("%d %d %d %d %d\n", x0,y0,relative_xstart, relative_ystart, w,h);
   for (int c = 0; c < source->nch; c++)
      for (int j = 0; j < h; j++)
	 for (int i = 0; i < w; i++)
	    _vpla(target, x0 + relative_xstart + i,
		  y0 + relative_ystart + j, c) =
		_vpla(source, relative_xstart + i, relative_ystart + j, c);
}



void stereoIntSlice(Fimage * in1, Fimage * in2,
		    Fimage * disp, Fimage * cost, Fimage * dispr,
		    Fimage * costr, int patch_sz[], int hor_range[],
		    int ver_range[], char *method)
{
   int nc1 = in1->ncol, nr1 = in1->nrow, nch = in1->nch;
   int nc2 = in2->ncol;

   //const int MIN_SLICEH = 32; // old 64
   //I've observed that with slices 2~3x h are the best performing
   int MIN_SLICEH = SLICED();
   int w = patch_sz[0];
   int h = patch_sz[1];
   int sliceh = fmax(2 * h, MIN_SLICEH);
   int halfh_up = h / 2;
   int halfh_down = h - (int) (h / 2) - 1;

   printf("sliced:%d window=(%d,%d) range=[%d,%d] x [%d,%d]:",
	  sliceh, w, h, hor_range[0], hor_range[1], ver_range[0],
	  ver_range[1]);

   // NO NEED FOR SLICING: pass directly to stereoInt
   if (sliceh >= nr1) {
      stereoInt(in1, in2, disp, cost, dispr, costr,
		patch_sz, hor_range, ver_range, method);
      printf("\n");		// ...
      return;
   }

   if (ver_range[0] != 0 || ver_range[1] != 0) {
      fprintf(stderr,
	      "  ATTENTION: SLICED, only works for horizontal offs.\n");
      ver_range[0] = 0;
      ver_range[1] = 0;
   }
#ifdef _OPENMP
#pragma omp parallel for
#endif
   for (int s = 0; s < nr1; s += sliceh) {

      // compute slice and padding parameterr
      int ymin = fmax(0, s - halfh_up);
      int ymax = fmin(s + sliceh + halfh_down, nr1);
      int top_padding = s - ymin;
      int bottom_padding = ymax - fmin(s + sliceh, nr1);

      Fimage *slicein1 = malloc_fimage3pla(nc1, ymax - ymin, nch);
      Fimage *slicein2 = malloc_fimage3pla(nc2, ymax - ymin, nch);
      Fimage *slicedisp = malloc_fimage3pla(nc1, ymax - ymin, 2);
      Fimage *slicecost = malloc_fimage3pla(nc1, ymax - ymin, 1);
      Fimage *slicedispr = NULL, *slicecostr = NULL;
      if (dispr)
	 slicedispr = malloc_fimage3pla(nc2, ymax - ymin, 2);
      if (costr)
	 slicecostr = malloc_fimage3pla(nc2, ymax - ymin, 1);

      cropFimage(in1, slicein1, 0, ymin, nc1, ymax - ymin);
      cropFimage(in2, slicein2, 0, ymin, nc2, ymax - ymin);

      stereoInt(slicein1, slicein2, slicedisp, slicecost, slicedispr,
		slicecostr, patch_sz, hor_range, ver_range, method);

#ifdef _OPENMP
#pragma omp critical
#endif
      pasteFimage(slicedisp, disp,
		  0, ymin, 0, top_padding, nc1,
		  ymax - ymin - top_padding - bottom_padding);
#ifdef _OPENMP
#pragma omp critical
#endif
      pasteFimage(slicecost, cost,
		  0, ymin, 0, top_padding, nc1,
		  ymax - ymin - top_padding - bottom_padding);
#ifdef _OPENMP
#pragma omp critical
#endif
      if (slicedispr)
	 pasteFimage(slicedispr, dispr,
		     0, ymin, 0, top_padding, nc2,
		     ymax - ymin - top_padding - bottom_padding);
#ifdef _OPENMP
#pragma omp critical
#endif
      if (slicecostr)
	 pasteFimage(slicecostr, costr,
		     0, ymin, 0, top_padding, nc2,
		     ymax - ymin - top_padding - bottom_padding);

      free(slicein1);
      free(slicein2);
      free(slicedisp);
      free(slicecost);
      free(slicedispr);
      free(slicecostr);
   }
   printf("\n");		// end ...
}








/************ SUBPIXEL UTILS ***************************************/

#include "shear.h"


SMART_PARAMETER(SUBPIXEL, 1)



/* shift in1 horizontally by q pixels
 * (could a noninteger factor) save the translated image in out */
static void shift(Fimage * in1, Fimage * out, float q)
{
   int nc = in1->ncol, nr = in1->nrow, nch = in1->nch;

   float *ti = xmalloc(nc * nr * sizeof(float));
   float *to = xmalloc(nc * nr * sizeof(float));

   for (int c = 0; c < nch; c++) {
      for (int y = 0; y < nr; y++)
	 for (int x = 0; x < nc; x++)
	    ti[x + y * nc] = _vpla(in1, x, y, c);

#ifdef _OPENMP
#pragma omp critical		// some issue with fftw3 and thread safety
#endif
      image_shear(ti, to, nc, nr, 0., -q);

      for (int y = 0; y < nr; y++)
	 for (int x = 0; x < nc; x++)
	    _vpla(out, x, y, c) = to[x + y * nc];
   }
   free(ti);
   free(to);
}



void stereoIntSliceSubpix(Fimage * in1, Fimage * in2,
			  Fimage * disp, Fimage * cost, Fimage * dispr,
			  Fimage * costr, int patch_sz[], int hor_range[],
			  int ver_range[], char *method)
{
   int nc1 = in1->ncol, nr1 = in1->nrow, nch = in1->nch;
   int nc2 = in2->ncol, nr2 = in2->nrow;

   float sub = SUBPIXEL();

   // NO NEED FOR SUBPIXEL: pass directly to stereoIntSlice
   if (sub == 1) {
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

   /*shift right image j/sub  pixels to the left */
   tmpin = malloc_fimage3pla(nc2, nr2, nch);
   for (int j = 0; j < sub; j++) {
      float sh = 1. / sub * j;
      shift(in2, tmpin, sh);
      stereoIntSlice(in1, tmpin, tmpdisp, tmpcost, tmpdispR, tmpcostR,
		     patch_sz, hor_range, ver_range, method);

      /*accumulate */
      for (int y = 0; y < nr1; y++)
	 for (int x = 0; x < nc1; x++) {
	    if (_vpla(cost, x, y, 0) > _vpla(tmpcost, x, y, 0)) {
	       _vpla(cost, x, y, 0) = _vpla(tmpcost, x, y, 0);
	       _vpla(disp, x, y, 0) = _vpla(tmpdisp, x, y, 0) + sh;
	       _vpla(disp, x, y, 1) = _vpla(tmpdisp, x, y, 1);
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

   /*shift left image i/sub pixels to the left */
   if (dispr) {
      tmpin = malloc_fimage3pla(nc1, nr1, nch);
      // the first value (i=0) was already done in the previous loop
      for (int i = 1; i < sub; i++) {
	 float sh = 1. / sub * i;
	 shift(in1, tmpin, sh);
	 stereoIntSlice(tmpin, in2, tmpdisp, tmpcost, tmpdispR, tmpcostR,
			patch_sz, hor_range, ver_range, method);


	 //keep clang --analyze happy
	 if (NULL == tmpcostR || NULL == tmpdispR)
	    error("NULL pointers in stereoIntSliceSubpix()");
	 /*accumulate */
	 for (int y = 0; y < nr2; y++)
	    for (int x = 0; x < nc2; x++) {
	       if (_vpla(tmpcostR, x, y, 0) < _vpla(costr, x, y, 0)) {
		  _vpla(costr, x, y, 0) = _vpla(tmpcostR, x, y, 0);
		  _vpla(dispr, x, y, 0) = _vpla(tmpdispR, x, y, 0) + sh;
		  _vpla(dispr, x, y, 1) = _vpla(tmpdispR, x, y, 1);
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
