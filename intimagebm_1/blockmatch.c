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
#include "utils.inc.c"

#include "integralImage4.c"
#include "fimage2.h"

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

// This structure contains the data needed by the block matching algorithm
// Concretely: * pointers to functions that determine the patch costs
//             * references to the input images, patch size and current offset
//             * references to precomputed images (if needed by the method)
struct t_struct_BlockMatchingStuff;

//
// To functions determine which patch cost is used for the block matching.
// * The first one computes the content of the integral image: absolute
//   differences, squared differnces,  pointwise products.
//   It translates in2 by [off_x,off_y] and computes the difference with in1
//   The resulting "difference image" (s->diff), is stored in the
//   structure t_struct_BlockMatchingStuff, while the integral image
//   itself is computed from s->diff in the main function.
//   This function is also responsible of precomputing the images
//   needed for the evaluation of the cost and store them in the structure.
//   Also updates in s the values of off_x/y.
// * The second function evaluates the integral image at a certain window
//   determined by (x,y,w,h). The computation usually requires to evaluate
//   some precomputed image, so the structure must be passed to the function
//   as well as the integral image (iifloat_t).
typedef void (*t_func_integralImageDiff) (Fimage *, Fimage *,
					  int, int,
					  struct
					  t_struct_BlockMatchingStuff *);
typedef long double (*t_func_integralImageEval) (struct iimage_t *, int,
						 int, int, int, int,
						 struct
						 t_struct_BlockMatchingStuff
						 *);

struct t_struct_BlockMatchingStuff {
   // Pointers to function that implement a specific cost
   t_func_integralImageDiff iiDiff;	//
   t_func_integralImageEval iiEval;	//

   Fimage *in1;			// reference to the original images needed
   Fimage *in2;			//  to gemerate the precomputed
   int patch_sz[2];		// patch (of block) size
   int offx, offy;		// current offset, changed t_func_integralImageDiff
   iifloat_t *diff;		// durrent diff image,
   int numdiff;			// number of "channels" of the diff image

   // precomputed images needed by some of the methods
   Fimage *mu1;
   Fimage *mu2;
   Fimage *std1;
   Fimage *std2;
   Fimage *norm1;
   Fimage *norm2;

   Fimage *BTmax1;
   Fimage *BTmin1;
   Fimage *BTmax2;
   Fimage *BTmin2;
};


/* Contain t_func_integralImageDiff and t_func_integralImageEval functions:
 *  integralImageDiff_{SAD,SSD,BTSSD,BTSAD,ZSSDc,ZSSD,SSDNorm,NCCc,NCC,LIN,AFF}
 *  integralImageEval_{SSD,ZSSD,SSDNorm,NCC,ZSSDc,NCCc,LIN,AFF}*/
#include "blockmatch_dists.h"



/** prepares the structure t_struct_BlockMatchingStuff:
 *   inits pointers to null,
 *   sets iiDiff & iiEval: function pointrs determined from method
 *   sets the patch size and saves the reference to the input images in1, in2*/
static void setup_BlockMatchingStuff(struct t_struct_BlockMatchingStuff *s,
				     Fimage * in1, Fimage * in2,
				     int patch_sz[], char *method)
{
   // keep a reference to the input images and patch sizes
   s->in1 = in1;
   s->in2 = in2;
   s->patch_sz[0] = patch_sz[0];
   s->patch_sz[1] = patch_sz[1];

   // important: initialize the variables updated by iiDiff to NULL
   s->mu1 = s->mu2 = NULL;
   s->std1 = s->std2 = NULL;
   s->norm1 = s->norm2 = NULL;
   s->BTmax1 = s->BTmax2 = NULL;
   s->BTmin1 = s->BTmin2 = NULL;
   s->diff = NULL;

   //determine what to do and assign pointer functions
   if (strcmp(method, "SAD") == 0) {
      s->iiDiff = &integralImageDiff_SAD;
      s->iiEval = &integralImageEval_SSD;
   } else if (strcmp(method, "BTSAD") == 0) {
      s->iiDiff = &integralImageDiff_BTSAD;
      s->iiEval = &integralImageEval_SSD;
   } else if (strcmp(method, "BTSSD") == 0) {
      s->iiDiff = &integralImageDiff_BTSSD;
      s->iiEval = &integralImageEval_SSD;
   } else if (strcmp(method, "BT2DSAD") == 0) {
      s->iiDiff = &integralImageDiff_BT2DSAD;
      s->iiEval = &integralImageEval_SSD;
   } else if (strcmp(method, "BT2DSSD") == 0) {
      s->iiDiff = &integralImageDiff_BT2DSSD;
      s->iiEval = &integralImageEval_SSD;
   } else if (strcmp(method, "SSD") == 0) {
      s->iiDiff = &integralImageDiff_SSD;
      s->iiEval = &integralImageEval_SSD;
   } else if (strcmp(method, "ZSSD") == 0) {
      s->iiDiff = &integralImageDiff_ZSSD;
      s->iiEval = &integralImageEval_ZSSD;
   } else if (strcmp(method, "SSDNorm") == 0) {
      s->iiDiff = &integralImageDiff_SSDNorm;
      s->iiEval = &integralImageEval_SSDNorm;
   } else if (strcmp(method, "LIN") == 0) {
      s->iiDiff = &integralImageDiff_LIN;
      s->iiEval = &integralImageEval_LIN;
   } else if (strcmp(method, "AFF") == 0) {
      s->iiDiff = &integralImageDiff_AFF;
      s->iiEval = &integralImageEval_AFF;
   } else if (strcmp(method, "NCC") == 0) {
      s->iiDiff = &integralImageDiff_NCC;
      s->iiEval = &integralImageEval_NCC;
   } else if (strcmp(method, "ZSSDc") == 0) {
      s->iiDiff = &integralImageDiff_ZSSDc;
      s->iiEval = &integralImageEval_ZSSDc;
   } else if (strcmp(method, "NCCc") == 0) {
      s->iiDiff = &integralImageDiff_NCCc;
      s->iiEval = &integralImageEval_NCCc;
   } else
      // safety exit in case of wrong method string
      error("Wrong method in  setup_BlockMatchingStuff()");
}

/** clean the precomputed images stored in t_struct_BlockMatchingStuff */
static void clean_BlockMatchingStuff(struct t_struct_BlockMatchingStuff *s)
{
   free(s->mu1);
   free(s->mu2);
   free(s->std1);
   free(s->std2);
   free(s->norm1);
   free(s->norm2);
   free(s->BTmax1);
   free(s->BTmax2);
   free(s->BTmin1);
   free(s->BTmin2);
   free(s->diff);
}



/* pixel stereo correlation using integral image*/
/* (w,h) is the size of the window.
 * if odd then the window is symmetric with respect
 * to the center O  (ie 5:  xxOxx)
 * it even the window is asymmetric (ie 4: xxOx) */
void stereoInt(Fimage * in1,	// left image
	       Fimage * in2,	// right image
	       Fimage * disp,	// output disparities/offsets prealloc. 2ch
	       Fimage * cost,	// output cost prealloc. 1ch
	       Fimage * dispr,	// right-to-left offsets prealloc. 2ch
	       Fimage * costr,	// right-to-left costs prealloc. 1ch
	       int patch_sz[],	// patch size [w,h]
	       int hor_range[],	// horizontal offset range [min,max]
	       int ver_range[],	// vertical offset range [min,max]
	       char *method)	// name of the method SSD,SAD,ZSSD,SSD/Norm,NCC..
{
   int nc1 = in1->ncol, nr1 = in1->nrow, nch = in1->nch;
   int nc2 = in2->ncol, nr2 = in2->nrow;
   int w = patch_sz[0];
   int h = patch_sz[1];
   int halfh = h / 2;
   int halfw = w / 2;

   printf(".");
   fflush(stdout);

   // This Structure contains precomputed variables and status
   // needed for the integral image computation and evaluation
   struct t_struct_BlockMatchingStuff BMStuff;

   // Precompute the variables based on the chosen method
   // and store them in BMStuff
   setup_BlockMatchingStuff(&BMStuff, in1, in2, patch_sz, method);

   // Function pointers to the methods for computing the
   // integral image differences and evaluating the integral images
   t_func_integralImageDiff iiDiff = BMStuff.iiDiff;
   t_func_integralImageEval iiEval = BMStuff.iiEval;

   /* Initialize the correlation vectors */
   for (int i = 0; i < cost->N; i++)
      _vv(cost, i) = INFINITY;
   if (costr)
      for (int i = 0; i < costr->N; i++)
	 _vv(costr, i) = INFINITY;

   Fimage *tmpCosts = malloc_fimage2(nc1, nr1);
   for (int off_y = ver_range[0]; off_y <= ver_range[1]; off_y++) {
      for (int off_x = hor_range[0]; off_x <= hor_range[1]; off_x++) {
	 // temporary cost vector

	 // this function computes the difference of in1
	 // with the translated in2 by [off_x,off_y]
	 // and stores the result in BMStuff.diff.
	 // Also updates in BMStuff the values of off_x/y
	 (*iiDiff) (in1, in2, off_x, off_y, &BMStuff);

	 /* compute integral image */
	 struct iimage_t *ii =
	     integralImage(BMStuff.diff, nc1, nr1, BMStuff.numdiff);

	 /* compute costs for the current offsets (ie SSD, SAD, or correlation)
	  * Costs are stored at the position corresponding to the patch's center.
	  * For odd sized patches the center is obvious, for even sized patches,
	  * the center is shifted to the right:     x x O x      */
	 for (int y = halfh; y <= nr1 - h + halfh; y++)	/* explore only positions */
	    for (int x = halfw; x <= nc1 - w + halfw; x++) {	/* with valid disparities */
	       _vpla(tmpCosts, x, y, 0) =
		   (*iiEval) (ii, x - halfw, y - halfh, w, h, nch,
			      &BMStuff);
	    }

	 /* free the integral image for this offset */
	 free(ii);

	 /* keep the best cost from tmpCosts and corresponding offset */
	 for (int y = halfh; y <= nr1 - h + halfh; y++)	/* explore only positions */
	    for (int x = halfw; x <= nc1 - w + halfw; x++)	/* with valid disparities */
	       if (_vpla(cost, x, y, 0) > _vpla(tmpCosts, x, y, 0)) {
		  _vpla(cost, x, y, 0) = _vpla(tmpCosts, x, y, 0);
		  _vpla(disp, x, y, 0) = off_x;
		  _vpla(disp, x, y, 1) = off_y;
	       }

	 /* keep the best costs for the right image */
	 if (dispr && costr)
	    for (int y = halfh; y <= nr1 - h + halfh; y++)	/* explore only positions */
	       for (int x = halfw; x <= nc1 - w + halfw; x++) {	/* with valid disparities */
		  // is the patch center within the bounds of in2?
		  if ((x + off_x) >= halfw
		      && (x + off_x) <= nc2 - w + halfw
		      && (y + off_y) >= halfh
		      && (y + off_y) <= nr2 - h + halfh) {
		     if (_vpla(costr, x + off_x, y + off_y, 0) >
			 _vpla(tmpCosts, x, y, 0)) {
			_vpla(costr, x + off_x, y + off_y, 0) =
			    _vpla(tmpCosts, x, y, 0);
			_vpla(dispr, x + off_x, y + off_y, 0) = -off_x;
			_vpla(dispr, x + off_x, y + off_y, 1) = -off_y;
		     }
		  }
	       }
      }
   }
   free(tmpCosts);

   clean_BlockMatchingStuff(&BMStuff);
}
