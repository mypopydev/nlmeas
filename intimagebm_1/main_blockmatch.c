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
#include "iio.h"
#include "smartparameter.h"
#include "fimage2.h"

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



/***************************************************/


/** IIO wrap. Read the file content into a new multi-channel Fimage
 * returns the allocated Fimage*/
Fimage *fimageread(char *name)
{
   int nx, ny, nc;
   float *fdisp = iio_read_image_float_split(name, &nx, &ny, &nc);
   if (NULL == fdisp) {
      fprintf(stderr, "unreadable input files\n");
      exit(EXIT_FAILURE);
   }
   Fimage *res = malloc_fimage3pla(nx, ny, nc);
   for (int i = 0; i < res->N; i++)
      res->v[i] = fdisp[i];
   free(fdisp);
   return res;
}

/** IIO wrap. Write a file with the Fimage */
void fimagewrite(char *name, Fimage * i)
{
   if (i->planar == 1)
      iio_save_image_float_split(name, i->v, i->ncol, i->nrow, i->nch);
   else
      iio_save_image_float_vec(name, i->v, i->ncol, i->nrow, i->nch);
}

void removeInfNan(Fimage * im)
{
   int i;
   for (i = 0; i < im->N; i++)
      if (isnan(_vv(im, i)) || isinf(_vv(im, i)))
	 _vv(im, i) = 0;
}


void help(char **v)
{
   fprintf(stderr, "Block matching (2D offsets) using integral image\n");
   fprintf(stderr,
	   "usage: %s [-t method] [-w window_width] [-h window_height] \n\
[-r horizontal_min] [-R horizontal_max] [-s vertical_min] [-S vertical_max] \n\
leftIM rightIM disparity mincost [disparity_Right mincost_Right]\n", v[0]);
   fprintf(stderr, "method can be: SSD sum of squaerd differences\n");
   fprintf(stderr, "               SAD sum of absolute differences\n");
   fprintf(stderr, "               ZSSD zero mean SSD\n");
   fprintf(stderr, "               SSDNorm L^2-normalized SSD\n");
   fprintf(stderr, "               NCC normalized cross correlation\n");
   fprintf(stderr, "               AFF “affine” similarity measure\n");
   fprintf(stderr,
	   "               LIN linearized “affine” similarity measure\n");
   fprintf(stderr,
	   "               BT2DSSD SSD with Birchfield-Tomasi 2D cost\n");
   fprintf(stderr,
	   "               BT2DSAD SAD with Birchfield-Tomasi 2D cost\n");
// these are not ment for 2d block matching, only consider 1d sampling artifacts
//   fprintf(stderr,"               BTSSD SSD with Birchfield-Tomasi cost\n");
//   fprintf(stderr,"               BTSAD SAD with Birchfield-Tomasi cost\n");

   exit(1);
}



// c: pointer to original argc
// v: pointer to original argv
// o: option name after hyphen
// d: default value (if NULL, the option takes no argument)
static char *pick_option(int *c, char ***v, char *o, char *d)
{
   int argc = *c;
   char **argv = *v;
   int id = d ? 1 : 0;
   for (int i = 0; i < argc - id; i++)
      if (argv[i][0] == '-' && 0 == strcmp(argv[i] + 1, o)) {
	 char *r = argv[i + id] + 1 - id;
	 *c -= id + 1;
	 for (int j = i; j < argc - id; j++)
	    (*v)[j] = (*v)[j + id + 1];
	 return r;
      }
   return d;
}

int main(int argc, char **argv)
{
   int patch_sz[2], hor_range[2], ver_range[2];
   patch_sz[0] = atoi(pick_option(&argc, &argv, "w", "7"));
   patch_sz[1] = atoi(pick_option(&argc, &argv, "h", "7"));
   hor_range[0] = atoi(pick_option(&argc, &argv, "r", "-5"));
   hor_range[1] = atoi(pick_option(&argc, &argv, "R", "5"));
   ver_range[0] = atoi(pick_option(&argc, &argv, "s", "-5"));
   ver_range[1] = atoi(pick_option(&argc, &argv, "S", "5"));
   char *method = pick_option(&argc, &argv, "t", "ZSSD");
   if (argc < 5) {
      fprintf(stderr, "too few parameters\n");
      help(argv);
      return 1;
   }
   char *in1_file = argv[1];
   char *in2_file = argv[2];
   char *outdisp_file = argv[3];
   char *outcost_file = argv[4];
   char *outdispR_file = argc > 5 ? argv[5] : NULL;
   char *outcostR_file = argc > 6 ? argv[6] : NULL;

   /* program call */

   Fimage *in1 = fimageread(in1_file);
   Fimage *in2 = fimageread(in2_file);
   removeInfNan(in1);
   removeInfNan(in2);
   int nx = in1->ncol;
   int ny = in1->nrow;
   int nx2 = in2->ncol;
   int ny2 = in2->nrow;
   Fimage *disp = malloc_fimage3pla(nx, ny, 2);
   Fimage *cost = malloc_fimage3pla(nx, ny, 1);
   Fimage *dispR = malloc_fimage3pla(nx2, ny2, 2);
   Fimage *costR = malloc_fimage3pla(nx2, ny2, 1);

   // sequential
//   stereoInt(in1, in2, disp, cost, dispR, costR,
//                  patch_sz, hor_range, ver_range, method);
   // parallel
   bmIntSliceSubpix(in1, in2, disp, cost, dispR, costR,
		    patch_sz, hor_range, ver_range, method);

   fimagewrite(outdisp_file, disp);
   fimagewrite(outcost_file, cost);
   if (outcostR_file) {
      fimagewrite(outdispR_file, dispR);
      fimagewrite(outcostR_file, costR);
   }

   /* free memory (the program ends here) */
   free(in1);
   free(in2);
   free(disp);
   free(cost);
   free(dispR);
   free(costR);

   return 0;
}
