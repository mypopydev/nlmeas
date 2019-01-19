// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright 2013, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
// All rights reserved.

/* IN THIS FILE: 
 *
 * integralImageDiff_SAD
 * integralImageDiff_SSD
 * integralImageDiff_BTSSD
 * integralImageDiff_BTSAD
 * integralImageDiff_ZSSDc
 * integralImageDiff_ZSSD
 * integralImageDiff_PROD
 * integralImageDiff_SSDNorm
 * integralImageDiff_LIN
 * integralImageDiff_AFF
 * integralImageDiff_NCCc
 * integralImageDiff_NCC
 *
 *
 *  integralImageEval_SSD   // Eval sum of sq. diff. also: SAD, BTSSD & BTSAD
 *  integralImageEval_ZSSD  // Eval SSD - mean
 *  integralImageEval_SSDNorm // Eval SSD/Norm
 *  integralImageEval_LIN   // Eval LIN
 *  integralImageEval_AFF   // Eval AFF
 *  integralImageEval_NCC   // Eval normalized cross correlation
 *  integralImageEval_ZSSDc // channel by channel versions of ZSSD AND NCC
 *  integralImageEval_NCCc
 *
 * */

#include <assert.h>

// mute explicitly unused parameters
#define UNUSED(x) (void)(x)


/**************************************************************************
  SAD AND SSD
 **************************************************************************/


/** generates a vector containing the pixel-wise AD (absolute differences)
 * of in1 and in2, by shifting in2 (dx,dy)
 * The output vector is allocated and stored in the structure s, 
 * it has a single channel, that is: for color images the pixels are added as
 * |d[r]| + |d[g]| + |d[b]| 
 * This function is used for generating the integral image for SAD */
static void integralImageDiff_SAD(Fimage * in1, Fimage * in2,
				  int dx, int dy,
				  struct t_struct_BlockMatchingStuff *s)
{
   assert(in1->nch == in2->nch);	// check that nothing has gone horribly wrong

   int npix = in1->ncol * in1->nrow;

   // update offset in the integralImageStructure
   s->offx = dx;
   s->offy = dy;

   /** allocate the difference vector in the structure s (if needed) */
   iifloat_t *diffVector = s->diff;
   if (diffVector == NULL)
      diffVector = malloc(npix * sizeof *diffVector);
   assert(NULL != diffVector);
   s->diff = diffVector;
   s->numdiff = 1;

   /** translate in2 by (dx,dy) symmetrizing at the boundaries */
   /** and compute the differences */
   for (int c = 0; c < in1->nch; c++)
      for (int y = 0; y < in1->nrow; y++)
	 for (int x = 0; x < in1->ncol; x++) {
	    long double tmp =
		_vpla(in1, x, y, c) - _vsympla(in2, x + dx, y + dy, c);
	    int pos = x + y * in1->ncol;
	    if (c == 0)
	       diffVector[pos] = fabsl(tmp);	// overwrite
	    else
	       diffVector[pos] += fabsl(tmp);	// accumulate
	 }
}



/** generates a vector containing the pixel-wise SD (squared differences)
 * of in1 and in2, where in2 have been shifted by (dx,dy)
 * The output vector is allocated and stored in the structure s, 
 * it has a single channel, that is: for color images the pixels are added as
 * d[r]^2 + d[g]^2 + d[b]^2 
 * This function is used for generating the integral image for SSD and ZSSD*/
static void integralImageDiff_SSD(Fimage * in1, Fimage * in2,
				  int dx, int dy,
				  struct t_struct_BlockMatchingStuff *s)
{
   assert(in1->nch == in2->nch);	// check that nothing has gone horribly wrong

   int npix = in1->ncol * in1->nrow;

   // update offset in the integralImageStructure
   s->offx = dx;
   s->offy = dy;

   /** allocate the difference vector in the structure s (if needed) */
   iifloat_t *diffVector = s->diff;
   if (diffVector == NULL)
      diffVector = malloc(npix * sizeof *diffVector);
   assert(NULL != diffVector);
   s->diff = diffVector;
   s->numdiff = 1;

   /** translate in2 by (dx,dy) symmetrizing at the boundaries */
   /** and compute the differences */
   for (int c = 0; c < in1->nch; c++)
      for (int y = 0; y < in1->nrow; y++)
	 for (int x = 0; x < in1->ncol; x++) {
	    long double tmp =
		_vpla(in1, x, y, c) - _vsympla(in2, x + dx, y + dy, c);
	    int pos = x + y * in1->ncol;
	    if (c == 0)
	       diffVector[pos] = tmp * tmp;	// overwrite
	    else
	       diffVector[pos] += tmp * tmp;	// accumulate
	 }
}


// is the same as  evalIntImageSAD
// is the same as  evalIntImageSSDBT
// is the same as  evalIntImageSADBT
static long double integralImageEval_SSD(struct iimage_t *ii, int x, int y,
					 int w, int h, int nch,
					 struct t_struct_BlockMatchingStuff
					 *s)
{
   int xx = x + s->offx;
   int yy = y + s->offy;
   long double v = 1000000;	// TODO should be INFINITY, but is inefficient
   // check that (x,y)+off falls inside in2 otherwise return a big number
   if (xx >= 0 && yy >= 0 && xx + w - 1 < s->in2->ncol
       && yy + h - 1 < s->in2->nrow)
      v = evaluateIntegralImage(ii, 0, x, y, w, h);
   assert(v >= -0.001);		// check that nothing has gone horribly wrong
   return fmax(v / (w * h * nch), 0);
}




/**************************************************************************
  BIRCHFIELD AND TOMMASI SAMPLING INSENSITIVE DISTANCE BTSAD AND BTSSD
 **************************************************************************/

/** compute the min and max of the bilinearly interpolated 
 * input image in an neighborhood of [-1/2,1/2] around each point.
 * This reduces to evaluating the max and min of the 
 * interpolated function at three locations  */
void malloc_and_compute_BTminmax(Fimage * in, Fimage ** outmin,
				 Fimage ** outmax)
{
   int nc = in->ncol;
   int nr = in->nrow;
   int nch = in->nch;
   if (*outmin == NULL)
      *outmin = malloc_fimage3pla(nc, nr, nch);
   if (*outmax == NULL)
      *outmax = malloc_fimage3pla(nc, nr, nch);

#define min3(a,b,c) (((a)<(b))? (((a)<(c))?(a):(c)) : (((c)<(b))?(c):(b)) )
#define max3(a,b,c) (((a)>(b))? (((a)>(c))?(a):(c)) : (((c)>(b))?(c):(b)) )

   for (int c = 0; c < nch; c++)
      for (int j = 0; j < nr; j++)
	 for (int i = 0; i < nc; i++) {
	    long double I1 = _vpla(in, i, j, c);
	    long double I1p = (_vsympla(in, i + 1, j, c) + I1) / 2;
	    long double I1m = (_vsympla(in, i - 1, j, c) + I1) / 2;
	    long double I1min = min3(I1p, I1m, I1);
	    long double I1max = max3(I1p, I1m, I1);
	    _vpla(*outmin, i, j, c) = I1min;
	    _vpla(*outmax, i, j, c) = I1max;
	 }
}




/* generates a vector containing the SD (squared differences) 
 * of in1 and in2, by shifting in2 (dx,dy). 
 * The output vector have a single channel, 
 * that is: for color images the pixels are added d[r]^2 + d[g]^2 + d[b]^2 
 * This function is used for generating the integral image for BTSSD*/
static void integralImageDiff_BTSSD(Fimage * in1, Fimage * in2,
				    int dx, int dy,
				    struct t_struct_BlockMatchingStuff *s)
{
   assert(in1->nch == in2->nch);	// check that nothing has gone horribly wrong

   // check if the precomputed min and max are there otherwise compute them
   if (s->BTmin1 == NULL || s->BTmax1 == NULL)
      malloc_and_compute_BTminmax(in1, &s->BTmin1, &s->BTmax1);

   if (s->BTmin2 == NULL || s->BTmax2 == NULL)
      malloc_and_compute_BTminmax(in2, &s->BTmin2, &s->BTmax2);

   int npix = in1->ncol * in1->nrow;
   int nch = in1->nch;
   int nc = in1->ncol;
   int nr = in1->nrow;

   // update offset in the integralImageStructure
   s->offx = dx;
   s->offy = dy;

   /** allocate the difference vector in the structure s (if needed) */
   iifloat_t *diffVector = s->diff;
   if (diffVector == NULL)
      diffVector = malloc(npix * sizeof *diffVector);
   assert(NULL != diffVector);
   s->diff = diffVector;
   s->numdiff = 1;

   /** compute the differences */
   for (int c = 0; c < nch; c++)
      for (int y = 0; y < nr; y++)
	 for (int x = 0; x < nc; x++) {
	    int pos1 = _ppla(in1, x, y, c);
	    int pos2 = _pinpla(in2, x + s->offx, y + s->offy, c);
	    long double tmp = 0;

	    if (pos2 >= 0) {
	       long double I1, I1max, I1min, I2, I2max, I2min, tmp1, tmp2;
	       I1 = _vv(in1, pos1);
	       I1max = _vv(s->BTmax1, pos1);
	       I1min = _vv(s->BTmin1, pos1);

	       I2 = _vv(in2, pos2);
	       I2max = _vv(s->BTmax2, pos2);
	       I2min = _vv(s->BTmin2, pos2);

	       tmp1 = max3(0, I2 - I1max, I1min - I2);
	       tmp2 = max3(0, I1 - I2max, I2min - I1);
	       tmp = fmin(tmp1, tmp2);
	    } else {
	       tmp = _vv(in1, pos1);	// TODO, this should be handled by EVAL_SSD
	    }

	    if (c == 0)
	       diffVector[x + y * nc] = tmp * tmp;	// overwrite
	    else
	       diffVector[x + y * nc] += tmp * tmp;	// accumulate
	 }
}



/* generates a vector containing the AD (absolute differences) 
 * of in1 and in2, by shifting in2 (dx,dy). 
 * The output vector have a single channel, 
 * that is: for color images the pixels are added |d[r]| + |d[g]| + |d[b]| 
 * This function is used for generating the integral image for BTSAD*/
static void integralImageDiff_BTSAD(Fimage * in1, Fimage * in2,
				    int dx, int dy,
				    struct t_struct_BlockMatchingStuff *s)
{
   assert(in1->nch == in2->nch);	// check that nothing has gone horribly wrong

   // check if the precomputed min and max are there otherwise compute them
   if (s->BTmin1 == NULL || s->BTmax1 == NULL)
      malloc_and_compute_BTminmax(in1, &s->BTmin1, &s->BTmax1);

   if (s->BTmin2 == NULL || s->BTmax2 == NULL)
      malloc_and_compute_BTminmax(in2, &s->BTmin2, &s->BTmax2);

   int npix = in1->ncol * in1->nrow;
   int nch = in1->nch;
   int nc = in1->ncol;
   int nr = in1->nrow;

   // update offset in the integralImageStructure
   s->offx = dx;
   s->offy = dy;

   /** allocate the difference vector in the structure s (if needed) */
   iifloat_t *diffVector = s->diff;
   if (diffVector == NULL)
      diffVector = malloc(npix * sizeof *diffVector);
   assert(NULL != diffVector);
   s->diff = diffVector;
   s->numdiff = 1;

   /** compute the differences */
   for (int c = 0; c < nch; c++)
      for (int y = 0; y < nr; y++)
	 for (int x = 0; x < nc; x++) {
	    int pos1 = _ppla(in1, x, y, c);
	    int pos2 = _pinpla(in2, x + s->offx, y + s->offy, c);
	    long double tmp = 0;

	    if (pos2 >= 0) {
	       long double I1, I1max, I1min, I2, I2max, I2min, tmp1, tmp2;
	       I1 = _vv(in1, pos1);
	       I1max = _vv(s->BTmax1, pos1);
	       I1min = _vv(s->BTmin1, pos1);

	       I2 = _vv(in2, pos2);
	       I2max = _vv(s->BTmax2, pos2);
	       I2min = _vv(s->BTmin2, pos2);

	       tmp1 = max3(0, I2 - I1max, I1min - I2);
	       tmp2 = max3(0, I1 - I2max, I2min - I1);
	       tmp = fmin(tmp1, tmp2);
	    } else {
	       tmp = _vv(in1, pos1);	// TODO, this should be handled by EVAL_SSD
	    }

	    if (c == 0)
	       diffVector[x + y * nc] = tmp > 0 ? tmp : -tmp;	// overwrite
	    else
	       diffVector[x + y * nc] += tmp > 0 ? tmp : -tmp;	// accumulate
	 }
}


static float evaluate_bilinear(float a, float b, float c, float d,
			       float x, float y)
{
   float r = 0;
   r += a * (1 - x) * (1 - y);
   r += b * (x) * (1 - y);
   r += c * (1 - x) * (y);
   r += d * (x) * (y);
   return r;
}


/** compute the min and max of the bilinearly interpolated 
 * input image in an neighborhood of [-1/2,1/2]^2 around each point.
 * This reduces to evaluating the max and min of the 
 * interpolated function at nine locations  */
void malloc_and_compute_BT2Dminmax(Fimage * in, Fimage ** outmin,
				   Fimage ** outmax)
{
   int nc = in->ncol;
   int nr = in->nrow;
   int nch = in->nch;
   if (*outmin == NULL)
      *outmin = malloc_fimage3pla(nc, nr, nch);
   if (*outmax == NULL)
      *outmax = malloc_fimage3pla(nc, nr, nch);

   for (int ch = 0; ch < nch; ch++)
      for (int j = 0; j < nr; j++)
	 for (int i = 0; i < nc; i++) {
	    float I1 = _vpla(in, i, j, ch);
	    float I1min = I1;
	    float I1max = I1;
	    for (int l = -1; l <= 1; l++)
	       for (int k = -1; k <= 1; k++) {
		  int dx = (int) floor(((float) k) / 2.0);
		  int dy = (int) floor(((float) l) / 2.0);
		  // symmetric bondaries act as Neumann boundaries
		  float a, b, c, d, Icur;
		  a = _vsympla(in, i + dx, j + dy, ch);
		  b = _vsympla(in, i + dx + 1, j + dy, ch);
		  c = _vsympla(in, i + dx, j + dy + 1, ch);
		  d = _vsympla(in, i + dx + 1, j + dy + 1, ch);
		  Icur =
		      evaluate_bilinear(a, b, c, d, k / 2.0 - dx,
					l / 2.0 - dy);

		  I1min = fmin(I1min, Icur);
		  I1max = fmax(I1max, Icur);
	       }
	    _vpla(*outmin, i, j, ch) = I1min;
	    _vpla(*outmax, i, j, ch) = I1max;
	 }
}

/* generates a vector containing the AD (absolute differences) 
 * of in1 and in2, by shifting in2 (dx,dy). 
 * The output vector have a single channel, 
 * that is: for color images the pixels are added |d[r]| + |d[g]| + |d[b]| 
 * This function is used for generating the integral image for BT2DSAD
 * it intercepts the generation of the BTmax? and BTmin? variables 
 * then calls the BTSAD */
static void integralImageDiff_BT2DSAD(Fimage * in1, Fimage * in2,
				      int dx, int dy,
				      struct t_struct_BlockMatchingStuff
				      *s)
{
   // check if the precomputed min and max are there otherwise compute them
   if (s->BTmin1 == NULL || s->BTmax1 == NULL)
      malloc_and_compute_BT2Dminmax(in1, &s->BTmin1, &s->BTmax1);

   if (s->BTmin2 == NULL || s->BTmax2 == NULL)
      malloc_and_compute_BT2Dminmax(in2, &s->BTmin2, &s->BTmax2);

   integralImageDiff_BTSAD(in1, in2, dx, dy, s);

}

static void integralImageDiff_BT2DSSD(Fimage * in1, Fimage * in2,
				      int dx, int dy,
				      struct t_struct_BlockMatchingStuff
				      *s)
{
   // check if the precomputed min and max are there otherwise compute them
   if (s->BTmin1 == NULL || s->BTmax1 == NULL)
      malloc_and_compute_BT2Dminmax(in1, &s->BTmin1, &s->BTmax1);

   if (s->BTmin2 == NULL || s->BTmax2 == NULL)
      malloc_and_compute_BT2Dminmax(in2, &s->BTmin2, &s->BTmax2);

   integralImageDiff_BTSSD(in1, in2, dx, dy, s);

}


/**************************************************************************
  ZSSD or SSD-mean  
 **************************************************************************/


/** compute the mean of IN over patches of w*h*nch pixels, 
 * stores it in OUT at the coordinate of the upper left corner of the 
 * patch (for consistency with the integral image evaluation coordinate 
 * system). OUT is allocated of size IN->ncol*IN->nrow*/
Fimage *malloc_and_compute_mu(Fimage * in, int patch_sz[])
{
   int w = patch_sz[0];
   int h = patch_sz[1];
   int nc = in->ncol;
   int nr = in->nrow;
   int nch = in->nch;
   Fimage *out = malloc_fimage3pla(nc, nr, 1);
   long double const norm_factor = 1.0 / (w * h * nch);

   for (int j = 0; j < nr - h + 1; j++)
      for (int i = 0; i < nc - w + 1; i++) {
	 long double mean = 0;
	 // compute mean over a patch centered at i,j
	 for (int c = 0; c < nch; c++)
	    for (int ox = 0; ox < w; ox++)
	       for (int oy = 0; oy < h; oy++)
		  mean += _vpla(in, i + ox, j + oy, c);
	 _vpla(out, i, j, 0) = mean * norm_factor;
      }
   return out;
}


/* generates a vector containing the SD (squared differences)
 * the output vector have a single channel, 
 * that is: for color images the pixels are added d[r]^2 + d[g]^2 + d[b]^2 
 * This function is used for generating the integral image for ZSSD*/
static void integralImageDiff_ZSSD(Fimage * in1, Fimage * in2,
				   int dx, int dy,
				   struct t_struct_BlockMatchingStuff *s)
{
   // check if the precomputed means are there otherwise compute them
   if (s->mu1 == NULL)
      s->mu1 = malloc_and_compute_mu(in1, s->patch_sz);
   if (s->mu2 == NULL)
      s->mu2 = malloc_and_compute_mu(in2, s->patch_sz);

   // call to the other SSD
   integralImageDiff_SSD(in1, in2, dx, dy, s);
}



// computes ZSSD(p,q) = SSD(p,q) - ( N *(mu_1(p) - mu_2(q) )^2 )
static long double integralImageEval_ZSSD(struct iimage_t *ii, int x,
					  int y, int w, int h, int nch,
					  struct
					  t_struct_BlockMatchingStuff *s)
{
   long double N = w * h * nch;
   /** evaluate the integral image */
   long double v = evaluateIntegralImage(ii, 0, x, y, w, h);

   /** compute the patch origin of im2, by compensating the offset */
   int xx = x + s->offx;
   int yy = y + s->offy;

   // accumulate mu1 and mu2 
   assert(s->mu1 && s->mu2);
   int pidx1 = _ppla(s->mu1, x, y, 0);
   int pidx2 = _pinpla(s->mu2, xx, yy, 0);
   // check if pidx2 is inside the valid region
   if (xx < s->mu2->ncol - w + 1 && yy < s->mu2->nrow - h + 1
       && pidx2 >= 0) {
      long double tmp = _vv(s->mu1, pidx1) - _vv(s->mu2, pidx2);
      v -= N * tmp * tmp;
   } else {
      v += _vv(s->mu1, pidx1);	// TODO should be INFINITY, but is inefficient
   }

   //assert( v >=-0.1);  // check that nothing has gone horribly wrong
   return fmax(v / N, 0);
}







/**************************************************************************
  SSD/Norm
 **************************************************************************/



/** Computes the L2-norm of patches of size w*h*nch from the image IN,
 * stores it in OUT at the coordinate of the upper left corner of the 
 * patch (for consistency with the integral image evaluation coordinate 
 * system). OUT is allocated of size IN->ncol*IN->nrow 
 * For multi-channel images the nch*w*h pixels of each patch are 
 * considered as a single vector and normalized altogether.*/
Fimage *malloc_and_compute_norm(Fimage * in, int patch_sz[])
{
   int w = patch_sz[0];
   int h = patch_sz[1];
   int nc = in->ncol;
   int nr = in->nrow;
   int nch = in->nch;
   Fimage *out = malloc_fimage3pla(nc, nr, 1);

   for (int j = 0; j < nr - h + 1; j++)
      for (int i = 0; i < nc - w + 1; i++) {
	 long double norm = 0;
	 // compute L2-norm of the color patch 
	 for (int c = 0; c < nch; c++)
	    for (int ox = 0; ox < w; ox++)
	       for (int oy = 0; oy < h; oy++) {
		  long double tmp = _vpla(in, i + ox, j + oy, c);
		  norm += tmp * tmp;
	       }
	 _vpla(out, i, j, 0) = sqrtl(norm);
      }
   return out;
}



/******************* SSD/Norm ***************************/
/* generates a vector with pixel-wise products 
 * needed for computing the Covariances and Correlations 
 * and that is used by SSD/Norm:
 *    min sum (u /||u|| - v /||v||)^2 
 *    or equivalently min  1 - sum u(x) v(x) / (||u|| ||v|| )
 * */
static void integralImageDiff_PROD(Fimage * in1, Fimage * in2,
				   int dx, int dy,
				   struct t_struct_BlockMatchingStuff *s)
{
   assert(in1->nch == in2->nch);	// check that nothing has gone horribly wrong

   int npix = in1->ncol * in1->nrow;

   // update offset in the integralImageStructure
   s->offx = dx;
   s->offy = dy;

   /** allocate the products vector in the structure s (if needed) */
   iifloat_t *diffVector = s->diff;
   if (diffVector == NULL)
      diffVector = malloc(npix * sizeof *diffVector);
   assert(NULL != diffVector);
   s->diff = diffVector;
   s->numdiff = 1;

   /** translate in2 by (dx,dy) symmetrizing at the boundaries */
   /** and compute the pixel-wise product */
   for (int c = 0; c < in1->nch; c++)
      for (int y = 0; y < in1->nrow; y++)
	 for (int x = 0; x < in1->ncol; x++) {
	    long double tmp =
		_vpla(in1, x, y, c) * _vsympla(in2, x + dx, y + dy, c);
	    int pos = x + y * in1->ncol;
	    if (c == 0)
	       diffVector[pos] = tmp;	// overwrite
	    else
	       diffVector[pos] += tmp;	// accumulate
	 }
}



static void integralImageDiff_SSDNorm(Fimage * in1, Fimage * in2,
				      int dx, int dy,
				      struct t_struct_BlockMatchingStuff
				      *s)
{
   // check if the precomputed norm are there otherwise compute them
   if (s->norm1 == NULL)
      s->norm1 = malloc_and_compute_norm(in1, s->patch_sz);
   if (s->norm2 == NULL)
      s->norm2 = malloc_and_compute_norm(in2, s->patch_sz);

   // call to PROD 
   integralImageDiff_PROD(in1, in2, dx, dy, s);
}


// computes SSD/Norm(p,q) = 2 - 2*PROD(p,q)/(||u(p)|| * ||v(q)||)
static long double integralImageEval_SSDNorm(struct iimage_t *ii, int x,
					     int y, int w, int h, int nch,
					     struct
					     t_struct_BlockMatchingStuff
					     *s)
{
   UNUSED(nch);
   /** evaluate the integral image */
   long double v = evaluateIntegralImage(ii, 0, x, y, w, h);

   /** compute the patch origin of im2, by compensating the offset */
   int xx = x + s->offx;
   int yy = y + s->offy;

   long double res = 0;
   // read norm1 and norm2 at the computed locations
   assert(s->norm1 && s->norm2);
   int pidx1 = _ppla(s->norm1, x, y, 0);
   int pidx2 = _pinpla(s->norm2, xx, yy, 0);
   // check if pidx2 is inside the valid region
   if (xx < s->norm2->ncol - w + 1 && yy < s->norm2->nrow - h + 1
       && pidx2 >= 0) {
      long double n1 = _vv(s->norm1, pidx1);
      long double n2 = _vv(s->norm2, pidx2);
      if (n1 > 0 && n2 > 0)
	 res = v / (n1 * n2);
   }
   //   assert( 1-res >= -0.1);  // check that nothing has gone horribly wrong
   return (2. - 2 * fmin(res, 1.));
}


/**************************************************************************
   LIN : Linearized version of the affine similarity measure defined in:
   J. Delon and A. Desolneux. Stabilization of flicker-like effects 
   in image sequences through local contrast correction. 
   SIAM Journal on Imaging Sciences, 3(4):703–734, 2010.
   LIN is from: 
   Adrian Marques, Fast nearest-neighbour searches, Master M2 thesis, ENS-Cachan.
   Advisors: Julie Delon, Andres Almansa and Yann Gousseau September 2010.
 **************************************************************************/

/** Computes the squared L2-norm of patches of size w*h*nch from the image IN,
 * stores it in OUT at the coordinate of the upper left corner of the 
 * patch (for consistency with the integral image evaluation coordinate 
 * system). OUT is allocated of size IN->ncol*IN->nrow
 * For multi-channel images the nch*w*h pixels of each patch are 
 * considered as a single vector and normalized altogether.*/
Fimage *malloc_and_compute_sqnorm(Fimage * in, int patch_sz[])
{
   int w = patch_sz[0];
   int h = patch_sz[1];
   int nc = in->ncol;
   int nr = in->nrow;
   int nch = in->nch;
   Fimage *out = malloc_fimage3pla(nc, nr, 1);

   for (int j = 0; j < nr - h + 1; j++)
      for (int i = 0; i < nc - w + 1; i++) {
	 long double norm = 0;
	 // compute squared L2-norm of the color patch 
	 for (int c = 0; c < nch; c++)
	    for (int ox = 0; ox < w; ox++)
	       for (int oy = 0; oy < h; oy++) {
		  long double tmp = _vpla(in, i + ox, j + oy, c);
		  norm += tmp * tmp;
	       }
	 _vpla(out, i, j, 0) = (norm);
      }
   return out;
}

// the same as SSDNorm, but where the L2norms are squared
static void integralImageDiff_LIN(Fimage * in1, Fimage * in2,
				  int dx, int dy,
				  struct t_struct_BlockMatchingStuff *s)
{
   // check if the precomputed norm are there otherwise compute them
   if (s->norm1 == NULL)
      s->norm1 = malloc_and_compute_sqnorm(in1, s->patch_sz);
   if (s->norm2 == NULL)
      s->norm2 = malloc_and_compute_sqnorm(in2, s->patch_sz);

   // call to PROD
   integralImageDiff_PROD(in1, in2, dx, dy, s);
}

// computes LIN^2(p,q) = max( ||u(p)||^2, ||v(q)||^2) (1 - PROD^2(p,q)/(||u(p)||^2 * ||v(q)||^2)
static long double integralImageEval_LIN(struct iimage_t *ii, int x, int y,
					 int w, int h, int nch,
					 struct t_struct_BlockMatchingStuff
					 *s)
{
   /** evaluate the integral image */
   long double prod = evaluateIntegralImage(ii, 0, x, y, w, h);
   UNUSED(nch);

   /** compute the patch origin of im2, by compensating the offset */
   int xx = x + s->offx;
   int yy = y + s->offy;

   long double res = 0;
   // read norm1 and norm2 at the computed locations
   assert(s->norm1 && s->norm2);
   int pidx1 = _ppla(s->norm1, x, y, 0);
   int pidx2 = _pinpla(s->norm2, xx, yy, 0);
   // check if pidx2 is inside the valid region
   if (xx < s->norm2->ncol - w + 1 && yy < s->norm2->nrow - h + 1
       && pidx2 >= 0) {
      long double n1 = _vv(s->norm1, pidx1);
      long double n2 = _vv(s->norm2, pidx2);
      if (n1 > 0 && n2 > 0)
	 res = fmax(n1, n2) * (1 - (prod * prod) / (n1 * n2));
      else
	 res = fmax(n1, 100000);
   } else {			// otherwise set a large number..
      long double n1 = _vv(s->norm1, pidx1);
      res = fmax(n1, 100000);
   }
   //   assert( 1-res >= -0.1);  // check that nothing has gone horribly wrong
   return (fmax(0, res));
}




/**************************************************************************
  NCC 
 **************************************************************************/


/** compute the standard deviation of IN over patches of w*h*nch pixels, 
 * stores it in OUT at the coordinate of the upper left corner of the 
 * patch (for consistency with the integral image evaluation coordinate 
 * system). OUT is allocated of size IN->ncol*IN->nrow*/
Fimage *malloc_and_compute_std(Fimage * in, int patch_sz[])
{
   int w = patch_sz[0];
   int h = patch_sz[1];
   int nc = in->ncol;
   int nr = in->nrow;
   int nch = in->nch;
   Fimage *out = malloc_fimage3pla(nc, nr, 1);
   long double const norm_factor = 1.0 / (w * h * nch);

   for (int j = 0; j < nr - h + 1; j++)
      for (int i = 0; i < nc - w + 1; i++) {
	 // compute mean over a patch centered at i,j
	 long double mean = 0;
	 for (int c = 0; c < nch; c++)
	    for (int ox = 0; ox < w; ox++)
	       for (int oy = 0; oy < h; oy++)
		  mean += _vpla(in, i + ox, j + oy, c);
	 mean = mean * norm_factor;
	 // compute var over a patch centered at i,j
	 long double var = 0;
	 for (int c = 0; c < nch; c++)
	    for (int ox = 0; ox < w; ox++)
	       for (int oy = 0; oy < h; oy++) {
		  long double tmp = _vpla(in, i + ox, j + oy, c) - mean;
		  var += tmp * tmp;
	       }
	 _vpla(out, i, j, 0) = sqrtl(var * norm_factor);
      }
   return out;
}


static void integralImageDiff_NCC(Fimage * in1, Fimage * in2,
				  int dx, int dy,
				  struct t_struct_BlockMatchingStuff *s)
{
   // check if the precomputed norm are there otherwise compute them
   if (s->mu1 == NULL)
      s->mu1 = malloc_and_compute_mu(in1, s->patch_sz);
   if (s->mu2 == NULL)
      s->mu2 = malloc_and_compute_mu(in2, s->patch_sz);
   if (s->std1 == NULL)
      s->std1 = malloc_and_compute_std(in1, s->patch_sz);
   if (s->std2 == NULL)
      s->std2 = malloc_and_compute_std(in2, s->patch_sz);

   // call to PROD 
   integralImageDiff_PROD(in1, in2, dx, dy, s);
}


// computes NCC(p,q) =  
//       1-  ( 1/N *PROD(p,q) - mu(p) * mu(q) / (std(p) * std(q)) ) 
static long double integralImageEval_NCC(struct iimage_t *ii, int x, int y,
					 int w, int h, int nch,
					 struct t_struct_BlockMatchingStuff
					 *s)
{
   // compute the patch centers and compensate offset in im2
   int xx = x + s->offx;
   int yy = y + s->offy;

   long double res = 0;
   // read mu1, mu2, std1, std2 at the computed patch centers
   assert(s->mu1 && s->mu2 && s->std1 && s->std2);
   long double v = evaluateIntegralImage(ii, 0, x, y, w, h);
   int pidx1 = _ppla(s->mu1, x, y, 0);
   int pidx2 = _pinpla(s->mu2, xx, yy, 0);
   // check if pidx2 is inside the valid region
   if (xx < s->mu2->ncol - w + 1 && yy < s->mu2->nrow - h + 1
       && pidx2 >= 0) {
      long double s1 = _vv(s->std1, pidx1);
      long double s2 = _vv(s->std2, pidx2);
      long double mu1 = _vv(s->mu1, pidx1);
      long double mu2 = _vv(s->mu2, pidx2);
      if (s1 > 0 && s2 > 0)
	 res = (v / (w * h * nch) - mu1 * mu2) / (s1 * s2);
   }
   //   assert( res/nch >=-0.1);  // check that nothing has gone horribly wrong
   return 1. - fmin(res, 1.);
}


/**************************************************************************
   AFF: Affine similarity measure as defined in:
   J. Delon and A. Desolneux. Stabilization of flicker-like effects 
   in image sequences through local contrast correction. 
   SIAM Journal on Imaging Sciences, 3(4):703–734, 2010.
 **************************************************************************/
// The preprocessing is the same as for NCC
static void integralImageDiff_AFF(Fimage * in1, Fimage * in2,
				  int dx, int dy,
				  struct t_struct_BlockMatchingStuff *s)
{
   integralImageDiff_NCC(in1, in2, dx, dy, s);
}

// computes AFF(p,q) =  
//    max(std(p)^2, std(q)^2 ) min(1,  1- CORR(p,q))
//    with CORR(p,q) = ( 1/N *PROD(p,q) - mu(p) * mu(q) / (std(p) * std(q)) ) )
static long double integralImageEval_AFF(struct iimage_t *ii, int x, int y,
					 int w, int h, int nch,
					 struct t_struct_BlockMatchingStuff
					 *s)
{
   // compute the patch centers and compensate offset in im2
   int xx = x + s->offx;
   int yy = y + s->offy;

   long double res = 0;
   // read mu1, mu2, std1, std2 at the computed patch centers
   assert(s->mu1 && s->mu2 && s->std1 && s->std2);
   long double v = evaluateIntegralImage(ii, 0, x, y, w, h);
   int pidx1 = _ppla(s->mu1, x, y, 0);
   int pidx2 = _pinpla(s->mu2, xx, yy, 0);
   // check if pidx2 is inside the valid region
   if (xx < s->mu2->ncol - w + 1 && yy < s->mu2->nrow - h + 1
       && pidx2 >= 0) {
      long double corr = 0;
      long double s1 = _vv(s->std1, pidx1);
      long double s2 = _vv(s->std2, pidx2);
      long double mu1 = _vv(s->mu1, pidx1);
      long double mu2 = _vv(s->mu2, pidx2);
      if (s1 > 0 && s2 > 0)
	 corr = (v / (w * h * nch) - mu1 * mu2) / (s1 * s2);
      else
	 corr = 0;
      res = fmax(s1 * s1, s2 * s2) * fmin(1, 1 - corr * fabsl(corr));
   } else {
      long double s1 = _vv(s->std1, pidx1);
      res = 1 * fmax(s1 * s1, 10000);
   }
   return res;
}


/**************************************************************************
  ZSSDc or SSD-mean-c : channel by channel version of ZSSD
 **************************************************************************/

/** compute the mean for each channel of IN over patches of w*h pixels, 
 * stores it in OUT at the coordinate of the upper left corner of the 
 * patch (for consistency with the integral image evaluation coordinate 
 * system). OUT is allocated of size of IN->ncol*IN->nrow*IN->nch  */
Fimage *malloc_and_compute_muc(Fimage * in, int patch_sz[])
{
   int w = patch_sz[0];
   int h = patch_sz[1];
   int nc = in->ncol;
   int nr = in->nrow;
   int nch = in->nch;
   Fimage *out = malloc_fimage3pla(nc, nr, nch);

   for (int c = 0; c < nch; c++)
      for (int j = 0; j < nr - h + 1; j++)
	 for (int i = 0; i < nc - w + 1; i++) {
	    long double mean = 0;
	    // compute mean over a patch centered at i,j
	    for (int ox = 0; ox < w; ox++)
	       for (int oy = 0; oy < h; oy++)
		  mean += _vpla(in, i + ox, j + oy, c);
	    _vpla(out, i, j, c) = mean / (w * h);
	 }
   return out;
}







/* generates a vector containing the SD (squared differences)
 * the output vector have a single channel, 
 * that is: for color images the pixels are added d[r]^2 + d[g]^2 + d[b]^2 
 * This function is used for generating the integral image for ZSSDc*/
static void integralImageDiff_ZSSDc(Fimage * in1, Fimage * in2,
				    int dx, int dy,
				    struct t_struct_BlockMatchingStuff *s)
{
   // check if the precomputed means are there otherwise compute them
   if (s->mu1 == NULL)
      s->mu1 = malloc_and_compute_muc(in1, s->patch_sz);
   if (s->mu2 == NULL)
      s->mu2 = malloc_and_compute_muc(in2, s->patch_sz);

   // call to the other SSD
   integralImageDiff_SSD(in1, in2, dx, dy, s);
}

// computes ZSSDc(p,q) = 
//       SSD(p,q) - \sum_{c \in nch} ( N *(mu^c(p) - mu^c(q) )^2 )
static long double integralImageEval_ZSSDc(struct iimage_t *ii, int x,
					   int y, int w, int h, int nch,
					   struct
					   t_struct_BlockMatchingStuff *s)
{
   long double N = w * h;
   /** evaluate the integral image */
   long double v = evaluateIntegralImage(ii, 0, x, y, w, h);

   /** compute the patch origin of im2, by compensating the offset */
   int xx = x + s->offx;
   int yy = y + s->offy;

   // accumulate mu1 and mu2 
   assert(s->mu1 && s->mu2);
   for (int c = 0; c < nch; c++) {
      int pidx1 = _ppla(s->mu1, x, y, c);
      int pidx2 = _pinpla(s->mu2, xx, yy, c);
      // check if pidx2 is inside the valid region
      if (xx < s->mu2->ncol - w + 1 && yy < s->mu2->nrow - h + 1
	  && pidx2 >= 0) {
	 long double tmp = _vv(s->mu1, pidx1) - _vv(s->mu2, pidx2);
	 v -= N * tmp * tmp;
      } else {
	 v += _vv(s->mu1, pidx1);	// TODO should be INFINITY, but is inefficient
      }
   }

   //assert( v >=-0.1);  // check that nothing has gone horribly wrong
   return fmax(v / (N * nch), 0);
}



/**************************************************************************
  NCCc channel by channel version of NCC
 **************************************************************************/



/** compute the std.dev. for each channel of IN over patches of w*h pixels, 
 * stores it in OUT at the coordinate of the upper left corner of the 
 * patch (for consistency with the integral image evaluation coordinate 
 * system). OUT is allocated of size IN->ncol*IN->nrow*IN->nch */
Fimage *malloc_and_compute_stdc(Fimage * in, int patch_sz[])
{
   int w = patch_sz[0];
   int h = patch_sz[1];
   int nc = in->ncol;
   int nr = in->nrow;
   int nch = in->nch;
   Fimage *out = malloc_fimage3pla(nc, nr, nch);
   long double const norm_factor = 1.0 / (w * h);

   for (int c = 0; c < nch; c++)
      for (int j = 0; j < nr - h + 1; j++)
	 for (int i = 0; i < nc - w + 1; i++) {
	    // compute mean over a patch centered at i,j
	    long double mean = 0;
	    for (int ox = 0; ox < w; ox++)
	       for (int oy = 0; oy < h; oy++)
		  mean += _vpla(in, i + ox, j + oy, c);
	    mean = mean * norm_factor;
	    // compute var over a patch centered at i,j
	    long double var = 0;
	    for (int ox = 0; ox < w; ox++)
	       for (int oy = 0; oy < h; oy++) {
		  long double tmp = _vpla(in, i + ox, j + oy, c) - mean;
		  var += tmp * tmp;
	       }
	    _vpla(out, i, j, c) = sqrtl(var * norm_factor);
	 }
   return out;
}


static void integralImageDiff_NCCc(Fimage * in1, Fimage * in2,
				   int dx, int dy,
				   struct t_struct_BlockMatchingStuff *s)
{
   // check if the precomputed norm are there otherwise compute them
   if (s->mu1 == NULL)
      s->mu1 = malloc_and_compute_muc(in1, s->patch_sz);
   if (s->mu2 == NULL)
      s->mu2 = malloc_and_compute_muc(in2, s->patch_sz);
   if (s->std1 == NULL)
      s->std1 = malloc_and_compute_stdc(in1, s->patch_sz);
   if (s->std2 == NULL)
      s->std2 = malloc_and_compute_stdc(in2, s->patch_sz);

   // COMPUTE THE channel by channel PRODUCTS
   // same as PROD but computed over several channels
   assert(in1->nch == in2->nch);	// check that nothing has gone horribly wrong

   int npix = in1->ncol * in1->nrow;
   int nch = in1->nch;

   // update offset in the integralImageStructure
   s->offx = dx;
   s->offy = dy;

   // allocate the difference vector if needed
   iifloat_t *diffVector = s->diff;
   if (diffVector == NULL)
      diffVector = malloc(npix * nch * sizeof *diffVector);
   assert(NULL != diffVector);
   s->diff = diffVector;
   s->numdiff = nch;

   /** translate in2 by (dx,dy) symmetrizing at the boundaries */
   /** and compute the pixel-wise product*/
   for (int c = 0; c < in1->nch; c++)
      for (int y = 0; y < in1->nrow; y++)
	 for (int x = 0; x < in1->ncol; x++) {
	    long double tmp =
		_vpla(in1, x, y, c) * _vsympla(in2, x + dx, y + dy, c);
	    int pos = x + y * in1->ncol + npix * c;
	    diffVector[pos] = tmp;
	 }
}


// computes: NCCc(p,q) = \sum_{c \in nch}  
//     1-  ( 1/N *PROD(p,q) - mu^c(p) * mu^c(q) / (std^c(p) * std^c(q)) ) 
static long double integralImageEval_NCCc(struct iimage_t *ii, int x,
					  int y, int w, int h, int nch,
					  struct
					  t_struct_BlockMatchingStuff *s)
{
   // compute the patch centers and compensate offset in im2
   int xx = x + s->offx;
   int yy = y + s->offy;

   long double res = 0;
   // read mu1, mu2, std1, std2 at the computed patch centers
   assert(s->mu1 && s->mu2 && s->std1 && s->std2);
   for (int c = 0; c < nch; c++) {
      long double v = evaluateIntegralImage(ii, 0, x, y, w, h);
      int pidx1 = _ppla(s->mu1, x, y, c);
      int pidx2 = _pinpla(s->mu2, xx, yy, c);
      // check if pidx2 is inside the valid region
      if (xx < s->mu2->ncol - w + 1 && yy < s->mu2->nrow - h + 1
	  && pidx2 >= 0) {
	 long double s1 = _vv(s->std1, pidx1);
	 long double s2 = _vv(s->std2, pidx2);
	 long double mu1 = _vv(s->mu1, pidx1);
	 long double mu2 = _vv(s->mu2, pidx2);
	 if (s1 > 0 && s2 > 0)
	    res += (v / (w * h) - mu1 * mu2) / (s1 * s2);
	 else
	    res += 0;
      }
   }

   //   assert( res/nch >=-0.1);  // check that nothing has gone horribly wrong
   return 1. - fmin(res / nch, 1.);
}
