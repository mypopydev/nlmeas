#ifndef FIMAGE_H
#define FIMAGE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


typedef struct fimage {
   int sz[3];			//[col,row,ch]
   int ncol;
   int nrow;
   int nch;
   int N;
   int planar;
   float *gray;			// must mantain both 
   float *v;			// same as gray
//   float v[];    // flexible array is another possibility
} Fimage;

#define FIMAGE_PLANAR 1
#define FIMAGE_VEC 0

// multi channel 2d image
Fimage *malloc_fimage4(int nc, int nr, int nch, int planar);

// 1d vector
Fimage *malloc_fimage1(int numel);

// single channel 2d image
Fimage *malloc_fimage2(int nc, int nr);

// multi channel 2d image planar   rrrr ggggg bbbbb
Fimage *malloc_fimage3pla(int nc, int nr, int nch);

// multi channel 2d image vectorial   rgb rgb rgb rgb
Fimage *malloc_fimage3vec(int nc, int nr, int nch);

// this is very strange, the structure is allocate only 
Fimage *malloc_fimageWrap(float *v, int nc, int nr, int nch, int planar);




// in place transform from vec to planar form of the array
void fimage_vec_to_planar(Fimage * u);

// in place transform from planar to vec form of the array
void fimage_planar_to_vec(Fimage * u);

void fimage_to_vec(Fimage * u);

void fimage_to_planar(Fimage * u);




#endif
