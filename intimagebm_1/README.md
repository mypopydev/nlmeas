Integral Images for Block Matching
==================================
Subpixel stereo matching and 2D block matching using integral images.

Gabriele Facciolo, Nicolas Limare and Enric Meinhardt
<{facciolo,nicolas.limare,enric.meinhardt}@ens-cachan.fr>,
CMLA, ENS Cachan, France

Complete IPOL article available at: <http://www.ipol.im/pub/art/2013/57/> 

Current Version: 20170813


Build instructions
==================
To compile this code cmake, make and a C99 compiler (ie. gcc) are required.
The FFTW3 library must be installed in the system; libpng,libjpeg and libtiff
are optional but WE RECOMMEND THEM.

Compile the binaries by running the following line in the current directory: 

    $ mkdir build; cd build; cmake ..; make

the process will produce two binaryes: simplebm and simplestereo



Usage
=====
Running the binaries without parameters will print usage instructions as:

    $ build/simplestereo

    Block matching stereo using integral image
    usage: ./simplestereo [-t method] [-w window_width] [-h window_height] 
         [-r disparity_min] [-R disparity_max] leftIM rightIM disparity mincost 
         [disparity_Right mincost_Right]
       The SUBPIXEL environ. var. determines the subpixel precision (default 1)
       The SLICED environ. var. determines the width of the bands to 
       divide the computation ( default max(32,h*2) ).
       method can be: SSD sum of squaerd differences
                   SAD sum of absolute differences
                   ZSSD zero mean SSD
                   SSDNorm L^2-normalized SSD
                   NCC normalized cross correlation
                   AFF “affine” similarity measure
                   LIN linearized “affine” similarity measure
                   BTSSD SSD with Birchfield-Tomasi cost
                   BTSAD SAD with Birchfield-Tomasi cost
                   BT2DSSD SSD with Birchfield-Tomasi 2D cost
                   BT2DSAD SAD with Birchfield-Tomasi 2D cost



## Example calls:
Call the stereo algorithm and the block matching as in these examples: 

    $ SUBPIXEL=4 build/simplestereo -t ZSSD -r -45 -R 0 \
      testdata/{a,b}.png /tmp/st{dL,cL,dR,cR}.tif

    $ SUBPIXEL=1 build/simplebm -t ZSSD -r -5 -R 5 -s -5 -S 5 \
      testdata/{a,b}.png /tmp/bm{dL,cL,dR,cR}.tif

The output files *{dL,dR}.tif correspond to the disparity/flow fileds, 
while *{cL,cR}.tif are the corresponding matching costs. The extension 
of the output determines the format of the output file, we recommend tiff
as it permits to store floating point values. 
[vflip](https://github.com/gfacciol/omniflip) can be used to display real valued tiffs. 




Files and main functions 
========================
The blockmatch*.c files contain the implementation of four functions
that perform block matching with integral images.

    stereoInt               in blockmatch.c
    stereoIntSlice          in blockmatch_stereo.c wraps stereoInt
    stereoIntSliceSubpix    in blockmatch_stereo.c wraps stereoIntSlice and stereoInt
    bmIntSlice              in blockmatch_2d.c wraps stereoIntSliceSubpix
    bmIntSliceSubpix        in blockmatch_2d.c wraps bmIntSlice 

* stereoInt: All the functions end up calling stereoInt which implements block
   matching with integer disparities along the horizontal direction.
* stereoIntSlice: Slices the input image in horizontal bands and for each band
   calls one instance of stereoInt in parallel.
   This function uses openmp if available
* stereoIntSliceSubpix: Computes horizontal disparities, by calling stereoIntSlice
   several times with images shifted by a subpixel quantity <1, and combining
   the results of the different runs.
* bmIntSlice: Computes 2D block matching with integer vertical offsets by shifting the input
   images vertically and by calling stereoIntSliceSubpix for each shift.
   Then combines the results to obtain the 2d offset map.
   This function is just for internal use, better use: bmIntSliceSubpix
* bmIntSliceSubpix: Computes 2D block matching with subpixel offsets by shifting the input
   images vertically and by calling bmIntSlice for each shift.
   Then combines the results to obtain the 2d offset map.

