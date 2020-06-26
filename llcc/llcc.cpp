/**
 * 
 * Copyright (c) 2020, Jose-Luis Lisani, joseluis.lisani@uib.es
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "library/io_png.h"
#include "library/io_RGB.h"

//Use Marco Mondelli' implementation of MCM filter
//Available at: http://www.ipol.im/pub/art/2013/53/
#include "library/MCMMondelli/FDS_MCM.h"

//Use Pascal Getreuer' implementation of Gaussian filter
//Available at: http://www.ipol.im/pub/art/2013/87/
#include "library/gaussianGetreuer/gaussian_conv.h"

//Use Sylvain Paris and Frédo Durand fast implementation of bilateral filter
//Available at: http://people.csail.mit.edu/sparis/bf/
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include "library/TRUNCATED_KERNEL_BF/include/geom.h"
#include "library/TRUNCATED_KERNEL_BF/include/fast_lbf.h"

//Uncomment this line to generate extra information: weight maps, alpha image
//#define DEBUG_OUTPUT

using namespace std; 
typedef Array_2D<double> image_type;


//Read input image and compute intensity
unsigned char *read_input_image(const char *name,
                                unsigned char **R,
                                unsigned char **G,
                                unsigned char **B,
                                unsigned char **I,
                                int &w, int &h)
{
    size_t ww, hh;
    /* read the PNG input image */
    unsigned char *img = io_png_read_u8_rgb(name, &ww, &hh);
    w=(int) ww;
    h=(int) hh;
    input2RGB(img, R, G, B, w*h);
    RGBtoI(*R, *G, *B, I, w*h); //compute average of RGB values
    
    return img;
}

//Linear stretch from [min(A), max(A)] to [0, 255]
//Implements line 2 of Algorithm 1
void stretch_range(unsigned char *A, int w, int h)
{
    float min=(float) A[0];
    float max=(float) A[0];
    for (int n=0; n < w*h; n++) {
        min=(A[n] < min)?((float) A[n]):(min);
        max=(A[n] > max)?((float) A[n]):(max);
    }
    if (max == min) return;

    for (int n=0; n < w*h; n++) {
        A[n]=255.0f*((float) A[n] - min)/(max-min);
    }

}



//Contrast enhancement: implements lines 2 and 10 of Algorithm 1
//option == 1 -> Gaussian filter
//option == 2 -> MCM filter
//option == 3 -> Bilateral filter
void loglocal_correction(unsigned char *A, int w, int h, 
                         float sigmaI, float sigmaS, float gradth,
                         float gammalog, unsigned char option)
{
    stretch_range(A, w, h);

    float *Anorm=new float[w*h];

    //auxiliary variables to create weight map
    float *mask=new float[w*h];
    
    //auxiliary variables to create weight map image
    #ifdef DEBUG_OUTPUT
    unsigned char *Imask=new unsigned char[w*h];
    #endif
    
    //normalize range to [0, 1]
    for (int n=0; n < w*h; n++) Anorm[n]=A[n]/255.0f;
    
    if (option == 1) {
        Gaussian2D(Anorm, mask, w, h, sigmaS); 
        //creation of weight map image
        #ifdef DEBUG_OUTPUT
        for (int n=0; n < w*h; n++) Imask[n]=(int) (255*mask[n]+0.5f);
        io_png_write_u8("maskG.png", Imask, w, h, 1); //1 channel
        #endif
    } else {
        if (option == 2) {
            mcm_main(Anorm, mask, w, h, sigmaS, gradth/255.0f);
            //creation of weight map image
            #ifdef DEBUG_OUTPUT
            for (int n=0; n < w*h; n++) Imask[n]=(int) (255*mask[n]+0.5f);
            io_png_write_u8("maskM.png", Imask, w, h, 1); //1 channel
            #endif
        } else { //option=3
        
            image_type image(w, h);
            for (int n=0; n < w*h; n++) image(n%w, n/w)=Anorm[n];

            image_type filtered_image(w, h);
            Image_filter::fast_LBF(image, image, sigmaS, sigmaI/255.0f,
                                   false, &filtered_image,&filtered_image);

            for (int n=0; n < w*h; n++) mask[n]=filtered_image(n%w, n/w);
            
            //creation of weight map image
             #ifdef DEBUG_OUTPUT
            for (int n=0; n < w*h; n++) Imask[n]=(int) (255*mask[n]+0.5f);
            io_png_write_u8("maskB.png", Imask, w, h, 1); //1 channel
            #endif
                       
        }
    }
    
   
    //get alpha paremeters values for each pixel
    float *alphavalues=new float[w*h];
    for (int n=0; n < w*h; n++) {
        if (mask[n] <= 0.5) alphavalues[n]=0.5-0.5*pow(mask[n]/0.5, gammalog);
        else alphavalues[n]=-(0.5-0.5*pow((1-mask[n])/0.5, gammalog));
    }

    //creation of alpha values image
    #ifdef DEBUG_OUTPUT
    for (int n=0; n < w*h; n++) Imask[n]=(int) (255*alphavalues[n]+127.5f);
    io_png_write_u8("alphavalues.png", Imask, w, h, 1); //1 channel
    #endif

    //apply tone mappings
    float aag, aagmax;
    int aai;
    for (int n=0; n < w*h; n++) {   
        float alpha=alphavalues[n];
     
        aagmax=255.0f/log(fabs(alpha)*255+1);
        if (alpha > 0) {
            aag=aagmax*log(alpha* (float) A[n] + 1);
        }
        if (alpha < 0) {
            aag=255.0f-aagmax*log(-alpha* (float) (255.0f-A[n]) + 1);
        }
        if (alpha == 0) aag=A[n]; //identity
        
        aai=(int) (aag + 0.5f);
        //clip values outside [0, 255] range
        if (aai < 0) aai=0;
        A[n]=(aai > 255)?(255):(aai);
    }
    
    //save result for 1-channel image
    #ifdef DEBUG_OUTPUT
    for (int n=0; n < w*h; n++) Imask[n]=A[n];
    io_png_write_u8("gray.png", Imask, w, h, 1); //1 channel
    #endif

    //cleanup
    delete[] Anorm;
    delete[] alphavalues;
    delete[] mask;
    
    //delete auxiliary variable used to create weight map image
    #ifdef DEBUG_OUTPUT
    delete[] Imask;
    #endif
}


//Color processing: implements lines 11 and 12 of Algorithm 1
void color_processing(unsigned char *R, unsigned char *G, unsigned char *B,
                      unsigned char *I, unsigned char *Iout, int w, int h)
{
    float factorI;
    for (int n=0; n < w*h; n++) {
        if (I[n] != 0) {
            factorI= (float) Iout[n]/ (float) I[n];
            float rr= (float) R[n] * factorI;
            float gg= (float) G[n] * factorI;
            float bb= (float) B[n] * factorI;
            //check if value is out of range (assume 8-bits images)
            float outmax=(rr > gg)?((rr > bb)?(rr):(bb)):((gg > bb)?(gg):(bb));
            if (outmax > 255) { //out of range: rescale but keeping the R/G/B ratio
                rr*=255.0f/outmax;
                gg*=255.0f/outmax;
                bb*=255.0f/outmax;
            }
            R[n]=(int) (rr+0.5f);
            G[n]=(int) (gg+0.5f);
            B[n]=(int) (bb+0.5f);
        } else {
            R[n]=0;
            G[n]=0;
            B[n]=0;
        }
    }
    
}


//Main function: implements Algorithm 1
int main(int argc, const char **argv)
{

    if (argc < 3) {
        printf("Usage: llcc input.png output.png [option=1] [sigmaS=5] [sigmaI=70/gradth=10]\n");
        printf("       option == 1 --> Use Gaussian filter\n");
        printf("       option == 2 --> Use MCM filter\n");
        printf("       option == 3 --> Use Bilinear filter\n");
        printf("       Parameter of Gaussian filter: sigmaS (spatial scale)\n");
        printf("       Parameter of MCM filter: sigmaS (spatial scale)\n");
        printf("                              : gradth (gradient threshold)\n");
        printf("       Parameter of Bilinear filter: sigmaI (intensity scale)\n");
        printf("                                     sigmaS (spatial scale)\n");
        return EXIT_FAILURE;
    }
    
    const char *namein=argv[1];
    const char *nameout=argv[2];
               
    //Parameters
    int option=1; //default value
    if (argc > 3) option=atoi(argv[3]);
    if ((option != 1) && (option != 2) && (option != 3)) option=1; //use default value if invalid entry
    
    float sigmaS; //scale of Gaussian and mcm filters (if option == 1  or  option == 2)
                  //or spatial scale of bilinear filter (if option == 3)
    float sigmaI;//intensity scale of bilinear filter (only used if option == 3)
    float gradth;//gradient threshold of MCM filter (only used if option == 2)

    //default values
    sigmaS=5.0f;
    sigmaI=70.0f;
    gradth=10.0f;
    if (argc > 4) sigmaS=atof(argv[4]);
    if (argc > 5) {
        if (option == 2) gradth=atof(argv[5]); //5th parameter is gradient threshold
        if (option == 3) sigmaI=atof(argv[5]); //5th parameter is intensity scale
    }
    
    float gammalog=0.05f; //parameter gamma in Equation (2)
               
    //Read input and compute Intensity
    int w, h; //input image dimensions
    unsigned char *img; //image data
    unsigned char *R, *G, *B, *I; //image data organized in R, G and B channels
    img=read_input_image(namein, &R, &G, &B, &I, w, h);
       
    //apply algorithm to intensity channel
    unsigned char *Iout=new unsigned char[w*h];
    memcpy(Iout, I, w*h);
    loglocal_correction(Iout, w, h, sigmaI, sigmaS, gradth, gammalog, option);

    //recover original chrominance
    color_processing(R, G, B, I, Iout, w, h);
    
    //save result
    RGB2output(R, G, B, img, w*h);
    io_png_write_u8(nameout, img, w, h, 3); //3 channels
    
    //cleanup
    free(img);
    delete[] R;
    delete[] G;
    delete[] B;
    delete[] I;
    delete[] Iout;
  
    return EXIT_SUCCESS;

}
               
               
               
