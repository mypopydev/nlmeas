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

#include <stdlib.h>

/**
 * \brief Convert color data from 1 array with format 
 * RR...RRGG...GGBB...BB to 3 arrays with formats RRR..., GGG..., BBB... 
 *
 *
 * @param[in] input input array
 * @param[out] RR, GG, BB output arrays 
 * @param[in] size number of colors
 *
 */
void input2RGB(unsigned char *input, 
               unsigned char **RR, unsigned char **GG, unsigned char **BB, 
               int size)
{
    unsigned char *R, *G, *B;
    int n;
    
    R=new unsigned char[size];
    G=new unsigned char[size];
    B=new unsigned char[size];
    for (n=0; n < size; n++) {
        R[n]=input[n];
        G[n]=input[size+n];
        B[n]=input[2*size+n];
    }
    
    *RR=R;
    *GG=G;
    *BB=B;
}

/**
 * \brief Convert color data from 3 arrays with formats 
 * RRR..., GGG..., BBB...  to 1 array with format RR...RRGG...GGBB...BB 
 *
 *
 * @param[in] R, G, B input arrays 
 * @param[out] output input array
 * @oaram[in] size number of colors
 *
 */
void RGB2output(unsigned char *R, unsigned char *G, unsigned char *B, 
                unsigned char *output, int size)
{
    int n;
    
    for (n=0; n < size; n++) {
        output[n]=R[n];
        output[size+n]=G[n];
        output[2*size+n]=B[n];
    }
}

/**
 * \brief Compute Intensity channel as average of RGB channels
 *
 */
void RGBtoI(unsigned char *R, unsigned char *G, unsigned char *B,
            unsigned char **II, int size)
{
    unsigned char *I=new unsigned char[size];
    for (int n=0; n < size; n++) 
        I[n]=(unsigned char) (0.5f + ((float) R[n] + (float) G[n] + (float) B[n])/3.0f);
    
    *II=I;
}

