/* DCT and IDCT - listing 2
 * Copyright (c) 2001 Emil Mikulic.
 * http://unix4lyfe.org/dct/
 *
 * Feel free to do whatever you like with this code.
 * Feel free to credit me.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "targa.h"



typedef uint8_t byte;

#define pixel(i,x,y) ( (i)->image_data[((y)*( (i)->width ))+(x)] )

#define DONTFAIL(x) do { tga_result res;  if ((res=x) != TGA_NOERR) { \
	printf("Targa error: %s\n", tga_error(res)); \
	exit(EXIT_FAILURE); } } while(0)



void load_tga(tga_image *tga, const char *fn)
{
	DONTFAIL( tga_read(tga, fn) );

	printf("Loaded %dx%dx%dbpp targa (\"%s\").\n",
		tga->width, tga->height, tga->pixel_depth, fn);

	if (!tga_is_mono(tga)) DONTFAIL( tga_desaturate_rec_601_1(tga) );
	if (!tga_is_top_to_bottom(tga)) DONTFAIL( tga_flip_vert(tga) );
	if (tga_is_right_to_left(tga)) DONTFAIL( tga_flip_horiz(tga) );

	if ((tga->width % 8 != 0) || (tga->height % 8 != 0))
	{
		printf("Width and height must be multiples of 8\n");
		exit(EXIT_FAILURE);
	}
}



#ifndef PI
 #ifdef M_PI
  #define PI M_PI
 #else
  #define PI 3.14159265358979
 #endif
#endif



/* DCT rows and columns separately
 *
 * C(i) = a(i)/2 * sum for x=0 to N-1 of
 *   s(x) * cos( pi * i * (2x + 1) / 2N )
 *
 * if x = 0, a(x) = 1/sqrt(2)
 *      else a(x) = 1
 */

void dct_1d(double *in, double *out, const int count)
{
	int x, u;

	for (u=0; u<count; u++)
	{
		double z = 0;

		for (x=0; x<count; x++)
		{
			z += in[x] * cos(PI * (double)u * (double)(2*x+1)
				/ (double)(2*count));
		}

		if (u == 0) z *= 1.0 / sqrt(2.0);
		out[u] = z/2.0;
	}
}

void dct(const tga_image *tga, double data[8][8],
	const int xpos, const int ypos)
{
	int i,j;
	double in[8], out[8], rows[8][8];

	/* transform rows */
	for (j=0; j<8; j++)
	{
		for (i=0; i<8; i++)
			in[i] = (double) pixel(tga, xpos+i, ypos+j);
		dct_1d(in, out, 8);
		for (i=0; i<8; i++) rows[j][i] = out[i];
	}

	/* transform columns */
	for (j=0; j<8; j++)
	{
		for (i=0; i<8; i++)
			in[i] = rows[i][j];
		dct_1d(in, out, 8);
		for (i=0; i<8; i++) data[i][j] = out[i];
	}
}



/* play with this bit */
void quantize(double dct_buf[8][8])
{
	int x,y;

	for (y=0; y<8; y++)
	for (x=0; x<8; x++)
	if (x > 3 || y > 3) dct_buf[y][x] = 0.0;
}



#define COEFFS(Cu,Cv,u,v) { \
	if (u == 0) Cu = 1.0 / sqrt(2.0); else Cu = 1.0; \
	if (v == 0) Cv = 1.0 / sqrt(2.0); else Cv = 1.0; \
	}

void idct(tga_image *tga, double data[8][8], const int xpos, const int ypos)
{
	int u,v,x,y;

	/* iDCT */
	for (y=0; y<8; y++)
	for (x=0; x<8; x++)
	{
		double z = 0.0;

		for (v=0; v<8; v++)
		for (u=0; u<8; u++)
		{
			double S, q;
			double Cu, Cv;
			
			COEFFS(Cu,Cv,u,v);
			S = data[v][u];

			q = Cu * Cv * S *
				cos((double)(2*x+1) * (double)u * PI/16.0) *
				cos((double)(2*y+1) * (double)v * PI/16.0);

			z += q;
		}

		z /= 4.0;
		if (z > 255.0) z = 255.0;
		if (z < 0) z = 0.0;

		pixel(tga, x+xpos, y+ypos) = (uint8_t) z;
	}
}



int main()
{
	tga_image tga;
	double dct_buf[8][8];
	int i, j, k, l;

	load_tga(&tga, "in.tga");

	k = 0;
	l = (tga.height / 8) * (tga.width / 8);
	for (j=0; j<tga.height/8; j++)
	for (i=0; i<tga.width/8; i++)
	{
		dct(&tga, dct_buf, i*8, j*8);
		quantize(dct_buf);
		idct(&tga, dct_buf, i*8, j*8);
		printf("processed %d/%d blocks.\r", ++k,l);
		fflush(stdout);
	}
	printf("\n");

	DONTFAIL( tga_write_mono("out.tga", tga.image_data,
		tga.width, tga.height) );

	tga_free_buffers(&tga);
	return EXIT_SUCCESS;
}
