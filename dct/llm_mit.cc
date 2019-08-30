#include <cmath>

// This code was adapted from:
// http://web.mit.edu/seven/src/AFNI/dct8.c
//
// In particular, the /sqrt(8) scaling was added.
//
void llm_mit(const double in[8], double out[8]) {
  double x0,x1,x2,x3,x4,x5,x6,x7,x8;

  static const double c1=cos(M_PI/16.);
  static const double s1=sin(M_PI/16.);
  static const double c3=cos(3.*M_PI/16.);
  static const double s3=sin(3.*M_PI/16.);
  static const double r2c6=sqrt(2.)*cos(6.*M_PI/16.);
  static const double r2s6=sqrt(2.)*sin(6.*M_PI/16.);
  static const double r2=sqrt(2.);

  x0 = in[0];
  x1 = in[1];
  x2 = in[2];
  x3 = in[3];
  x4 = in[4];
  x5 = in[5];
  x6 = in[6];
  x7 = in[7];

  /*Stage 1*/
  x8=x7+x0;
  x0-=x7;
  x7=x1+x6;
  x1-=x6;
  x6=x2+x5;
  x2-=x5;
  x5=x3+x4;
  x3-=x4;

  /* Stage 2 */
  x4=x8+x5; x8-=x5;  x5=x7+x6; x7-=x6;
  x6=c1*(x1+x2); x2=(-s1-c1)*x2+x6; x1=(s1-c1)*x1+x6;
  x6=c3*(x0+x3); x3=(-s3-c3)*x3+x6; x0=(s3-c3)*x0+x6;

  /* Stage 3 */
  x6=x4+x5; x4-=x5;
  x5=r2c6*(x7+x8); x7=(-r2s6-r2c6)*x7+x5; x8=(r2s6-r2c6)*x8+x5;
  x5=x0+x2;x0-=x2; x2=x3+x1; x3-=x1;

  /* Stage 4 */
  out[0]=(x6) / sqrt(8);
  out[4]=(x4) / sqrt(8.);
  out[2]=(x8) / sqrt(8.);
  out[6]=(x7) / sqrt(8.);
  out[7]=(x2-x5) / sqrt(8.);
  out[1]=(x2+x5) / sqrt(8.);
  out[3]=(x3*r2) / sqrt(8.);
  out[5]=(x0*r2) / sqrt(8.);
}

// vim:set ts=2 sw=2 et:
