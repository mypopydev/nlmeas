#include <cmath>

void rot1(const double x1, const double y1, const double theta,
          double* x2, double* y2) {
  const double s = sin(theta);
  const double c = cos(theta);
  *x2 = x1 * c - y1 * s;
  *y2 = y1 * c + x1 * s;
}

void rot2(const double x1, const double y1, const double theta,
          double* x2, double* y2) {
  const double s = sin(theta);
  const double c = cos(theta);
  const double z = x1 * (s + c);
  *x2 = z - (x1 + y1) * s;
  *y2 = z + (y1 - x1) * c;
}

void rot3(const double x1, const double y1,
          double* x2, double* y2) {
  const double s = sin(1. * M_PI / 16.);
  const double c = cos(1. * M_PI / 16.);
  const double z = c * (x1 + y1);
  *x2 = (-s-c) * y1 + z;
  *y2 = ( s-c) * x1 + z;
}

// $ g++ -v
// gcc version 4.9.0
// $ g++ -S -m64 -O3 -march=core-avx-i -masm=intel rot.cc
// $ cat rot.s | awk '{print $1}' | grep 'v....d' | sort | uniq -c
//
// rot1:
//       1 vaddsd
//       8 vmovsd
//       4 vmulsd
//       1 vsubsd
// 
// rot2:
//       3 vaddsd
//       8 vmovsd
//       3 vmulsd
//       2 vsubsd
// 
// rot3:
//       3 vaddsd
//       2 vmovsd
//       3 vmulsd

// vim:set ts=2 sw=2 et:
