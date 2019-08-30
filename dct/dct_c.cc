// https://unix4lyfe.org/dct-1d/
#include <cmath>

static const double c0 = 1. / sqrt(2.) * sqrt(2. / 8.);
static const double c1 = cos(M_PI * 1. / 16.) * sqrt(2. / 8.);
static const double c2 = cos(M_PI * 2. / 16.) * sqrt(2. / 8.);
static const double c3 = cos(M_PI * 3. / 16.) * sqrt(2. / 8.);
static const double c4 = cos(M_PI * 4. / 16.) * sqrt(2. / 8.);
static const double c5 = cos(M_PI * 5. / 16.) * sqrt(2. / 8.);
static const double c6 = cos(M_PI * 6. / 16.) * sqrt(2. / 8.);
static const double c7 = cos(M_PI * 7. / 16.) * sqrt(2. / 8.);

#define a x[0]
#define b x[1]
#define c x[2]
#define d x[3]
#define e x[4]
#define f x[5]
#define g x[6]
#define h x[7]

// Matrix multiplication in full.
void dct_ii_8a(const double x[8], double X[8]) {
  X[0] = a*c0 + b*c0 + c*c0 + d*c0 + e*c0 + f*c0 + g*c0 + h*c0;
  X[1] = a*c1 + b*c3 + c*c5 + d*c7 - e*c7 - f*c5 - g*c3 - h*c1;
  X[2] = a*c2 + b*c6 - c*c6 - d*c2 - e*c2 - f*c6 + g*c6 + h*c2;
  X[3] = a*c3 - b*c7 - c*c1 - d*c5 + e*c5 + f*c1 + g*c7 - h*c3;
  X[4] = a*c4 - b*c4 - c*c4 + d*c4 + e*c4 - f*c4 - g*c4 + h*c4;
  X[5] = a*c5 - b*c1 + c*c7 + d*c3 - e*c3 - f*c7 + g*c1 - h*c5;
  X[6] = a*c6 - b*c2 + c*c2 - d*c6 - e*c6 + f*c2 - g*c2 + h*c6;
  X[7] = a*c7 - b*c5 + c*c3 - d*c1 + e*c1 - f*c3 + g*c5 - h*c7;
}

// Factoring.
void dct_ii_8b(const double x[8], double X[8]) {
  double ah = a - h;
  double bg = b - g;
  double cf = c - f;
  double de = d - e;
  double adeh = a - d - e + h;
  double bcfg = b - c - f + g;
  X[0] = (a + b + c + d + e + f + g + h)*c0;
  X[1] = ah*c1 + bg*c3 + cf*c5 + de*c7;
  X[2] = adeh*c2 + bcfg*c6;
  X[3] = ah*c3 - bg*c7 - cf*c1 - de*c5;
  X[4] = (a - b - c + d + e - f - g + h)*c4;
  X[5] = ah*c5 - bg*c1 + cf*c7 + de*c3;
  X[6] = adeh*c6 - bcfg*c2;
  X[7] = ah*c7 - bg*c5 + cf*c3 - de*c1;
}

// Trading multiplies for adds in "rotation" operations:
//   y0 =  a * x0 + b * x1     =  (b - a) * x1 + a * (x0 + x1)                          
//   y1 = -b * x0 + a * x1     = -(a + b) * x0 + a * (x0 + x1)
void dct_ii_8c(const double x[8], double X[8]) {
  const double adehp = a + d + e + h;
  const double bcfgp = b + c + f + g;

  X[0] = (adehp + bcfgp)*c0;
  X[4] = (adehp - bcfgp)*c4;

  const double adeh = a - d - e + h;
  const double bcfg = b - c - f + g;
  const double adeh26 = adeh * (c2 + c6);
  X[2] = (bcfg - adeh)*c6 + adeh26;
  X[6] = adeh26 - (bcfg + adeh)*c2;

  const double ah = a - h;
  const double bg = b - g;
  const double cf = c - f;
  const double de = d - e;
  const double cf53 = cf * (c5 + c3);
  const double ah17 = ah * (c1 + c7);
  const double ah35 = ah * (c3 + c5);
  const double bg71 = bg * (c7 + c1);
  const double deah = de - ah;
  //const double bgmcf = bg - cf;  // Not worth it.
  X[1] = deah*c7 + ah17 + (bg - cf)*c3 + cf53;
  X[5] = deah*c3 + ah35 + (bg + cf)*c7 - bg71;

  //const double ahde = ah + de;  // Not worth it.
  X[3] = ah35 - (ah + de)*c5 + (bg - cf)*c1 - bg71;
  X[7] = cf53 - (bg + cf)*c5 + ah17 - (ah + de)*c1;
}

#undef a
#undef b
#undef c
#undef d
#undef e
#undef f
#undef g
#undef h
