// https://unix4lyfe.org/dct-1d/
#include <cmath>

void dct_ii_8(const double x[8], double X[8]) {
  for (int k = 0; k < 8; ++k) {
    double sum = 0.;
    double s = (k == 0) ? sqrt(.5) : 1.;
    for (int n = 0; n < 8; ++n) {
      sum += s * x[n] * cos(M_PI * (n + .5) * k / 8);
    }
    X[k] = sum * sqrt(2. / 8);
  }
}
