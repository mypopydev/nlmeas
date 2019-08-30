// https://unix4lyfe.org/dct-1d/
#include <cmath>
#include <cstdio>

const constexpr double c1  = cos(M_PI *  1. / 16.) * sqrt(2. / 8.);
const constexpr double c31 = cos(M_PI * 31. / 16.) * sqrt(2. / 8.);

int main() {
  printf("%a\t%.54f\n", c1, c1);
  printf("%a\t%.54f\n", c31, c31);
}

/* 0x1.f6297cff75cbp-2  0.490392640201615215289621119154617190361022949218750000
 * 0x1.f6297cff75cafp-2 0.490392640201615159778469887896790169179439544677734375
 */
