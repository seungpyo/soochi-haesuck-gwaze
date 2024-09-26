#include <math.h>
#include "solve.h"

void muller(
    F f,
    float *x0,
    float *x1,
    float *x2,
    unsigned int max_iter,
    float xacc)
{
  for (int i = 0; i < max_iter; ++i)
  {
    if (fabs(*x2 - *x1) < xacc)
    {
      break;
    }
    float p0 = *x0, p1 = *x1, p2 = *x2, p3;
    float a, b, c;
    c = f(p2);
    float d02 = p0 - p2;
    float d12 = p1 - p2;
    float d20 = p2 - p0;
    float d01 = p0 - p1;
    float b1 = d02 * d02 * (f(p1) - f(p2));
    float b2 = d12 * d12 * (f(p0) - f(p2));
    b = (b1 - b2) / (d02 * d12 * d01);
    float a1 = d12 * (f(p0) - f(p2));
    float a2 = d02 * (f(p1) - f(p2));
    a = (a1 - a2) / (d02 * d12 * d01);
    float D = b * b - 4 * a * c;
    if (D < 0)
    {
      break;
    }
    p3 = p2 - 2 * c / (b + (b >= 0 ? 1 : -1) * sqrt(D));
    *x0 = p1;
    *x1 = p2;
    *x2 = p3;
  }
}

float muller_solver(
    float (*func)(float),
    float x1,
    float x2,
    float xacc)
{
  float p0 = x1, p1 = (x1 + x2) / 2, p2 = x2;
  muller(func, &p0, &p1, &p2, 30, xacc);
  return p2;
}