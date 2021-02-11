#ifndef HELPERS_H
#define HELPERS_H

#include <stdlib.h>
#include <math.h>

double rd(){
  return drand48();
}

double randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
  {
    call = !call;
    return (mu + sigma * X2);
  }
  do
  {
    U1 = -1 + rd() * 2;
    U2 = -1 + rd() * 2;
    W = pow (U1, 2) + pow (U2, 2);
  }
  while (W >= 1. || W == 0.);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * X1);
}

#endif
