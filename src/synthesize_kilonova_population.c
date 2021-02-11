#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "helpers.h"
#include "physics.h"
#include "units.h"

int main(int argc, char ** argv){

  /* Results*/
  double d, tv, g0, r0, i0, z0, r1, r10;

  /* Distributions */
  int n_pop;
  int n_run;
  double d_max;

  /* Initialize random number generator */
  srand48(time(0));

  /* read input */
  d_max = atoi(argv[1]) * Mpc;
  n_pop = atoi(argv[2]);

  for (n_run = 0; n_run < n_pop; n_run++){
    d = d_max * pow(rd(), 1./3);
    tv = acos(rd());

    g0 = kngmagwall(tv, d);
    r0 = knrmagwall(tv, d);
    i0 = knimagwall(tv, d);
    z0 = knzmagwall(tv, d);
    r1 = peak_flux_kn(1., d, tv);
    r10 = peak_flux_kn(10., d, tv);
    fprintf(stdout, "%.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
        d / Mpc, tv, g0, 2 * rd() - 1, r0, 2 * rd() - 1, i0, 2 * rd() - 1, z0, 2 * rd() - 1, r1 / microJy, r10 / microJy);
  }
}
