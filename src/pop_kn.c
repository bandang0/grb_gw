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
  double d, tv, g, r, i, z, r1, r10;

  /* Distributions */
  int n_pop;
  int n_run;
  double d_max;

  /* Jet constants */
  srand(time(0));

  /* read input */
  d_max = atoi(argv[1]) * Mpc;
  n_pop = atoi(argv[2]);

  for (n_run = 0; n_run < n_pop; n_run++){
    d = d_max * pow(rd(), 1./3);
    tv = acos(rd());

    g = kngmagwall(tv, d);
    r = knrmagwall(tv, d);
    i = knimagwall(tv, d);
    z = knzmagwall(tv, d);
    r1 = peak_flux_kn(1., d, tv);
    r10 = peak_flux_kn(1., d, tv);
    fprintf(stdout, "%e %e %e %e %e %e %e %e\n",
        d / Mpc, tv, g, r, i, z, r1 / microJy, r10 / microJy);
  }
}
