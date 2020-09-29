#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <libconfig.h>
#include "helpers.h"
#include "physics.h"
#include "units.h"

int main(int argc, char ** argv){

  /* Jet configuration */
  double eiso, n, epsilon_e, epsilon_b, tv, tj, d;

  /* Physical data */
  double n_run;

  /* Distributions */
  int n_pop = 10000000;
  double ran, ran2;
  // G16
  double p1 = 0.53;
  double p2 = 3.4;
  // WP15
  //double p1 = 1.9;
  //double p2 = 3.;
  double e0_min;
  double e0_star = 2.0e52;
  double e0_max = 1.0e53;
  //double d_max = 143.0e0;
  double d_max;
  double fstar = 1. / (1. + pow(e0_star, p2 - p1) * (1. - p1)
  * (pow(e0_max, 1. - p2) - pow(e0_star, 1. - p2))
  / ((1. - p2) * (pow(e0_star, 1. - p1) - pow(e0_min, 1. - p1))));

  /* Jet constants */
  srand(time(0));
  epsilon_e = 0.1;
  tj = atof(argv[1]);
  d_max = atof(argv[2]);
  e0_min = pow(10.0, atof(argv[3]));

  for (n_run = 1; n_run < n_pop; n_run++){
    /* Generate random jet parameters */
    n = pow(10., -3.);
    epsilon_b = pow(10., randn(-3, 0.75));
    tv = acos(rd());
    d = d_max * pow( rd(), 1. / 3);
    ran = rd();
    if (ran < fstar){
       eiso = e0_min
      * pow(1 - (ran/fstar)*(1-pow(e0_star/e0_min,1 - p1)), 1/(1-p1));
    }
    else {
       eiso = e0_star
      * pow(1 - ((ran-fstar)/(1-fstar))*(1-pow(e0_max/e0_star, 1-p2)), 1/(1-p2));
    }

    /* print away */
    fprintf(stdout, "%e %e %e %e %e %e %e %e %e %e\n",
        d, n, eiso, epsilon_b,
        tv, tj, pl_peak_time_ex(d, n, epsilon_e, epsilon_b, eiso, tj, tv),
                pl_peak_time_no_ex(d, n, epsilon_e, epsilon_b, eiso, tj, tv),
                pl_peak_flux(d, n, epsilon_e, epsilon_b, eiso, tj, tv),
                pl_mu(d, n, epsilon_e, epsilon_b, eiso, tj, tv));
  }
}
