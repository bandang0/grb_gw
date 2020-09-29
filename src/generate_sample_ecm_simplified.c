#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <libconfig.h>
#include "helpers.h"
#include "physics.h"
#include "units.h"

int main(int argc, char ** argv){

  /* Jet configuration */
  double eiso, gamma_0, n, p, epsilon_e, epsilon_b, tv, tj, d, nu_obs;

  /* Results*/
  double fp, tp, dtvlbi, ppm, kng, knr, agpkn, fpkn, tpkn;

  /* Distributions */
  int n_pop;
  int n_run;
  double ran;
  double p1;
  double p2;
  double e0_min = 1.0e50;
  double e0_star = 2.0e52;
  double e0_max = 1.0e53;
  double d_max = 530.0e0 * Mpc;
  double fstar;

  /* Jet constants */
  srand(time(0));
  epsilon_e = 0.1;
  p = 2.2;
  tj = 0.1;
  nu_obs = 3.0e9;

  /* thresholds */
  double fvla = 15 * microJy;
  double ainc = 0.8;
  double adec = -p;

  /* read input */
  if (strcmp(argv[1], "g16") == 0){
    p1 = 0.53;
    p2 = 3.4;
  }
  else if (strcmp(argv[1], "wp15") == 0){
    p1 = 1.9;
    p2 = 3.;
  }
  else{
    fprintf(stderr, "Wrong energy function label.\n");
    exit(1);
  }
  n_pop = atoi(argv[2]);
  fstar = 1. / (1. + pow(e0_star, p2 - p1) * (1. - p1)
  * (pow(e0_max, 1. - p2) - pow(e0_star, 1. - p2))
  / ((1. - p2) * (pow(e0_star, 1. - p1) - pow(e0_min, 1. - p1))));

  for (n_run = 0; n_run < n_pop; n_run++){
    /* Generate random jet paradeceters */
    d = d_max * rd();
    n = pow(10., randn(-3, 0.75));
    epsilon_b = pow(10., randn(-3, 0.75));
    tv = acos(rd());
    ran = rd();
    if (ran < fstar){
       eiso = e0_min
      * pow(1 - (ran/fstar)*(1-pow(e0_star/e0_min,1 - p1)), 1/(1-p1));
    }
    else {
       eiso = e0_star
      * pow(1 - ((ran-fstar)/(1-fstar))*(1-pow(e0_max/e0_star, 1-p2)), 1/(1-p2));
    }

    fp = pl_peak_flux(d, n, epsilon_e, epsilon_b, eiso, tj, tv);
    tp = pl_peak_time_ex(d, n, epsilon_e, epsilon_b, eiso, tj, tv);
    if((tv > tj) && (tv < 1.)){
      ppm = proper_motion(1/tv, tv, d);
    }else{
      ppm = 0.;
    }
    if (fp > fvla){
      dtvlbi = tp * (pow(fvla/fp, 1/adec)-pow(fvla/fp, 1/ainc));
    }
    else{
      dtvlbi = 0.;
    }
    //tpkn = knagtp(n);
    //fpkn = knagfp(epsilon_b, n, d);
    kng = kngmagwall(tv, d);
    knr = knrmagwall(tv, d);
    //agpkn = fp * pow(tpkn/tp, (tpkn > tp ? adec : ainc));
    fprintf(stdout, "%e %e %e %e %e %e %e %e %e %e %e\n",
        d, n, eiso, epsilon_b,
        tv, tp, fp, ppm,
        dtvlbi,
        kng, knr);
  }
}
