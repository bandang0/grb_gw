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
  double eiso, gamma_0, n, p, epsilon_e, epsilon_b, tv, tj, d, nu_obs;
  double r_dec, r_nr, miso, omega, dr0, inc_r;

  /* Physical data */
  double r, dr, last_r;
  int n_step, n_run;
  double t_obs, pt, tinvlbi, toutvlbi;
  double t;
  double texz;
  double Niso, nz;
  double b;
  double rho_s;
  double eps_s;
  double beta, ppm;
  double gam, pg;
  double F_nu_obs, max_f_nu_obs, pf;
  double tpkn, fpkn, agpkn, kng, knr;
  double Pmax;
  double nu_a, nu_i, nu_c;
  double gam_i;
  double gam_c;
  double dopp, pdopp;

  /* Distributions */
  int n_pop = 2000000;
  double ran;
  //G16
  double p1 = 0.53;
  double p2 = 3.4;
  // WP15
  //double p1 = 1.9;
  //double p2 = 3.;
  double e0_min = 1.0e50;
  double e0_star = 2.0e52;
  double e0_max = 1.0e53;
  double d_max = 500.0e0;
  double fstar = 1. / (1. + pow(e0_star, p2 - p1) * (1. - p1)
  * (pow(e0_max, 1. - p2) - pow(e0_star, 1. - p2))
  / ((1. - p2) * (pow(e0_star, 1. - p1) - pow(e0_min, 1. - p1))));

  /* Jet constants */
  srand(time(0));
  gamma_0 = 100.0;
  epsilon_e = 0.1;
  p = 2.2;
  tj = 0.1;
  nu_obs = 3.0e9;

  /* thresholds */
  double fpminvlbi = 15 * microJy;

  for (n_run = 1; n_run < n_pop; n_run++){
    /* Generate random jet parameters */
    d = d_max * rd() * Mpc;
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

    /* KNAG */
    tpkn = knagtp(n);
    fpkn = knagfp(epsilon_b, n, d);
    kng = kngmag(tv, d);
    knr = knrmag(tv, d);


    /* Global variables */
    r_dec = pow(3 * eiso
      / (4 * PI * n * mProton * pow(gamma_0, 2) * pow(cLight, 2)), 1. /3 );
    omega = 2 * PI * (1 - cos(tj));

    /* Spatial steps */
    dr0 = r_dec / 2048;
    inc_r = pow(2048, 1. / 4096);

    /* Initial conditions */
    max_f_nu_obs = 0.0;
    texz = 0.0;
    t_obs = 0.0;
    t = 0.0;
    n_step = 1;
    r = 0.0;
    last_r = 0.0;
    tinvlbi = 2.0e7 * day;
    toutvlbi = -1.;

    /* Principal loop */
    while (t_obs / day < 2.0e7){
      /* Dynamics */
      last_r = r;
      r = dr0 * pow(inc_r, n_step - 1);
      dr = r - last_r;
      gam = rel_dec_gamma(r_dec, gamma_0, r);
      beta = sqrt(1. - 1. / pow(gam, 2));
      Niso = 4. * PI * n * pow(r, 3) / 3;
      texz = r / (gam * beta * cLight);
      t = t + dr / (beta * cLight);

      /* Shock frame conditions */
      nz = (3 * gam + 4) * n;
      rho_s = mProton * (3 * gam + 4) * n;
      eps_s = (gam - 1) * pow(cLight, 2);
      b = sqrt(8 * PI * epsilon_b * rho_s * eps_s);

      /* Lorentz factors and frequencies */
      gam_i = (p - 2.) * mProton * epsilon_e * (gam - 1)
      / ((p - 1.) * mElectron);
      gam_c = 6 * PI * mElectron * cLight
      / (SigmaThomson * pow(b, 2) * texz);
      nu_i = eElectron * pow(gam_i, 2) * b
      / (2 * PI * mElectron * cLight);
      nu_c = eElectron * pow(gam_c, 2) * b
      / (2 * PI * mElectron * cLight);
      nu_a = calc_nu_a(texz, nu_i, nu_c, nz, p, gam_c, gam_i);

      Pmax = mElectron * pow(cLight, 2) * SigmaThomson * b
      / (4 * PI * pow(d, 2) * 3 * eElectron);

      if (tv < tj){
        t_obs = t - r / cLight;
        dopp = 1 / (gam * (1 - beta));
        F_nu_obs = Pmax
        * (omega * Niso / (4 * PI)) * pow(dopp, 2)
        * elem_fluxz(nu_i, nu_c, p, nu_obs / dopp)
        * fmin(taz(texz, nu_i, nu_c, nz, p, gam_c, gam_i, nu_obs / dopp) / texz, 1.)
        * fmin(2 * PI * (1 - cos(1/gam))/omega, 1.);
      }
      else{
        t_obs = t - r * cos(tv)/ cLight;
        dopp = 1 / (gam * (1 - cos(tv) * beta));
        F_nu_obs = Pmax
        * (omega * Niso / (4 * PI)) * pow(dopp, 2)
        * elem_fluxz(nu_i, nu_c, p, nu_obs / dopp)
        * fmin(taz(texz, nu_i, nu_c, nz, p, gam_c, gam_i, nu_obs / dopp) / texz, 1.);
      }

      /* update maximum values */
      if (max_f_nu_obs < F_nu_obs){
        max_f_nu_obs = F_nu_obs;
        pt = t_obs;
        pf = F_nu_obs;
        pg = gam;
        pdopp = dopp;
      }

      /* check VLBI detectability */
      if (F_nu_obs > fpminvlbi){
        if(t_obs < tinvlbi){
          tinvlbi = t_obs;
        }
        if (t_obs > toutvlbi){
          toutvlbi = t_obs;
        }
      }

      /* get jet AG value at KNAG peak */
      if (t_obs < tpkn){
        agpkn = F_nu_obs;
      }
      n_step = n_step + 1;
    }
    if (tv > tj){
      ppm = proper_motion(pg, tv, d);
    }
    else{
      ppm = 0.;
    }
    fprintf(stdout, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
        d / Mpc, n, eiso, epsilon_b,
        tv, pt / day, pf / microJy, ppm / (mas/month),
        tinvlbi / day, toutvlbi / day,
        tpkn / day, fpkn / microJy, agpkn / microJy, kng, knr);
  }
}
