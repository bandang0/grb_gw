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
  double eiso, gamma_0, n, p, epsilon_e, epsilon_b, tv, tj, d, logn;
  double nu_radio, nu_r, nu_x;
  double r_dec, r_nr, miso, omega, dr0, inc_r;

  /* Physical data */
  double r, dr, last_r;
  int n_step, n_run;
  double t_obs, pt;
  double t;
  double texz;
  double Niso, nz;
  double b;
  double rho_s;
  double eps_s;
  double beta;
  double gam, pg;
  double F_nu_obs_radio, max_f_nu_obs_radio;
  double F_nu_obs_r, max_f_nu_obs_r;
  double F_nu_obs_x, max_f_nu_obs_x;
  double Pmax;
  double nu_a, nu_i, nu_c;
  double gam_i;
  double gam_c;
  double dopp;

  /* Distributions */
  int n_pop = 100000;
  double n_gw, n_radio, n_r, n_x, n_m;
  double ran;
  //G16
  //double p1 = 0.53;
  //double p2 = 3.4;
  // WP15
  double p1 = 1.9;
  double p2 = 3.;
  double e0_min = 1.0e50;
  double e0_star = 2.0e52;
  double e0_max = 1.0e53;
  //
  double d_max = 143.0e0;
  double fstar = 1. / (1. + pow(e0_star, p2 - p1) * (1. - p1)
  * (pow(e0_max, 1. - p2) - pow(e0_star, 1. - p2))
  / ((1. - p2) * (pow(e0_star, 1. - p1) - pow(e0_min, 1. - p1))));

  /* Jet constants */
  srand(time(0));
  gamma_0 = 100.0;
  epsilon_e = 0.1;
  p = 2.2;
  epsilon_b = 1.0e-3;
  tj = 0.1;
  nu_radio = 3.0e9;
  nu_r = 4.56e14;
  nu_x = 2.42e17;

  for (logn =2.; logn < 4.5; logn = logn + 0.25){
    n_gw = 0.0;
    n_radio = 0.0;
    n_r = 0.0;
    n_x = 0.0;
    n_m = 0.0;

    for (n_run = 1; n_run < n_pop; n_run++){
      /* Generate random jet parameters */
      n = pow(10.0, logn);
      tv = acos(rd());
      d = d_max * pow(rd(), 1. / 3) * Mpc;

      if (1. + 6. * pow(cos(tv), 2) + pow(cos(tv), 4)
      < 8 * pow(d / Mpc, 2) / pow(d_max, 2)){
          continue;
      }
      else{
            n_gw = n_gw + 1.;
      }
      ran = rd();
      if (ran < fstar){
         eiso = e0_min
        * pow(1 - (ran/fstar)*(1-pow(e0_star/e0_min,1 - p1)), 1/(1-p1));
      }
      else {
         eiso = e0_star
        * pow(1 - ((ran-fstar)/(1-fstar))*(1-pow(e0_max/e0_star, 1-p2)), 1/(1-p2));
      }

      /* Global variables */
      r_dec = pow(3 * eiso
        / (4 * PI * n * mProton * pow(gamma_0, 2) * pow(cLight, 2)), 1. /3 );
      omega = 2 * PI * (1 - cos(tj));

      /* Spatial steps */
      dr0 = r_dec / 2048;
      inc_r = pow(2048, 1. / 4096);

      /* Initial conditions */
      max_f_nu_obs_radio = 0.0;
      max_f_nu_obs_r = 0.0;
      max_f_nu_obs_x = 0.0;
      texz = 0.0;
      t_obs = 0.0;
      t = 0.0;
      n_step = 1;
      r = 0.0;
      last_r = 0.0;

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
          F_nu_obs_radio = Pmax
          * (omega * Niso / (4 * PI)) * pow(dopp, 2)
          * elem_fluxz(nu_i, nu_c, p, nu_radio / dopp)
          * fmin(taz(texz, nu_i, nu_c, nz, p, gam_c, gam_i, nu_radio / dopp) / texz, 1.)
          * fmin(2 * PI * (1 - cos(1/gam))/omega, 1.);

          F_nu_obs_r = Pmax
          * (omega * Niso / (4 * PI)) * pow(dopp, 2)
          * elem_fluxz(nu_i, nu_c, p, nu_r / dopp)
          * fmin(taz(texz, nu_i, nu_c, nz, p, gam_c, gam_i, nu_r / dopp) / texz, 1.)
          * fmin(2 * PI * (1 - cos(1/gam))/omega, 1.);

          F_nu_obs_x = Pmax
          * (omega * Niso / (4 * PI)) * pow(dopp, 2)
          * elem_fluxz(nu_i, nu_c, p, nu_x / dopp)
          * fmin(taz(texz, nu_i, nu_c, nz, p, gam_c, gam_i, nu_x / dopp) / texz, 1.)
          * fmin(2 * PI * (1 - cos(1/gam))/omega, 1.);
        }
        else{
          t_obs = t - r * cos(tv)/ cLight;
          dopp = 1 / (gam * (1 - cos(tv) * beta));
          F_nu_obs_radio = Pmax
          * (omega * Niso / (4 * PI)) * pow(dopp, 2)
          * elem_fluxz(nu_i, nu_c, p, nu_radio / dopp)
          * fmin(taz(texz, nu_i, nu_c, nz, p, gam_c, gam_i, nu_radio / dopp) / texz, 1.);

          F_nu_obs_r = Pmax
          * (omega * Niso / (4 * PI)) * pow(dopp, 2)
          * elem_fluxz(nu_i, nu_c, p, nu_r / dopp)
          * fmin(taz(texz, nu_i, nu_c, nz, p, gam_c, gam_i, nu_r / dopp) / texz, 1.);

          F_nu_obs_x = Pmax
          * (omega * Niso / (4 * PI)) * pow(dopp, 2)
          * elem_fluxz(nu_i, nu_c, p, nu_x / dopp)
          * fmin(taz(texz, nu_i, nu_c, nz, p, gam_c, gam_i, nu_x / dopp) / texz, 1.);
        }

        if (max_f_nu_obs_radio < F_nu_obs_radio){
          max_f_nu_obs_radio = F_nu_obs_radio;
        }
        if (max_f_nu_obs_r < F_nu_obs_r){
          max_f_nu_obs_r = F_nu_obs_r;
        }
        if (max_f_nu_obs_x < F_nu_obs_x){
          max_f_nu_obs_x = F_nu_obs_x;
        }
        n_step = n_step + 1;
      }
      if(max_f_nu_obs_radio / microJy > 15.){
        n_radio = n_radio + 1.;
      }
      if(max_f_nu_obs_r / microJy > 9.0e-1){
        n_r = n_r + 1.;
      }
      if(max_f_nu_obs_x / microJy > 1.0e-4){
        n_x = n_x + 1.;
      }
      if((max_f_nu_obs_x / microJy > 1.0e-4)
        || (max_f_nu_obs_r / microJy > 9.0e-1)
        || (max_f_nu_obs_radio / microJy > 15.)){
        n_m = n_m + 1.;
      }
    }
    fprintf(stdout, "%e %i %e %e %e %e %e %e %e %e\n",
        logn, n_pop, n_gw, n_radio, n_r, n_x, n_radio / n_gw, n_r / n_gw, n_x / n_gw, n_m / n_gw);
  }
}
