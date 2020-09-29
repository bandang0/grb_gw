#include <stdio.h>
#include <libconfig.h>
#include "helpers.h"
#include "physics.h"
#include "units.h"

int main(int argc, char ** argv){

  /* Jet configuration */
  double eiso, gamma_0, n, p, epsilon_e, epsilon_b, tv, tj, d, last_day,
          first_day, nu_obs;

  /* Read configuration */
  config_t cfg;
  config_setting_t *setting;

  config_init(&cfg);
  config_read_file(&cfg, argv[1]);

  config_lookup_float(&cfg, "eiso", &eiso);
  config_lookup_float(&cfg, "gamma_0", &gamma_0);
  config_lookup_float(&cfg, "nu_obs", &nu_obs);
  config_lookup_float(&cfg, "n", &n);
  config_lookup_float(&cfg, "p", &p);
  config_lookup_float(&cfg, "epsilon_e", &epsilon_e);
  config_lookup_float(&cfg, "epsilon_b", &epsilon_b);
  config_lookup_float(&cfg, "tv", &tv);
  config_lookup_float(&cfg, "tj", &tj);
  config_lookup_float(&cfg, "d", &d);
  config_lookup_float(&cfg, "last_day", &last_day);
  config_lookup_float(&cfg, "first_day", &first_day);


  /* Global variables */
  double r_dec = pow(3 * eiso
    / (4 * PI * n * mProton * pow(gamma_0, 2) * pow(cLight, 2)), 1. /3 );
  double r_nr = r_dec * pow(gamma_0, 2. / 3);
  double miso = eiso / (pow(cLight, 2) * gamma_0);
  double omega = 2 * PI * (1 - cos(tj));
//  FILE * output = fopen(argv[1]);

  /* Spatial steps */
  double dr0 = r_dec / 2048;
  double inc_r = pow(2048, 1. / 4096);

  /* Physical data */
  double r, dr, last_r;
  int n_step;
  double t_obs;
  double t;
  double texz;
  double Niso, nz;
  double b;
  double rho_s;
  double eps_s;
  double beta;
  double gam;
  double F_nu_obs;
  double Pmax;
  double nu_a, nu_i, nu_c;
  double gam_i;
  double gam_c;
  double dopp;

  /* Initial conditions */
  texz = 0.0;
  t_obs = 0.0;
  t = 0.0;
  n_step = 1;
  r = 0.0;
  last_r = 0.0;

  /* Principal loop */
  while (r < 10 * r_nr && t_obs / day < last_day){

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

        if (t_obs / day > first_day){
          fprintf(stdout,
            "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
            r, t, t_obs, gam, rho_s, eps_s, Niso,
            epsilon_e, epsilon_b, p, beta, nu_obs, F_nu_obs,
            nu_a, nu_i, nu_c, gam_i, gam_c, dopp);
        }

        n_step = n_step + 1;
      }
  //  fclose(output)
}
