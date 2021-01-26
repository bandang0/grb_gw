#ifndef PHYSICS_H
#define PHYSICS_H

#include <math.h>
#include "units.h"

/* Doppler factor in spherical coordinates */
double doppler(double g, double theta, double phi, double tv){
  double b = sqrt(1. - 1./pow(g, 2));
  double cx = cos(theta) * cos(tv) + sin(tv) * sin(theta) * cos(phi);

  return 1./(g * (1. - b * cx));
}

/* Power-law peak fluxes etc. */
double pl_peak_flux(double d, double n, double ee, double eb, double e,
double tj, double tv){
  return 8600. * (e/1.0e52)*pow(tj/0.1, 2)*pow(n/1.0e-3, 0.8)*pow(ee/0.1, 1.2)
               * pow(eb/1.0e-3, 0.8)*pow(d/(100. * Mpc), -2)
               * pow(fmax(tj/0.1, tv/0.1), -4.4) * microJy;
}

double pl_peak_time_ex(double d, double n, double ee, double eb, double e,
double tj, double tv){
  return 60. * pow(e/1.0e52, 1./3)*pow(tj/0.1, 2./3)*pow(n/1.0e-3, -1./3)
               * pow(tv/0.35, 2.) * day;
}

double pl_peak_time_no_ex(double d, double n, double ee, double eb, double e,
double tj, double tv){
  return 137. * pow(e/1.0e52, 1./3)*pow(tv/0.35, 8./3)*pow(n/1.0e-3, -1./3) * day;
}

double pl_mu(double d, double n, double ee, double eb, double e,
double tj, double tv){
  if (tv > 1.) {return 0.;}
  return 9.15e-2 * pow(d/50, -1.)*sqrt(1 - tv*tv)
         * sin(tv) / (1 - sqrt(1 - tv * tv) * cos(tv));
}

/* Relativistic deceleration */
double rel_dec_gamma(double r_dec, double gamma_0, double r){
  double eta = pow(r / r_dec, 3);
  return 0.5 * gamma_0 * (-1.0 + sqrt(1 + 4 * eta * (1 + eta / pow(gamma_0, 2))))
          / eta;
}

/* Comoving spectrum from a single electron */
double elem_fluxz(double nu_i, double nu_c, double p, double nu){
  double ef;
  if (nu_c < nu_i){
    if (nu < nu_c){
      ef = pow(nu / nu_c, 1. / 3);
    }
    else if (nu_c < nu && nu < nu_i){
      ef = pow(nu / nu_c, - 0.5);
    }
    else {
      ef = pow(nu / nu_i, - p / 2) * pow(nu_c / nu_i, 0.5);
    }
  }
  else {
      if (nu < nu_i){
        ef = pow(nu / nu_i, 1. / 3);
      }
      else if (nu_i < nu && nu < nu_c){
        ef = pow(nu / nu_i, - (p - 1.) / 2);
      }
      else {
        ef = pow(nu / nu_c, - p / 2) * pow(nu_c / nu_i, - (p - 1) / 2.);
      }
    }
    return ef;
  }

/* Shock frame self-absorption time scale 3.150, 3.151*/
double taz(double texz, double nu_i, double nu_c, double nz, double p,
            double gc, double gi, double nu){
              double pl;
              double s;
              if (gc < gi){
                s = 4 * PI * pow(nu_c, 3) / (pow(cLight, 3) * nz);
                if (nu < nu_c){
                  pl = pow(nu / nu_c, 5. / 3);
                }
                else if (nu_c < nu && nu < nu_i){
                  pl = pow(nu / nu_c, 3);
                }
                else{
                  pl = pow(nu_i / nu_c, 3) * pow(nu / nu_i, (p + 5) / 2);
                }
                return texz * s * pl;
              }
              else{
                s = (3 * p + 2) * (8 * PI) * pow(nu_i, 3) * gc
                    / (4 * (p + 2) * nz * pow(cLight, 3) * gi);
                if (nu < nu_i){
                  pl = pow(nu / nu_i, 5. / 3);
                }
                else if (nu_i < nu && nu < nu_c){
                  pl = pow(nu / nu_i, (p + 4) / 2);
                }
                else{
                  pl = pow(nu_c / nu_i, (p + 4) / 2)
                      * pow(nu / nu_c, (p + 5) / 2);
                }
                return texz * s * pl;
              }
            }

/* Shock frame self-absorption frequency 3.152, 3.153*/
double calc_nu_a(double texz, double nu_i, double nu_c, double nz, double p,
  double gc, double gi){
  double d;
  if (gc < gi){
    d = taz(texz, nu_i, nu_c, nz, p, gc, gi, nu_c) / texz;
    if (d > 1.){
      return nu_c * pow(d, - 3. / 5);
    }
    else if (1. > d && d > pow(nu_c / nu_i, 3)){
      return nu_c * pow(d, - 1. / 3);
    }
    else{
      return nu_c * pow(d, - 2 / (p + 5)) * pow(nu_i / nu_c, (p - 1) / (p + 5));
    }
  }
  else{
    d = taz(texz, nu_i, nu_c, nz, p, gc, gi, nu_i) / texz;
    if (d > 1.){
      return nu_i * pow(d, - 3. / 5);
    }
    else if (1. > d && d > pow(nu_i / nu_c, (p + 4) / 2)){
      return nu_i * pow(d, - 2. / (p + 4));
    }
    else{
      return nu_i * pow(d, - 2 / (p + 5)) * pow(nu_c / nu_i, 1 / (p + 5));
    }
  }
}

/* KNAG peak time */
double knagtp(double n){
  return 8.5 * yr * pow(n / 0.01, - 1./3.);
}

/* KNAG peak flux */
double knagfp(double eb, double n, double d){
  return 115. * microJy * pow(eb * n / 1.e-5, 0.8) * pow(d / (32 * Mpc), -2);
}

/* KN peak magnitudes */
double knrmag(double tv, double d){
  return -17.8 + sin(tv) * 4.25 + 5. * log10(d / (10 * pc));
}

double kngmag(double tv, double d){
  return -19.0 + sin(tv) * 7. + 5. * log10(d / (10 * pc));

}

double kngmagwall(double tv, double d){
  if (tv < 60 * Deg){
    return -17.6 + 7 * (1 - cos(tv)) / 0.5 + 5. * log10(d / (10 * pc));
  }
  else{
    return -17.6 + 7 + 5. * log10(d / (10 * pc));
  }
}

double knrmagwall(double tv, double d){
  if (tv < 60 * Deg){
    return -16.9 + 4 * (1 - cos(tv)) / 0.5 + 5. * log10(d / (10 * pc));
  }
  else{
    return -16.9 + 4 + 5. * log10(d / (10 * pc));
  }
}

double knimagwall(double tv, double d){
  if (tv < 60 * Deg){
    return -17.0 + 3.5 * (1 - cos(tv)) / 0.5 + 5. * log10(d / (10 * pc));
  }
  else{
    return -17.0 + 3.5 + 5. * log10(d / (10 * pc));
  }
}

double knzmagwall(double tv, double d){
  if (tv < 60 * Deg){
    return -16.8 + 2.5 * (1 - cos(tv)) / 0.5 + 5. * log10(d / (10 * pc));
  }
  else{
    return -16.8 + 2.5 + 5. * log10(d / (10 * pc));
  }
}

double peak_flux_kn(double phi, double d, double tv){
  return 8.6 * phi * pow(d / (100 * Mpc), -2.)
   * pow(fmax(1, tv / 0.1), -4.4) * mJy;
}

/* apparent proper motion */
double proper_motion(double g, double tv, double d){
  double b = sqrt(1 - pow(g, -2));
  return cLight * b * sin(tv) / (d * (1 - b * cos(tv)));
}

#endif
