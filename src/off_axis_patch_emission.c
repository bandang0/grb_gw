#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <libconfig.h>
#include "helpers.h"
#include "physics.h"
#include "units.h"

#define DL (1000. * Mpc)
#define G0 (400.)
#define B0 (sqrt(1. - 1./pow(G0, 2)))
#define D0 (1/(G0 * (1 - B0)))
#define TJ (0.1)
#define ALPHA (8.)
#define BETA (10.)
#define EPON (300. * keV)
#define LPON (1.e52)
#define TOBS1ON (10.)
#define TPON (1.0)
#define LPP (LPON / (D0 * D0))
#define EPP (EPON / D0)
#define TVAR (B0 * TOBS1ON / (G0 * G0 * (1 - B0)))
#define T1 (TOBS1ON + G0 * G0 * TVAR)
#define T2 (T1 + TPON/(1 - B0))
#define EG (300.0 * keV)
#define EX (5.0 * keV)
#define NUG (EG/hPlanck)
#define NUX (EX/hPlanck)
#define NUOBS (NUX)
#define NT (1.e5)
#define NP (1.e3)

double gam(double theta){
  if (theta > PI / 2.){
    return gam(PI - theta);
  }
  if (theta < TJ){
    return G0;
  }
  else{
    return 1. + (G0 - 1.) * pow(theta / TJ, - BETA);
  }
}

double dLdO(double nup, double theta, double t){
  if (theta > PI / 2.){
    return dLdO(nup, PI - theta, t);
  }
  return (LPP / EPP)
  * (theta < TJ ? 1. : pow(theta / TJ, - ALPHA))
  * pow(hPlanck * nup / EPP, hPlanck * nup < EPP ? 1.1 : -2.3)
  * (fabs(t - 0.5 * (T1 + T2)) < 0.5 * fabs(T2 - T1) ? 1. : 0.);

}

double flux(double tobs, double nuobs, double tv){
  double theta, phi, g, d;
  double flu = 0.;
  int i, j;

  for (i = 0; i < NT; i++){
    theta = i * PI / NT;
    for (j = 0; j < NP; j++){
      phi = j * 2 * PI / NP;
      g = gam(theta);
      d = doppler(g, theta, phi, tv);
      flu = flu + sin(theta)
      * dLdO(nuobs/d, theta, g * d * tobs)
      * pow(d, 2);
    }
  }
  return flu * 2. * PI * PI / (NP * NT);
}

double fluxopti(double tobs, double nuobs, double tv){
  double theta, phi, dtheta, g, d;
  double flux = 0.;
  double dphi = 2. * PI / NP;

  while (theta < PI){
     if (theta < TJ){
       dtheta = acos(cos(theta) - (1. - cos(TJ))/NT) - theta;
     }
     else{
       dtheta = fmin(PI - theta,
          (1. - cos(TJ))/(NT * sin(theta) * pow(theta/TJ, - ALPHA)));
     }
    for (int j = 0; j < NP; j++){
      phi = j * 2 * PI / NP;
      g = gam(theta);
      d = doppler(g, theta, phi, tv);
      flux = flux + dtheta * dphi * sin(theta)
      * dLdO(nuobs/d, theta, g * d * tobs)
      * pow(d, 2);
    }
    theta = theta + dtheta;
  }
  return flux;
}

int main(int argc, char ** argv){
  double tv = atof(argv[1]);
  double tobsmin = atof(argv[2]);
  double tobsmax = atof(argv[3]);
  double Npt = atof(argv[4]);

  MPI_Init(NULL, NULL);

  // Get the number of processes and rank
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank > 0){
    double tobscalc;
    double fluxcalc, flugcalc;
    while(1){
      // receive time from master
      MPI_Recv(&tobscalc, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
        MPI_STATUS_IGNORE);
      if (tobscalc == -1.){
        // termination signal
        break;
      }
#ifdef OPTI
      fluxcalc = fluxopti(tobscalc, NUX, tv);
      flugcalc = fluxopti(tobscalc, NUG, tv);
#else
      fluxcalc = flux(tobscalc, NUX, tv);
      flugcalc = flux(tobscalc, NUG, tv);
#endif
      // send time and flux back
      MPI_Send(&tobscalc, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&fluxcalc, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&flugcalc, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
  }
  else{
    double tobssend, tobsrec, tobswrite, fluxwrite, flugwrite;
    int rts = 1;
    int rtg;
    double terminate = -1.;
    double fluxmaster, flugmaster;
    for (tobssend = tobsmin;
      tobssend < tobsmax;
      tobssend = tobssend * pow(tobsmax/tobsmin, 1/Npt)){
      if (rts == world_size){
        // all slaves have received work
#ifdef OPTI
        fluxmaster = fluxopti(tobssend, NUX, tv);
        flugmaster = fluxopti(tobssend, NUG, tv);
#else
        fluxmaster = flux(tobssend, NUX, tv);
        flugmaster = flux(tobssend, NUG, tv);
#endif
        for (rtg = 1; rtg < world_size; rtg++){
          // retreive work from slaves
          MPI_Recv(&tobswrite, 1, MPI_DOUBLE, rtg, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);
          MPI_Recv(&fluxwrite, 1, MPI_DOUBLE, rtg, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);
          MPI_Recv(&flugwrite, 1, MPI_DOUBLE, rtg, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);
          fprintf(stdout, "%e %e %e\n", tobswrite, fluxwrite, flugwrite);
        }
        fprintf(stdout, "%e %e %e\n", tobssend, fluxmaster, flugmaster);
        rts = 1;
        continue;
      }
      MPI_Send(&tobssend, 1, MPI_DOUBLE, rts, 0, MPI_COMM_WORLD);
      rts = rts + 1;
    }
    for (rts = 1; rts < world_size; rts++){
      // send termination signal
      MPI_Send(&terminate, 1, MPI_DOUBLE, rts, 0, MPI_COMM_WORLD);
    }
  }
  MPI_Finalize();
}
