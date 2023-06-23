#include "decs.h"
#include "bl_coord.h"
void cool_electrons_1zone(struct FluidState *S, int i, int j)

void cool_electrons(struct FluidState *S)
{
  #pragma omp parallel for collapse(2)
  ZLOOP {
    cool_electrons_1zone(S, i, j);
  }
}

void cool_electrons_1zone(struct FluidState *S, int i, int j)
{
  //to find r and ut:
  // double r = G->r
  double X[NDIM], r, th, ucon[NDIM], trans[NDIM][NDIM], tmp[NDIM];
  double AA, BB, CC, discr;
  double alpha, gamma, beta[NDIM];
  struct blgeom;
  struct of_geom blgeom;
  coord(i, j, CENT, X);
  bl_coord(X, &r, &th);
  blgset(i, j, &blgeom);
  ucon[1] = S->P[U1][j][i];
  ucon[2] = S->P[U2][j][i];
  ucon[3] = S->P[U3][j][i];
  AA = blgeom.gcov[0][0];
  BB = 2.*(blgeom.gcov[0][1]*ucon[1] +
           blgeom.gcov[0][2]*ucon[2] +
           blgeom.gcov[0][3]*ucon[3]);
  CC = 1. +
      blgeom.gcov[1][1]*ucon[1]*ucon[1] +
      blgeom.gcov[2][2]*ucon[2]*ucon[2] +
      blgeom.gcov[3][3]*ucon[3]*ucon[3] +
      2.*(blgeom.gcov[1][2]*ucon[1]*ucon[2] +
          blgeom.gcov[1][3]*ucon[1]*ucon[3] +
          blgeom.gcov[2][3]*ucon[2]*ucon[3]);
  discr = BB*BB - 4.*AA*CC;
  ucon[0] = (-BB - sqrt(discr))/(2.*AA);
  double ut = ucon[0];

  //m is arbitrary and we can just set is to 3 because it just needs tobe larger than 1
  double m = 3.0;
  
  //dt is a global variable so we don't even need to initialize it

  //now we just need to update the internal energy of the fluid at (i,j):
  S->P[UU][j][i] += -1*S->P[UU][j][i]/(m*pow(r, 1.5)*ut)*dt;
}
