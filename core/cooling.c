#include "decs.h"
void cool_electrons_1zone(struct FluidState *S, int i, int j);

void cool_electrons(struct GridGeom *G, struct FluidState *S)
{
  #pragma omp parallel for collapse(2)
  ZLOOP {
    cool_electrons_1zone(G, S, i, j);
  }
}

void cool_electrons_1zone(struct GridGeom *G, struct FluidState *S, int i, int j)
{
  //to find ut:
  ucon_calc(G, S, i, j, CENT);
  double ut = S->ucon[0][j][i];

  // to find r:
  double X[NDIM];
  coord(i, j, CENT, X);
  double r, th;
  bl_coord(X, &r, &th);

  //m is arbitrary
  double m = 3.0;
  
  //dt is a global variable so we don't even need to initialize it

  //to fing uel:
  double uel = pow(S->P[RHO][j][i],game)*exp(S->P[KEL0][j][i]*(game-1));

  //update the internal energy of the electrons at (i,j):
  uel -= uel/(m*pow(r, 1.5)*ut)*dt/2;

  //update the entropy with the new internal energy
  S->P[KEL0][j][i] = log(uel/pow(S->P[RHO][j][i], game))/(game-1);
}
