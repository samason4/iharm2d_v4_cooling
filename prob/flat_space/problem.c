/*---------------------------------------------------------------------------------

  PROBLEM.C
 
  -Read and store problem-specific parameters from parameter file
  -Initialize flat space stuff

-----------------------------------------------------------------------------------*/

#include "decs.h"
#include <complex.h>

// Local declarations
#define FML_DBL_OUT "%28.18e"
#define FML_INT_OUT "%10d"
#define STRING_OUT "%15s"

static int tau;

// Assign problem param names and pointers
void set_problem_params()
{
	set_param("tau", &tau);
}

// Save problem specific details
void save_problem_data(FILE *fp)
{
	fprintf(fp, FML_INT_OUT, tau);
}

// Initializing mean state and perturbations based on nmode
void init(struct GridGeom *G, struct FluidState *S)
{
	double X[NDIM];

	// Mean (background) state
	double rho0 = 1.;
	double u0 = 1.;
	double B10 = 0.;
	double B20 = 0.;
	double B30 = 0.;
	double U10 = 0.;
	double U20 = 0.;
	double U30 = 0.;

	// Set grid
	set_grid(G);
	LOG("Grid set");

	ZLOOP
	{
		coord(i, j, CENT, X);
		S->P[RHO][j][i] = rho0;
		S->P[UU][j][i] = u0;
		S->P[U1][j][i] = U10;
		S->P[U2][j][i] = U20;
		S->P[U3][j][i] = U30;
		S->P[B1][j][i] = B10;
		S->P[B2][j][i] = B20;
		S->P[B3][j][i] = B30;
	}

	// Enforce boundary conditions
	set_bounds(G, S);

	LOG("Finished init()");
}
