#ifndef SIMUL_H
#define SIMUL_H

#include "mkl.h"
#include "mkl_vsl.h"

#define N_OBS 15

#define BATCH_SIZE 256
// See https://software.intel.com/en-us/mkl-developer-reference-c-brng-parameter-definition
#define CUSTOM_RNG VSL_BRNG_SFMT19937

#define MIN(a,b) (((a)<(b))?(a):(b))

typedef struct {
	long flux;
	long *positions;
	int *occupations;
	long winding;
} State;

/*void run(const long n_sites, const long n_parts, const double prob,
		 const double tfin, const long n_simuls, const long n_times,
		 long **observables, const unsigned long seed, const int determ);*/
void run(const long n_sites, double rhop, double rhom, const double prob,
		 const long n_simuls, const int n_times, const double *tfins, 
		 long **observablesX, 
#ifdef FLUX		 
		 long **observablesQ, 
#endif
		 const unsigned long seed,
		 const int determ);

long init_random(State *state, const long n_sites, const double rhop, const double rhom,
		         VSLStreamStatePtr stream);
long init_determ(State *state, const long n_sites, const double rhop, const double rhom);
double evolve(State *state, const long n_sites, long n_parts,
		      const double prob, const double biaiseq, const double tini, const double tfin,
			  VSLStreamStatePtr stream);
void updateObsX(State *state, const long n_sites, long **observables);
#ifdef FLUX		 
void updateObsQ(State *state, const long n_sites, long **observables);
#endif

#endif
