#ifndef SIMUL_H
#define SIMUL_H

#include "mkl.h"
#include "mkl_vsl.h"

#define N_OBS 4

#define BATCH_SIZE 256
// See https://software.intel.com/en-us/mkl-developer-reference-c-brng-parameter-definition
#define CUSTOM_RNG VSL_BRNG_SFMT19937

#define MIN(a,b) (((a)<(b))?(a):(b))

typedef struct {
	long n_agg;
	long *aggregate;
} State;

/*void run(const long n_sites, const long n_parts, const double prob,
		 const double tfin, const long n_simuls, const long n_times,
		 long **observables, const unsigned long seed, const int determ);*/
void run(const long k_max, const long n_i, const long a_fin, const long n_simuls, const long n_times, const double *tfins, long **observables, const unsigned long local_seed);

double evolve(State *state, const long n_i, const double tini, const double tfin, VSLStreamStatePtr stream);

void updateObs(State *states, long **observables, const long a_fin, const long k_max, const long n_times);

#endif
