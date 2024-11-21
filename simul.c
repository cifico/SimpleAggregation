#include <stdio.h>
#include <math.h>
#include "simul.h"

void run(const long k_max, const long n_i, const long a_fin, const long n_simuls, const long n_times, const double *tfins, long **observables, const unsigned long local_seed) {
	
	const int n_obs_tot = N_OBS * n_times;

	VSLStreamStatePtr stream;
	vslNewStream(&stream, CUSTOM_RNG, local_seed);

	State *states;
	states = (State *) malloc(n_times * sizeof(State));	
	for (int k = 0 ; k < n_times ; ++k) 
	{
		states[k].n_agg = n_i;
		states[k].aggregate = (long *) calloc(n_i, sizeof(long));
	}
	

	for (long i = 0 ; i < n_simuls ; ++i) {
		states[0].n_agg = n_i;
		
		for (int j = 0 ; j < n_i ; ++j) 
		{
			states[0].aggregate[j] = 1;
		}
		
		int k = 0;
		double ti = 0.0, tf = 0.0;

		tf = tfins[k];
		ti = evolve(&states[k], ti, tf, stream);	

		for (k = 1 ; k < n_times ; ++k) {
			states[k].n_agg = states[k-1].n_agg;

			for (int j = 0 ; j < states[k-1].n_agg; ++j) 
			{	
				states[k].aggregate[j] = states[k-1].aggregate[j] ;
			}

			tf = tfins[k];
			ti = evolve(&states[k], ti, tf, stream);	
		}
		
		updateObs(states, observables, a_fin, k_max, n_times); 
	}
	
	for (int k = 0 ; k < n_times ; ++k) 
	{
		free(states[k].aggregate);
	}
	free(states);

	vslDeleteStream(&stream);
}

double evolve(State *state, const double tini, const double tfin, VSLStreamStatePtr stream) {
	double times[BATCH_SIZE];
	double indices[2 * BATCH_SIZE];

	double t = tini;
	long agg1 = 0;
	long agg2 = 0;

	while (1) {

		// Generation of the random numbers in batches
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, BATCH_SIZE, times, 0.0, 1.0);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 2 * BATCH_SIZE, indices, 0.0, 1.0);
		
		for (int i = 0 ; i < BATCH_SIZE ; ++i) {
			if (state->n_agg == 1) return tfin;

			t += -2 * log(times[i]) / (state->n_agg - 1);
			if (t > tfin) return tfin;

			agg1 = indices[2 * i] * state->n_agg;
			
			agg2 = indices[2 * i + 1] * (state->n_agg - 1);

			if (agg1 == state->n_agg - 1)
			{
				state->aggregate[agg2] += state->aggregate[agg1];
			}
			else
			{
				agg2 += (agg2 >= agg1); //we want two distinct aggregates 

				state->aggregate[agg1] += state->aggregate[agg2];
				state->aggregate[agg2] = state->aggregate[state->n_agg - 1];
			}

			state->n_agg -= 1;
			//At this point aggregate contain all current aggregates between indices 0 and n_agg - 1 included
		}
	}
	return tfin;
}


void updateObs(State *states, long **observables, const long a_fin, const long k_max, const long n_times) {
	long naT = 0;
	
	for (long i = 0 ; i < states[n_times - 1].n_agg ; ++i) 
	{
		if (states[n_times - 1].aggregate[i] == a_fin) naT += 1;
	}

	long k = 0;
	long naT2 = naT * naT;
	long naT3 = naT2 * naT;

	#pragma ivdep
	for (long t = 0 ; t < n_times ; ++t) {
		for (long i = 0 ; i < states[t].n_agg ; ++i) 
		{
			k = MIN(k_max, states[t].aggregate[i]) - 1;

			//this stores the number of aggregates of size k+1 
			observables[t * N_OBS][k] += 1;
			observables[t * N_OBS + 1][k] += naT;
			observables[t * N_OBS + 2][k] += naT2;
			observables[t * N_OBS + 3][k] += naT3;
		}
	}
}

