#include <stdio.h>
#include <math.h>
#include "simul.h"

void run(const long n_sites, double rhop, double rhom, const double prob,
		 const long n_simuls, const int n_times, const double *tfins, 
		 long **observablesX, 
#ifdef FLUX		 
		 long **observablesQ, 
#endif
		 const unsigned long seed, const int determ) {
	
	const int n_obs_tot = N_OBS * n_times;

	VSLStreamStatePtr stream;
	vslNewStream(&stream, CUSTOM_RNG, seed);

	long n_parts = 0;
	//equilibrate the densities at junction on the ring
	double biaiseq = rhom*(1-rhop)/(rhop+rhom-2*rhom*rhop);

	State state;
	//n_sites front + n_sites back
	state.occupations = malloc(2 * n_sites * sizeof(int));
	state.positions = malloc(2 * n_sites * sizeof(long));

	for (long i = 0 ; i < n_simuls ; ++i) {
		// Deterministic / random initialization
		if (determ) {		
			n_parts = init_determ(&state, n_sites, rhop, rhom);
		} 
		else {
			n_parts = init_random(&state, n_sites, rhop, rhom, stream);
		}
		double ti = 0.0, tf = 0.0;
		for (int k = 0 ; k < n_times ; ++k) {
			tf = tfins[k];
			ti = evolve(&state, n_sites, n_parts, prob, biaiseq, ti, tf, stream);
			// Shift observables depending on index of final time
			updateObsX(&state, n_sites, &observablesX[k * N_OBS]);
#ifdef FLUX	
			updateObsQ(&state, n_sites, &observablesQ[k * N_OBS]);
#endif
		}
	}

	free(state.positions);
	free(state.occupations);

	vslDeleteStream(&stream);
}

long init_random(State *state, const long n_sites, const double rhop, const double rhom,
		         VSLStreamStatePtr stream) {
	
	int bern[n_sites];
	//the index of the current particle
	long numpart = 0;
	
#ifdef FLUX	
	state->flux = 0;
#endif

	// The TP
	//state->positions[0] = 0;
	//state->occupations[0] = 1;
	state->winding = 0;

	// Random Bernoulli of parameter rhop
	viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF, stream, n_sites, bern, rhop);

	for (long i = 0 ; i < n_sites ; ++i) {			 
		if (bern[i] == 1)
		{
			state->occupations[i] = 1;
			state->positions[numpart] = i;			
			numpart += 1;
		}
		else
		{
			state->occupations[i] = 0;
		}
	}

	// Random Bernoulli of parameter rhom
	viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF, stream, n_sites, bern, rhom);

	for (long i = n_sites ; i < 2*n_sites ; ++i) {			 
		if (bern[i-n_sites] == 1)
		{
			state->occupations[i] = 1;
			state->positions[numpart] = i;			
			numpart += 1;
		}
		else
		{
			state->occupations[i] = 0;
		}
	}

	//the total number of particles (tracer included)
	return numpart;	
}

long init_determ(State* state, const long n_sites, const double rhop, const double rhom) {
	// Different strategies depending whether the system is more or less than
	// half full.
	long n_parts = n_sites * rhop;

	if (n_sites >= 2 * n_parts) {
		for (long i = 0; i < n_sites; ++i) {
			state->occupations[i] = 0;
		}

		// Place the particles equidistantly
		const long step = n_sites / n_parts; // At least 2
		for (long i = 0; i < n_parts; ++i) {
			long s = i * step;
			state->positions[i] = s;
			state->occupations[s] = 1;
		}

	}
	else {
		for (long i = 0; i < n_sites; ++i) {
			state->occupations[i] = 1;
		}

		const long n_vacs = n_sites - n_parts;
		const long step = n_sites / n_vacs; // At least 2
		const long hstep = step / 2;
		for (long i = 0; i < n_vacs; ++i) {
			long s = i * step + hstep;
			state->occupations[s] = 0;
		}

		long j = 0;
		for (long i = 0; i < n_sites; ++i) {
			if (state->occupations[i])
				state->positions[j++] = i;
		}
	}

	long n_parts2 = n_sites * rhom;

	if (n_sites >= 2 * n_parts2) {
		for (long i = n_sites; i < 2 * n_sites; ++i) {
			state->occupations[i] = 0;
		}

		// Place the particles equidistantly
		const long step = n_sites / n_parts2; // At least 2
		for (long i = 0; i < n_parts2; ++i) {
			long s = i * step;
			state->positions[n_parts + i] = 2 * n_sites - 1 - s;
			state->occupations[2 * n_sites - 1 - s] = 1;
		}

	}
	else {
		for (long i = n_sites; i < 2 * n_sites; ++i) {
			state->occupations[i] = 1;
		}

		const long n_vacs = n_sites - n_parts2;
		const long step = n_sites / n_vacs; // At least 2
		const long hstep = step / 2;
		for (long i = 0; i < n_vacs; ++i) {
			long s = i * step + hstep;
			state->occupations[2 * n_sites - 1 - s] = 0;
		}

		long j = n_parts;
		for (long i = n_sites; i < 2 * n_sites; ++i) {
			if (state->occupations[i])
				state->positions[j++] = i;
		}
	}
	state->winding = 0;
#ifdef FLUX	
	state->flux = 0;
#endif
	return n_parts + n_parts2;
}

double evolve(State *state, const long n_sites, long n_parts,
		      const double prob, const double biaiseq, const double tini, const double tfin,
			  VSLStreamStatePtr stream) {
	const double biais = prob - 0.5;
	//potentially the total jump rate is greater than n_parts because of the biais at junction
	const double lambda = (tfin - tini) * (n_parts + 1);
	const long endsite = 2*n_sites - 1;
	int indices[BATCH_SIZE];
	int bern[BATCH_SIZE];
	double us[BATCH_SIZE];

	long part = 0;
	long pos = 0;
	long pos_next = 0;

	int n_iters = 0; // Number of time iterations (Poisson distributed)
	viRngPoisson(VSL_RNG_METHOD_POISSON_PTPE, stream, 1, &n_iters, lambda);

	while (n_iters) {
		int n_i = MIN(n_iters, BATCH_SIZE);

		// Generation of the random numbers in batches
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_i,
					     us, 0.0, 1.0);
		viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_i, indices,
					 0, n_parts+1);
		
		for (int i = 0 ; i < n_i ; ++i) {
			if(indices[i] == n_parts)
			{
				//where to add jump rate at junction
				pos = n_sites - (biaiseq > 0.5);
				if (state->occupations[pos] == 0)
				{
					//if there is no particle at junction, go to next step
					continue;
				}
				else if (us[i] > fabs(biaiseq - 0.5))
				{
					//we add the part of the jump rate greater than 0.5
					//up to 0.5 it is taken into account normally
					continue;
				}
				part = 0; 
				while(state->positions[part] != pos)
				{
					part++;
				}			
				pos_next = n_sites + (biaiseq > 0.5) - 1;
			}
			else
			{
				part = indices[i];
				pos = state->positions[part];

				// Next position
				
				double pp = 0.5 + (part == 0) * biais;
				pos_next = pos + 2 * (us[i] < pp) - 1;


				if (pos_next == 2*n_sites)
				{
					pos_next = 0;
				}
				else if (pos_next == -1)
				{
					pos_next = endsite;
				}
				//At the junction, rate biaiseq to the right and 1-biaiseq to the left, biased flux
				else if (pos_next == n_sites && pos == n_sites - 1 && us[i] > biaiseq)
				{
					//rate biaiseq
					continue;				
				}
				else if (pos_next == n_sites - 1 && pos == n_sites && us[i] <= biaiseq)
				{
					//rate 1-biaiseq
					continue;
				}
			}

			int o = state->occupations[pos_next];
			long dpos = (1 - o) * (pos_next - pos);
			state->positions[part] += dpos; // Move if next site not occupied
			state->occupations[pos] = o; // Occupied only if doesn't move
			state->occupations[pos_next] = 1; // Next site is occupied anyways

			// Add or substract one if the TP did one turn
			long dw = (part == 0) * (1 - o) * (
			           (pos == endsite && pos_next == 0)
					   - (pos == 0 && pos_next == endsite) 
					  );
			state->winding += dw;
			/*if (dw) {
				printf("%ld: %ld, %ld, %ld\n", dw, part, pos, pos_next);
			}*/


#ifdef FLUX	
			long df = dpos * ((pos == 0 && pos_next == 1) + (pos == 1 && pos_next == 0));
			state->flux += df;
#endif
				
		} 

		n_iters -= n_i;
	}
	// This may introduce a small error
	return tfin;
}


void updateObsX(State *state, const long n_sites, long **observables) {
	long X = state->positions[0];
	// long n_sites_2 = n_sites / 2;
	// long xPer = (X + n_sites_2) % n_sites - n_sites_2;
	long xPer = X + (state->winding * n_sites * 2);

	long ep = state->occupations[(X + 1) % (2*n_sites)];
	long em = state->occupations[(X + 2*n_sites - 1) % (2*n_sites)];
	long xPerEp = xPer * ep;
	long xPerEm = xPer * em;

	long xPer2 = xPer * xPer;
	long xPer2Ep = xPer2 * ep;
	long xPer2Em = xPer2 * em;
	long xPer3 = xPer2 * xPer;
	long xPer3Ep = xPer3 * ep;
	long xPer3Em = xPer3 * em;

	#pragma ivdep // Go icc, you can vectorize it (no dependancy)
	for (long i = 0 ; i < n_sites*2 ; ++i) {
		long k = (X + i) % (2*n_sites);
		long er = state->occupations[k];

		observables[0][i] += er;
		observables[1][i] += ep * er;
		observables[2][i] += em * er;
		observables[3][i] += xPer * er;
		observables[4][i] += xPerEp * er;
		observables[5][i] += xPerEm * er;
		observables[6][i] += xPer;
		observables[7][i] += xPer2 * er;
		observables[8][i] += xPer2Ep * er;
		observables[9][i] += xPer2Em * er;
		observables[10][i] += xPer2;
		observables[11][i] += xPer3 * er;
		observables[12][i] += xPer3Ep * er;
		observables[13][i] += xPer3Em * er;
		observables[14][i] += xPer3;
	}
}

#ifdef FLUX	
void updateObsQ(State *state, const long n_sites, long **observables) {
	long xPer = state->flux;
	// long n_sites_2 = n_sites / 2;
	// long xPer = (X + n_sites_2) % n_sites - n_sites_2;

	long ep = state->occupations[1];
	long em = state->occupations[0];
	long xPerEp = xPer * ep;
	long xPerEm = xPer * em;

	long xPer2 = xPer * xPer;
	long xPer2Ep = xPer2 * ep;
	long xPer2Em = xPer2 * em;
	long xPer3 = xPer2 * xPer;
	long xPer3Ep = xPer3 * ep;
	long xPer3Em = xPer3 * em;

	#pragma ivdep // Go icc, you can vectorize it (no dependancy)
	for (long i = 0 ; i < n_sites*2 ; ++i) {
		long k = (i + 1) % (2*n_sites);
		long er = state->occupations[k];

		observables[0][i] += er;
		observables[1][i] += ep * er;
		observables[2][i] += em * er;
		observables[3][i] += xPer * er;
		observables[4][i] += xPerEp * er;
		observables[5][i] += xPerEm * er;
		observables[6][i] += xPer;
		observables[7][i] += xPer2 * er;
		observables[8][i] += xPer2Ep * er;
		observables[9][i] += xPer2Em * er;
		observables[10][i] += xPer2;
		observables[11][i] += xPer3 * er;
		observables[12][i] += xPer3Ep * er;
		observables[13][i] += xPer3Em * er;
		observables[14][i] += xPer3;
	}
}
#endif
