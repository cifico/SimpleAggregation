#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>
#include "simul.h"

int parseArgs(int argc, char **argv, long *k_max, double *n_i, double *rhom,
		      double *prob, long *n_simuls, char *fname,
			  int *determ, int *n_repet, int *n_times, double **tfins) {
	if (argc < 9) {
		fprintf(stderr,
				"Usage: %s k_max n(1,0) n_simuls fname n_repet tfin1"
			    " [tfin2 ...]\n",
		        argv[0]);
		return 1;
	}
	*k_max = atol(argv[1]);
	*n_i = atof(argv[2]);
	*n_simuls = atol(argv[3]);
	strcpy(fname, argv[4]);
	*n_repet = atoi(argv[5]);

	*n_times = argc - 6;
	*tfins = malloc(*n_times * sizeof(double));
	for (long k = 0 ; k < *n_times ; ++k) {
		(*tfins)[k] = atof(argv[6 + k]);
	}
	return 0;
}

void export(long * observables[N_OBS], char *fname, long k_max,
		    long n_simul_tot, char *header, int R) {
	FILE *file = NULL;
	file = fopen(fname, "w+");
	//assert(file);
	fprintf(file, "%sR=%d \t", header, R);
	fprintf(file, "# r er e+*er e-*er");
	fprintf(file, " X*er X*e+*er X*e-*er X");
	fprintf(file, " X^2*er X^2*e+*er X^2*e-*er X^2");
	fprintf(file, " X^3*er X^3*e+*er X^3*e-*er X^3");
	fprintf(file, "\n");

	double fac = 1.0 / n_simul_tot / R;

	for (long i = 0 ; i < k_max ; ++i) {
		fprintf(file, "%ld ", i);
		for (int k = 0 ; k < N_OBS ; ++k) {
			fprintf(file, " %.10e", fac * observables[k][i]);
		}
		fprintf(file, "\n");
	}
	fclose(file);
}

void exportCums(long * obs_sum[N_OBS], char *fname, long n_simul_tot,
                long n_times, double *tfins, char *header, int R) {
	FILE *file = NULL;
	file = fopen(fname, "w+");

	fprintf(file, "%s\n", header);
	fprintf(file, "# t X X^2 X^3 X^4 X^5 X^6\n");
	double fac = 1.0 / n_simul_tot / R;

	for (long i = 0 ; i < n_times ; ++i) {
		fprintf(file, "%lf ", tfins[i]);
		for (int k = 0 ; k < 3 ; ++k) {
			fprintf(file, " %.10e", fac * obs_sum[i * N_OBS + 6 + 4*k][0]);
		}
		fprintf(file, "\n");
	}

	fclose(file);
}


// From https://arxiv.org/abs/1005.4117
/*unsigned long seedgen1(void)  {
    unsigned long s, seed, pid;
    pid = getpid();
    s = time(NULL);
    seed = abs(((s*181)*((pid-83)*359))%104729);
    return seed;
}*/

// F*** it's hard to find reliable advice on how to seed parallel generators
// Let's read from /dev/urandom...
unsigned long seedgen2(void)  {
	int fd = open("/dev/urandom", O_RDONLY);
    unsigned long seed;
	read(fd, &seed, sizeof(unsigned long));
	close(fd);
    return seed;
}

int main(int argc, char **argv) {
	long k_max, n_simuls;
	int n_times;
	double prob, n_i, rhom;
	double *tfins;
	char fname[200];
	int determ, status, n_repet;

	// Initializations
	status = parseArgs(argc, argv, &k_max, &n_i, &rhom, &prob, &n_simuls,
			           fname, &determ, &n_repet, &n_times, &tfins);
	if (status) {
		return status;
	}

	int world_rank, world_size;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	//world_rank = 0; 

	char header[500];
	if (world_rank == 0) {
		sprintf(header,
				"# k_max=%ld, n_i=%lf, rhom=%lf, prob=%lf, "
				"n_simuls=Rx%dx%ld, fname=%s, determ=%d, n_repet=%d",
			    k_max, n_i, rhom, prob, world_size, n_simuls,
				fname, determ, n_repet);
		printf("%s\nTimes: ", header);
		for (int k = 0 ; k < n_times ; ++k) {
			printf("%g ", tfins[k]);
		}
		printf("\nNumber of processors: %d\n", world_size);
	}

	const int n_obs_tot = N_OBS * n_times;
	long **observablesX = (long **) malloc(n_obs_tot * sizeof(long *));
#ifdef FLUX
	long **observablesQ = (long **) malloc(n_obs_tot * sizeof(long *));
#endif
	// It is crucial to allocate the array of pointers on all processes
	long **obs_sum = (long **) malloc(n_obs_tot * sizeof(long *));
#ifdef FLUX
	assert(observablesQ);
#endif
	assert(observablesX);
	assert(obs_sum);

	long n_sites_eff = 2*k_max;

	for (int k = 0 ; k < n_obs_tot ; ++k) {
		observablesX[k] = (long *) calloc(n_sites_eff, sizeof(long));
		assert(observablesX[k]);

#ifdef FLUX
		observablesQ[k] = (long *) calloc(n_sites_eff, sizeof(long));
		assert(observablesQ[k]);
#endif
		
		if (world_rank == 0) {
			obs_sum[k] = (long *) calloc(n_sites_eff, sizeof(long));
			assert(obs_sum[k]);
		}
	}
	// End of initializations


	for (int r = 0 ; r < n_repet ; r++)
	{
		// Run the simulation
		// Generate the seed
		unsigned long local_seed = seedgen2();

		run(k_max, n_i, rhom, prob, n_simuls, n_times, tfins, observablesX, 
#ifdef FLUX
		observablesQ,
#endif		
		local_seed, determ);

		// We sum all the observables from all the processes for position
		for (int k = 0 ; k < n_obs_tot ; ++k) {
			MPI_Reduce(observablesX[k], obs_sum[k], n_sites_eff, MPI_LONG, MPI_SUM,
					0, MPI_COMM_WORLD);
		}

		if (world_rank == 0) {
#ifdef ONLY_CUMS
			exportCums(obs_sum, fname, n_simuls * world_size, n_times,
		  	         tfins, header, r+1);
#else
			for (int k = 0 ; k < n_times ; ++k) {
				char hd[500], ff[100];
				sprintf(hd, "%s, tfin=%g\n", header, tfins[k]); // Header
				sprintf(ff, "%s_%g_X.dat", fname, tfins[k]); // Filename
				export(&obs_sum[k * N_OBS], ff, n_sites_eff, n_simuls * world_size, hd, r+1);
			}
#endif
	}

#ifdef FLUX
		// We sum all the observables from all the processes for flux
		for (int k = 0 ; k < n_obs_tot ; ++k) {
			MPI_Reduce(observablesQ[k], obs_sum[k], n_sites_eff, MPI_LONG, MPI_SUM,
					0, MPI_COMM_WORLD);
		}

		// Export
		if (world_rank == 0) {	
			for (int k = 0 ; k < n_times ; ++k) {
				char hd[500], ff[100];
				sprintf(hd, "%s, tfin=%g\n", header, tfins[k]); // Header
				sprintf(ff, "%s_%g_Q.dat", fname, tfins[k]); // Filename
				export(&obs_sum[k * N_OBS], ff, n_sites_eff, n_simuls * world_size, hd, r+1);
			}
	
		}
#endif

		MPI_Barrier(MPI_COMM_WORLD);
	}
	

	MPI_Finalize();

	// Free memory
	free(tfins);
	for (int k = 0 ; k < n_obs_tot ; ++k) {
		free(observablesX[k]);
#ifdef FLUX
		free(observablesQ[k]);
#endif
		if (world_rank == 0)
			free(obs_sum[k]);
	}
	free(observablesX);
#ifdef FLUX
	free(observablesQ);
#endif
	free(obs_sum);

	return 0;
}
