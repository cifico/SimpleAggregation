#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>
#include "simul.h"

int parseArgs(int argc, char **argv, long *k_max, long *n_i, long *a_fin, long *n_simuls, char *fname, int *n_repet, int *n_times, double **tfins) {
	if (argc < 7) {
		fprintf(stderr,
				"Usage: %s k_max n(1,0) a_fin n_simuls fname n_repet tfin1"
			    " [tfin2 ...]\n",
		        argv[0]);
		return 1;
	}
	*k_max = atol(argv[1]);
	*n_i = atol(argv[2]);
	*a_fin = atol(argv[3]);
	*n_simuls = atol(argv[4]);
	strcpy(fname, argv[5]);
	*n_repet = atoi(argv[6]);

	*n_times = argc - 7;
	*tfins = malloc((*n_times) * sizeof(double));
	for (long k = 0 ; k < *n_times ; ++k) {
		(*tfins)[k] = atof(argv[7 + k]);
	}
	return 0;
}

void export(long * observables[N_OBS], char *fname, long k_max, long n_simul_tot, char *header, int R) {
	FILE *file = NULL;
	file = fopen(fname, "w+");
	//assert(file);
	fprintf(file, "%sR=%d \t", header, R);
	fprintf(file, "# k n(k,t)");
	fprintf(file, "n(k ,t) n(a,T) ");
	fprintf(file, "n(k ,t) n(a,T)^2 ");
	fprintf(file, "n(k ,t) n(a,T)^3 ");
	fprintf(file, "\n");

	double fac = 1.0 / n_simul_tot / R;

	for (long k = 0 ; k < k_max ; ++k) {
		fprintf(file, "%ld ", k);
		for (int i = 0 ; i < N_OBS ; ++i) {
			fprintf(file, " %.10e", fac * observables[i][k]);
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
	long n_i, a_fin;
	double *tfins;
	char fname[200];
	int status, n_repet;

	// Initializations
	status = parseArgs(argc, argv, &k_max, &n_i, &a_fin, &n_simuls, fname, &n_repet, &n_times, &tfins);
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
		sprintf(header, "# k_max=%ld, n_i=%ld, a=%ld, T=%lf"
				"n_simuls=Rx%dx%ld, fname=%s, R=%d",
			    k_max, n_i, a_fin, tfins[n_times - 1], world_size, n_simuls, fname, n_repet);

		printf("%s\nTimes: ", header);
		for (int k = 0 ; k < n_times ; ++k) {
			printf("%g ", tfins[k]);
		}
		printf("\nNumber of processors: %d\n", world_size);
	}

	const int n_obs_tot = N_OBS * n_times;
	long **observables = (long **) malloc(n_obs_tot * sizeof(long *));

	// It is crucial to allocate the array of pointers on all processes
	long **obs_sum = (long **) malloc(n_obs_tot * sizeof(long *));

	assert(observables);
	assert(obs_sum);

	for (int k = 0 ; k < n_obs_tot ; ++k) {
		observables[k] = (long *) calloc(k_max, sizeof(long));
		assert(observables[k]);
		
		if (world_rank == 0) {
			obs_sum[k] = (long *) calloc(k_max, sizeof(long));
			assert(obs_sum[k]);
		}
	}
	// End of initializations


	for (int r = 0 ; r < n_repet ; r++)
	{
		// Run the simulation
		// Generate the seed
		unsigned long local_seed = seedgen2();

		run(k_max, n_i, a_fin, n_simuls, n_times, tfins, observables, local_seed);

		// We sum all the observables from all the processes for position
		for (int k = 0 ; k < n_obs_tot ; ++k) {
			MPI_Reduce(observables[k], obs_sum[k], k_max, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		}

		if (world_rank == 0) {
			for (int k = 0 ; k < n_times ; ++k) {
				char hd[500], ff[100];
				sprintf(hd, "%s, t=%g\n", header, tfins[k]); // Header
				sprintf(ff, "%s_%g.dat", fname, tfins[k]); // Filename
				export(&obs_sum[k * N_OBS], ff, k_max, n_simuls * world_size, hd, r+1);
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}
	

	MPI_Finalize();

	// Free memory
	free(tfins);
	for (int k = 0 ; k < n_obs_tot ; ++k) {
		free(observables[k]);
		if (world_rank == 0)
			free(obs_sum[k]);
	}

	free(observables);
	free(obs_sum);

	return 0;
}
