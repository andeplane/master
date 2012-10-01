#include <iostream>
#include "math.h"
#include "time.h"
#include <armadillo>
#include <fstream>
#include "System.h"
#include "StatisticsSampler.h"
#include <omp.h>
#include "defines.h"

using namespace std;

int main(int args, char* argv[]) {
	if(args == 2 && atoi(argv[1]) == -1) {
		printf("Run with:\n\n");
		printf("dsmc printPositions[0|1] N[int] T[double] timesteps[int]\n");
		printf("\n");
		return 0;
	}

	bool printPositions = args > 1 ? atoi(argv[1]) : false;
	int N = args > 2 ? atoi(argv[2]) : 10000;
	int T = args > 3 ? atof(argv[3]) : 3;
	int timesteps = args > 4 ? atof(argv[4]) : 1000;

	System *system = new System(N,T);
	StatisticsSampler *sampler = new StatisticsSampler(system, timesteps);
	
	FILE *positions = 0;
	if(printPositions) {
		positions = fopen("pos.xyz","w");
		system->printPositionsToFile(positions);
	}
	
	time_t t0 = clock();
	for(int i=0;i<timesteps;i++) {
		if(timesteps >= 100 && !(i%(timesteps/100))) {
			printf("%d%%..",(100*i)/timesteps);
			fflush(stdout);
		}
		system->step();

		sampler->sample();

		if(printPositions) system->printPositionsToFile(positions);
	}

	sampler->finish();

	printf("100%%\n\n");
	printf("Summary:\n");
	printf("Collisions: %d\n",system->collisions);
	printf("Average temperature: %.3f\n",sampler->temperatureSum/sampler->temperatureSamples);
	
	return 0;
}