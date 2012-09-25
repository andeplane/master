#include <iostream>
#include "math.h"
#include "time.h"
#include <armadillo>
#include <fstream>
#include "System.h"
#include "StatisticsSampler.h"

using namespace std;

int main(int args, char* argv[]) {
	int N = args > 1 ? atoi(argv[1]) : 1000;
	int T = args > 2 ? atof(argv[2]) : 1;
	int timesteps = args > 4 ? atof(argv[4]) : 1000;

	System *system = new System(N);
	StatisticsSampler *sampler = new StatisticsSampler(system);

	// ofstream *file = new ofstream;
	// file->open("pos.xyz");
	FILE *positions = fopen("pos.xyz","w");
	double t = 0;

	system->printPositionsToFile(positions);

	for(int i=0;i<timesteps;i++) {
		if(!(i%(timesteps/100))) {
			printf("%d%%..",(100*i)/timesteps);
			fflush(stdout);
		}

		system->step();
		// sampler->sample(t);

		system->printPositionsToFile(positions);
	}

	return 0;
}