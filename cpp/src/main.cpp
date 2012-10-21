#include <iostream>
#include "math.h"
#include "time.h"
#include <armadillo>
#include <fstream>
#include <System.h>
#include <StatisticsSampler.h>
#include <defines.h>
#include <CIniFile.h>

using namespace std;

int main(int args, char* argv[]) {

    CIniFile ini;

    ini.load("dsmc.ini");
    bool printPositions = ini.getbool("printPositions");


    int N = ini.getint("N");
    double T = ini.getdouble("T");
    int timesteps = ini.getint("timesteps");
    int threads = ini.getint("threads");



    System system;
    system.initialize(N,T,threads);

    StatisticsSampler *sampler = new StatisticsSampler(&system, timesteps);

    FILE *positions = 0;
    if(printPositions) {
        positions = fopen("pos.xyz","w");
        system.printPositionsToFile(positions);
    }

    for(int i=0;i<timesteps;i++) {
        if(timesteps >= 100 && !(i%(timesteps/100))) {
            printf("%d%%..",(100*i)/timesteps);
            fflush(stdout);
        }
        system.step();

        sampler->sample();

        if(printPositions) system.printPositionsToFile(positions);
    }

    sampler->finish();

    printf("100%%\n\n");
    printf("Time consumption: \n");
    printf("Sorting: %.2f s\n",system.time_consumption[SORT]);
    printf("Moving: %.2f s\n",system.time_consumption[MOVE]);
    printf("Collisions: %.2f s\n",system.time_consumption[COLLIDE]);
    printf("Sampling: %.2f s\n",system.time_consumption[SAMPLE]);
    printf("Summary:\n");
    printf("Collisions: %d\n",system.collisions);
    printf("Average temperature: %.3f\n",sampler->temperatureSum/sampler->temperatureSamples);


    return 0;
}
