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
    bool print_positions = ini.getbool("print_positions");
    int timesteps = ini.getint("timesteps");

    System system;
    system.initialize(ini);
    StatisticsSampler sampler(&system, ini);
    // StatisticsSampler *sampler = new StatisticsSampler(&system, timesteps);

    FILE *positions = 0;
    if(print_positions) {
        positions = fopen("pos.xyz","w");
        system.printPositionsToFile(positions);
    }

    for(int i=0;i<timesteps;i++) {
        if(timesteps >= 100 && !(i%(timesteps/100))) {
            printf("%d%%..",(100*i)/timesteps);
            fflush(stdout);
        }
        system.step();

        sampler.sample();

        if(print_positions) system.printPositionsToFile(positions);
    }

    sampler.finish();

    printf("100%%\n\n");
    printf("Time consumption: \n");
    printf("Sorting: %.2f s\n",system.time_consumption[SORT]);
    printf("Moving: %.2f s\n",system.time_consumption[MOVE]);
    printf("Collisions: %.2f s\n",system.time_consumption[COLLIDE]);
    printf("Sampling: %.2f s\n",system.time_consumption[SAMPLE]);
    printf("Summary:\n");
    printf("Collisions: %d\n",system.collisions);
    if(ini.getbool("print_temperature"))
        printf("Average temperature: %.3f\n",sampler.temperature_sum/sampler.temperature_samples);


    return 0;
}
