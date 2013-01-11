#include <iostream>
#include "math.h"
#include "time.h"
#include <armadillo>
#include <fstream>
#include <system.h>
#include <statisticssampler.h>
#include <defines.h>
#include <CIniFile.h>
#include <unitconverter.h>

using namespace std;

int main(int args, char* argv[]) {
    CIniFile ini;

    ini.load("../dsmc.ini");
    bool print_positions = ini.getbool("print_positions");
    int timesteps = ini.getint("timesteps");
    int print_every_n_step = ini.getint("print_every_n_step");

    System system;
    system.initialize(ini);
    StatisticsSampler sampler(&system, &ini);
    // StatisticsSampler *sampler = new StatisticsSampler(&system, timesteps);

    FILE *positions = 0;
    if(print_positions) {
        positions = fopen("pos.xyz","w");
        system.printPositionsToFile(positions);
    }

    UnitConverter uc;
    FILE *pressure_file = fopen("pressure.dat","w");

    for(int i=0;i<timesteps;i++) {
        if(timesteps >= 100 && !(i%(timesteps/100))) {
            printf("%d%%..",(100*i)/timesteps);
            fflush(stdout);
        }

        system.step();
        sampler.sample();

        if(print_positions && !(i%print_every_n_step)) system.printPositionsToFile(positions);

        mat p = sampler.calculate_global_pressure_tensor();
        vector<mat> pressure_tensors = sampler.calculate_local_pressure_tensor();

        for(int n=0;n<pressure_tensors.size();n++) {
            mat pressure_tensor = pressure_tensors[n];
            fprintf(pressure_file,"%.5f ",pressure_tensor(0,0));
        }
        fprintf(pressure_file,"\n");

    }
    sampler.calculate_pressure();
    sampler.calculate_velocity_field();

    sampler.finish();

    printf("100%%\n\n");
    printf("Time consumption: \n");
    printf("Sorting: %.2f s\n",system.time_consumption[SORT]);
    printf("Moving: %.2f s\n",system.time_consumption[MOVE]);
    printf("Collisions: %.2f s\n",system.time_consumption[COLLIDE]);
    printf("Sampling: %.2f s\n",system.time_consumption[SAMPLE]);
    printf("Summary:\n");
    printf("Collisions: %d\n",system.collisions);

    printf("System volume: %f\n",system.volume);
    printf("Average temperature: %.3f\n",sampler.get_temperature());
    printf("Average pressure: %.3f\n",sampler.get_pressure());
    printf("V*P: %.3f\n",sampler.get_pressure()*system.volume);

    return 0;
}
