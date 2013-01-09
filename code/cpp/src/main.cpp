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

    for(int i=0;i<timesteps;i++) {
        if(timesteps >= 100 && !(i%(timesteps/100))) {
            printf("%d%%..",(100*i)/timesteps);
            fflush(stdout);

            double diffusion_constant = sampler.calculate_diffusion_constant();
            cout << "Diffusion constant: " << uc.diffusion_to_SI(diffusion_constant) << endl;
        }
        system.step();

        sampler.sample();

        if(print_positions && !(i%print_every_n_step)) system.printPositionsToFile(positions);
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
