#include <iostream>
#include "math.h"
#include "time.h"
#include <armadillo>
#include <fstream>
#include <system.h>
#include <statisticssampler.h>
#include <defines.h>
#include <unitconverter.h>
#include <settings.h>
#include <dsmc_io.h>


using namespace std;

int main(int args, char* argv[]) {
    Settings *settings = new Settings("../dsmc.ini");
    System system;

    system.initialize(settings);

    for(int i=0;i<settings->timesteps;i++) {
        if(settings->timesteps >= 100 && !(i%(settings->timesteps/100))) {
            printf("%d%%..",(100*i)/settings->timesteps);
            fflush(stdout);
        }

        system.io->save_state_to_movie_file();
        system.step();
    }
    system.io->save_state_to_file_binary();
    system.io->finalize();

    return 0;

    // sampler.calculate_velocity_field();

    // sampler.finish();

    printf("100%%\n\n");
    printf("Time consumption: \n");
    printf("Sorting: %.2f s\n",system.time_consumption[SORT]);
    printf("Moving: %.2f s\n",system.time_consumption[MOVE]);
    printf("Collisions: %.2f s\n",system.time_consumption[COLLIDE]);
    printf("Sampling: %.2f s\n",system.time_consumption[SAMPLE]);
    printf("Summary:\n");
    printf("Collisions: %d\n",system.collisions);

    printf("System volume: %f\n",system.volume);
    // printf("Average temperature: %.3f\n",sampler.get_temperature());

    return 0;
}
