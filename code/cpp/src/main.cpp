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
#include <mpi.h>

using namespace std;

int main(int args, char* argv[]) {
    int numprocs, myid;

    MPI_Init(&args,&argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    Settings *settings = new Settings("../dsmc.ini");
    System system;

    system.initialize(settings, myid);

    // StatisticsSampler *sampler = new StatisticsSampler(&system);

    for(int i=0;i<settings->timesteps;i++) {
        system.io->save_state_to_movie_file();
        system.step();

        // sampler->sample();
    }

    // system.io->save_state_to_file_binary();
    // system.io->finalize();

    MPI_Finalize();

    return 0;
}
