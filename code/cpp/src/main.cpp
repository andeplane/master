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
#include <dsmctimer.h>

using namespace std;

int main(int args, char* argv[]) {
    int numprocs, myid;

    MPI_Init(&args,&argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    double t_start = MPI_Wtime();

    Settings *settings = new Settings("../dsmc.ini");
    System system;

    int num_nodes = settings->nodes_x*settings->nodes_y*settings->nodes_z;
    if(numprocs != num_nodes) {
        if(myid==0) cout << "Wrong number of processors. " << endl << "Config files says " << num_nodes << ". MPI started with " << numprocs << "." << endl;
        MPI_Finalize();
        return(0);
    }

    system.initialize(settings, myid);
    StatisticsSampler *sampler = new StatisticsSampler(&system);

    for(int i=0;i<settings->timesteps;i++) {
        system.io->save_state_to_movie_file();
        system.step();

        sampler->sample();
    }

    system.io->save_state_to_file_binary();
    system.io->finalize();
    system.timer->gather_all_nodes(&system);
    if(myid==0) {
            double system_initialize_percentage = system.timer->fraction_system_initialize();
            double fraction_moving = system.timer->fraction_moving();
            double fraction_colliding = system.timer->fraction_colliding();
            double fraction_io = system.timer->fraction_io();
            double fraction_mpi = system.timer->fraction_mpi();
            double fraction_total = fraction_moving + fraction_colliding + fraction_io + fraction_mpi;
            double time_total = system.timer->system_initialize + system.timer->moving + system.timer->colliding + system.timer->io + system.timer->mpi;

            double total_time = MPI_Wtime() - t_start;
            cout.precision(2);
            cout << endl << "Program finished after " << total_time << " seconds. Time analysis:" << endl;
            cout << fixed
                 << "      System initialize : " << system.timer->system_initialize << " s ( " << 100*system_initialize_percentage << "% )" <<  endl
                 << "      Moving            : " << system.timer->moving << " s ( " << 100*fraction_moving << "% )" <<  endl
                 << "      Colliding         : " << system.timer->colliding << " s ( " << 100*fraction_colliding << "% )" <<  endl
                 << "      Disk IO           : " << system.timer->io << " s ( " << 100*fraction_io << "% )" <<  endl
                 << "      MPI communication : " << system.timer->mpi << " s ( " << 100*fraction_mpi << "% )" <<  endl << endl
                 << "      TOTAL             : " << time_total << " s ( " << 100*fraction_total << "% )" <<  endl;
            cout << endl << settings->timesteps / total_time << " timesteps / second. " << endl;
            cout << system.num_molecules_global*settings->timesteps / (1000*total_time) << "k atom-timesteps / second. " << endl;
            cout << system.num_molecules_global*settings->timesteps / (1000*total_time*numprocs) << "k atom-timesteps / second (per node). " << endl;
        }


    MPI_Finalize();

    return 0;
}
