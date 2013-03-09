#include <dsmctimer.h>
#include <mpi.h>
#include <system.h>
#include <threadcontrol.h>

DSMCTimer::DSMCTimer() {
    moving = 0;
    colliding = 0;
    mpi = 0;
}

void DSMCTimer::start_moving() {
    moving_t0 = MPI_Wtime();
}

void DSMCTimer::end_moving() {
    moving += MPI_Wtime() - moving_t0;
}

void DSMCTimer::start_colliding() {
    colliding_t0 = MPI_Wtime();
}

void DSMCTimer::end_colliding() {
    colliding += MPI_Wtime() - colliding_t0;
}

void DSMCTimer::start_mpi() {
    mpi_t0 = MPI_Wtime();
}

void DSMCTimer::end_mpi() {
    mpi += MPI_Wtime() - mpi_t0;
}

void DSMCTimer::start_system_initialize() {
    system_initialize_t0 = MPI_Wtime();
}

void DSMCTimer::end_system_initialize() {
    system_initialize += MPI_Wtime() - system_initialize_t0;
}

void DSMCTimer::gather_all_nodes(System *system) {
    colliding_global = 0;
    moving_global = 0;
    mpi_global = 0;
    system_initialize_global = 0;

    MPI_Reduce(&moving,&moving_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&colliding,&colliding_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&mpi,&mpi_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&system_initialize,&system_initialize_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    colliding_global/=system->thread_control.num_nodes;
    mpi_global/=system->thread_control.num_nodes;
    moving_global/=system->thread_control.num_nodes;
    system_initialize_global/=system->thread_control.num_nodes;
}
