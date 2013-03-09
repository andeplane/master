#pragma once

class System;

class DSMCTimer
{
public:
    DSMCTimer();
    double moving_t0;
    double moving;
    double moving_global;
    void start_moving();
    void end_moving();

    double colliding_t0;
    double colliding;
    double colliding_global;
    void start_colliding();
    void end_colliding();

    double mpi_t0;
    double mpi;
    double mpi_global;
    void start_mpi();
    void end_mpi();

    double system_initialize_t0;
    double system_initialize;
    double system_initialize_global;
    void start_system_initialize();
    void end_system_initialize();

    void gather_all_nodes(System *system);
};
