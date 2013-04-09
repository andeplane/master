#pragma once
class Cell;
class System;
class Random;
class Grid;
class ThreadControl;

class MoleculeMover
{
public:
    System *system;
    ThreadControl *thread_control;
    unsigned char *voxels;
    Grid *grid;
    long count_periodic_x;
    long count_periodic_y;
    long count_periodic_z;

    MoleculeMover();
    void initialize(System *system_);
    void move_molecules(double dt, Random *rnd);
    void move_molecule(const int &idx);
    inline void do_move(double *r, double *v, double *r0, const double &dt);
    void move_molecule(int &idx, double dt, Random *rnd, int depth);
};
