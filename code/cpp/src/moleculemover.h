#pragma once
class Cell;
class System;
class Random;
class Grid;

class MoleculeMover
{
public:
    Cell *current_cell;
    System *system;

    unsigned char *voxels;
    Grid *grid;

    MoleculeMover();
    void initialize(System *system_);
    void move_molecules(Cell *cell, double dt, Random *rnd);
    void move_molecule(const int &idx);
    void do_move(double *r, double *v, double *r0, const double &dt);
    void move_molecule(int &idx, double dt, Random *rnd, int depth);
};
