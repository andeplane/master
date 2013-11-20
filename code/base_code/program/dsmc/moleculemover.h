#pragma once
#include <grid.h>
#include <colliderbase.h>

class Cell;
class System;
class Random;
class CVector;
class MoleculeMover
{
public:
    System *system;
    unsigned char *voxels;
    Grid *grid;
    long count_periodic[3];
    ColliderBase *surface_collider;

    MoleculeMover();
    void initialize(System *system_, ColliderBase *surface_collider_);
    void move_molecules(double dt, Random *rnd);
    void move_molecule(const int &idx);
    void do_move(data_type *r, data_type *v, double dt);
    void move_molecule(int &idx, double dt, Random *rnd, int depth);
    void move_molecule_cylinder(int &idx, double dt, Random *rnd, int depth);
    void move_molecule_box(int &molecule_index, double dt, Random *rnd, int depth);
    void apply_periodic_boundary_conditions(int &molecule_index, vector<data_type> &r, const CVector &system_length);
};


