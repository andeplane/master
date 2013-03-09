#pragma once
class System;
class Random;
class Cell;

#include <iostream>

using namespace std;

class Molecule {
public:
    char   *type;
    int    atoms;
    bool   active;
    double mass;
    int       cell_index;
    int       index_in_cell;

    double *r;
    double *v;
    double *r_initial;

    System *system;

    Molecule(System *system);
    inline void fixR();
    inline void do_move(const double &dt);
    void move(double dt, Random *rnd, int depth = 0);
    void collide_with(Molecule *m, Random *rnd, const double &cr);
    double squared_distance_from_initial_position();
};
