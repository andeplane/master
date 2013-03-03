#ifndef MOLECULE_H
#define MOLECULE_H
class System;

#include <iostream>
#include <armadillo>
#include <system.h>
#include <random.h>

using namespace arma;
using namespace std;

class Molecule {
public:
    char   *type;
	int    atoms;
    int    index;
    bool   active;
    int    information_carrier;
    double mass;

    double *r;
    double *v;
    double *initial_r;
    vec temp_vector;
    vec vrel;
    vec vcm;

    System *system;

    Molecule(System *system);
    inline void fixR();
    inline void do_move(const double &dt);
    void move(double dt, Random *rnd, int depth = 0);
    void collide_with(Molecule *m, Random *rnd, const double &cr);
    double squared_distance_from_initial_position();
};

#endif
