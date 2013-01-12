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

    vec r;
    vec initial_r;
    vec v;

    System *system;

    Molecule(System *system);
    inline void fixR();
    inline void addR(vec dr);
    void move(double dt, Random *rnd, int depth = 0);
    void move_old(double dt, Random *rnd);
    vec collide_with(Molecule *m, Random *rnd, double cr);
    double squared_distance_from_initial_position();
};

#endif
