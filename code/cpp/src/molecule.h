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

    vec r;
    vec initial_r;
    vec v;

    System *system;

    Molecule(System *system);
    inline void fixR();
    inline void addR(vec dr);
    void move(double dt, Random *rnd, int depth = 0);
    void move_old(double dt, Random *rnd);
    double r_squared_from_initial();
};

#endif
