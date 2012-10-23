#ifndef MOLECULE_H
#define MOLECULE_H
class System;

#include <iostream>
#include <armadillo>
#include "System.h"
#include <Random.h>

using namespace arma;
using namespace std;

class Molecule {
public:
	char   type;
	int    atoms;
    int    index;
    bool   active;

    vec r;
    vec v;

    System *system;

    Molecule(System *system);
    inline void fixR();
	void addR(vec dr);
    void move(double dt, Random *rnd, int depth = 0);
    void move_old(double dt, Random *rnd);
};

#endif
