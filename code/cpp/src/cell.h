#pragma once

class System;
#include <random.h>
#include <armadillo>
class System;
class Molecule;

using namespace arma;
using namespace std;

class Cell {
public:
	double volume;
    int pixels; // Used to calculate volume in a cell
    int total_pixels;

	double vr_max;
    System *system;
    double collision_coefficient;

    int collision_pairs;
    int i,j,k;
    int particles;
    int particle_capacity;
    unsigned int *particle_indices;
    Molecule *first_molecule;

    Cell(System *system);
	void reset();
    int prepare();
    void resize(int n);
    int collide(Random *rnd);

    void update_volume();
    void add_molecule(Molecule *m);
    void remove_molecule(Molecule *m);

    static bool cmp(Cell *c1, Cell *c2);
};
