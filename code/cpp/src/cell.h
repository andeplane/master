#pragma once

class System;
#include <random.h>
#include <armadillo>
class System;
class Molecule;
class Cell;

using namespace arma;
using namespace std;

class DummyCell {
public:
    int node_id;
    int index;

    vector<Molecule*> new_molecules;
    Cell *real_cell;

    DummyCell() { real_cell = NULL; }
};

class Cell {
public:
	double volume;
    int pixels; // Used to calculate volume in a cell
    int total_pixels;
    int index;
    int collision_pairs;

    double x0, y0, z0;
    double Lx, Ly, Lz;
	double vr_max;
    double collision_coefficient;

    DummyCell *dummy_cell;
    System *system;
    vector<Molecule*> molecules;
    int num_molecules;

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
