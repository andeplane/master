#pragma once

#include <random.h>
#include <armadillo>
#include <structs.h>

class System;
class System;
class Molecule;
class Cell;

using namespace arma;
using namespace std;

class DummyCell {
public:
    int node_id;
    int index;
    int index_vector[3];

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

    int num_molecules;
    int new_molecules;

    double *r;
    double *v;
    double *r0;
    bool *atom_moved;

    Cell(System *system);
	void reset();
    int prepare();
    void resize(int n);
    int collide(Random *rnd);
    void collide_molecules(const int &ip0, const int &ip1, const double &v_rel, Random *rnd);
    void update_volume();
    void add_molecule(Molecule &m);
    void add_molecule(double *r_, double *v_, double *r0_);
    void add_molecule(double *r_, double *v_);
    void update_molecule_cells(int dimension);
    void update_molecule_arrays();

    static bool cmp(Cell *c1, Cell *c2);
};
