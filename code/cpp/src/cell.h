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
    int node_index_vector[3]; // <node_x, node_y, node_z>
    int node_delta_index_vector[3]; // Same as above, but with this node as origo

    int index;
    int test_value;

    vector<Molecule*> new_molecules;
    Cell *real_cell;

    DummyCell() { real_cell = NULL; test_value = 1337; }
};

class Cell {
public:
	double volume;
    int pixels; // Used to calculate volume in a cell
    int total_pixels;
    int index;
    int collision_pairs;
    int test_value;

    double x0, y0, z0;
    double Lx, Ly, Lz;
	double vr_max;
    double collision_coefficient;

    DummyCell *dummy_cell;
    System *system;
    vector<int> molecules;
    int num_molecules;

    Cell(System *system);
	void reset();
    int prepare();
    void resize(int n);
    int collide(Random *rnd);
    void collide_molecules(double *v0, double *v1, const double &v_rel, Random *rnd);
    void update_volume();
    void add_molecule(const int &molecule_index, unsigned long *index_in_cell, unsigned long *cell_index);
    void remove_molecule(const int &molecule_index, unsigned long *index_in_cell);

    static bool cmp(Cell *c1, Cell *c2);
};
