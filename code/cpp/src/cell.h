#pragma once
#include <vector>

class Random;
class System;
class Cell;

using namespace std;

class Cell {
public:
    int node_id;
    int node_index_vector[3]; // <node_x, node_y, node_z>
    int node_delta_index_vector[3]; // Same as above, but with this node as origo

	double volume;
    int pixels; // Used to calculate volume in a cell
    int total_pixels;
    int index;
    int collision_pairs;

    double x0, y0, z0;
    double Lx, Ly, Lz;
	double vr_max;
    double collision_coefficient;

    bool is_dummy_cell;
    System *system;
    vector<int> molecules;
    int num_molecules;

    Cell(System *system);
	void reset();
    unsigned long prepare();
    void resize(int n);
    int collide(Random *rnd);
    void collide_molecules(double *v0, double *v1, const double &v_rel, Random *rnd);
    void update_volume();
    void add_molecule(const int &molecule_index, unsigned long *index_in_cell, unsigned long *cell_index);
    void remove_molecule(const int &molecule_index, unsigned long *index_in_cell);

    static bool cmp(Cell *c1, Cell *c2);
};
