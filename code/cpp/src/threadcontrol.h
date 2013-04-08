#pragma once
#include <vector>
#include <armadillo>

class System;
class Molecule;
class Settings;
class DummyCell;
class Cell;

using namespace arma;
using namespace std;

class ThreadControl
{
public:
    int myid;
    int idx, idy, idz;
    int num_nodes;
    double porosity;
    int num_molecules;
    vec3 origo;
    System *system;
    Settings *settings;

    vector<Cell*> cells;
    vector<DummyCell*> dummy_cells;
    double *mpi_receive_buffer;
    double **molecules_to_be_moved;
    int *num_molecules_to_be_moved;

    double *new_molecules;
    bool *molecule_moved;
    int num_new_molecules;

    double *r;
    double *v;
    double *r0;
    unsigned long *molecule_index_in_cell;
    unsigned long *molecule_cell_index;

    ThreadControl();
    void setup(System *system);
    void setup_molecules();
    void setup_cells();
    void update_cells();
    void update_cell_volume();
    void update_mpi();
    void update_local_cells();
    void calculate_porosity();
    int cell_index_from_position(double *r);
    inline int cell_index_from_ijk(const int &i, const int &j, const int &k);
    inline void find_position(double *r);
};
