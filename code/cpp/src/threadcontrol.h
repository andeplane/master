#pragma once
#include <vector>

#define MAX_MOLECULE_NUM 100000

class System;
class Settings;
class DummyCell;
class Cell;

using namespace std;

class ThreadControl
{
public:
    int myid;
    int my_vector_index[3];
    int num_processors[3];
    int neighbor_nodes[6];

    int num_nodes;
    double porosity;
    int num_molecules;
    double *origo;
    System *system;
    Settings *settings;

    vector<Cell*> all_cells;
    vector<Cell*> my_cells;
    double *mpi_receive_buffer;
    double **molecules_to_be_moved;
    int *num_molecules_to_be_moved;

    double *new_molecules;
    bool *new_molecule_moved;
    int num_new_molecules;

    double *r;
    double *v;
    double *r0;

    unsigned long *molecule_index_in_cell;
    unsigned long *molecule_cell_index;
    bool *molecule_moved;

    ThreadControl();
    void setup(System *system);
    void setup_molecules();
    void setup_topology();
    void setup_cells();
    void update_cells();
    void update_cell_volume();
    void update_local_cells();
    void calculate_porosity();
    int cell_index_from_position(double *r);
    inline int cell_index_from_ijk(const int &i, const int &j, const int &k);
    inline void find_position(double *r);
    void update_mpi(int dimension);
    void update_molecule_cells_local();
    void update_molecule_cells(int dimension);
    void update_new_molecules(int dimension);
    void update_molecule_arrays();
};
