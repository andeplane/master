#pragma once
#include <vector>
#include <armadillo>
#include <structs.h>

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
    int my_vector_index[3];
    int num_processors[3];
    int neighbor_nodes[6];

    int num_nodes;
    double porosity;
    int num_particles;
    vec3 origo;
    System *system;
    Settings *settings;

    vector<Cell*> cells;
    vector<DummyCell*> dummy_cells;

    double *mpi_receive_buffer;
    double *positions;
    double *velocities;
    double *initial_positions;

    double **molecules_to_be_moved;
    int *num_molecules_to_be_moved;

    double *new_molecules;
    bool *molecule_moved;
    int num_new_molecules;


    ThreadControl();
    void setup(System *system);
    void setup_molecules();
    void setup_topology();
    void setup_cells();
    void update_cells();
    void update_cell_volume();
    void update_mpi(int dimension);
    void calculate_porosity();
    void add_molecule_to_cell(Molecule &molecule, int cell_index);
    int cell_index_from_position(double *r);
    inline int cell_index_from_ijk(const int &i, const int &j, const int &k);
    inline void find_position(double *r);
    void update_new_molecules(int dimension);
    void add_new_molecules_to_cells();
};
