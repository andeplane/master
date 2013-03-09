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
    int num_particles;
    vec3 origo;
    System *system;
    Settings *settings;

    vector<Cell*> cells;
    vector<DummyCell*> dummy_cells;
    vector< vector<Molecule*> > nodes_new_atoms_list;
    vector<Molecule*> free_molecules;
    vector<Molecule*> all_molecules;

    double *mpi_data;
    double *positions;
    double *velocities;
    double *initial_positions;
    int allocated_particle_data;

    ThreadControl();
    void setup(System *system);
    void setup_molecules();
    void setup_cells();
    void update_cells();
    void update_cell_volume();
    void update_mpi();
    void update_local_cells();
    void calculate_porosity();
    int cell_index_from_molecule(Molecule *m);
    inline int cell_index_from_ijk(const int &i, const int &j, const int &k);
    inline void find_position(Molecule *m);
};
