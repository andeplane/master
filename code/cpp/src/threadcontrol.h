#pragma once
#include <vector>
#include <armadillo>

class System;
class Molecule;
class Settings;
class DummyCell;

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
    double *mpi_data;

    ThreadControl();
    void setup(System *system);
    void setup_molecules();
    void update_mpi();
    int cell_index_from_molecule(Molecule *m);
    int cell_index_from_ijk(const int &i, const int &j, const int &k);
    inline void find_position(Molecule *m);

};
