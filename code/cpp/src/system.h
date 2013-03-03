#pragma once

class Molecule;
class Cell;
class Sorter;
class Wall;
class Grid;

#include <fstream>
#include "molecule.h"
#include "cell.h"
#include "sorter.h"
#include "random.h"
#include <vector>
#include <CIniFile.h>
#include <grid.h>
#include <settings.h>

using namespace std;
using namespace arma;

class System {
private:
    void init_positions();
    void init_velocities();
    void init_molecules();
    void init_cells();
	void move();
    void init_randoms();
	int  collide();
	void accelerate();
public:
    vector<Molecule*>molecules;
    Cell ***cells;
    Cell ***load_balanced_cell_list;
    int  *cells_in_list;

    Grid *world_grid;
    Grid *initial_world_grid;
    Settings *settings;

    Random *rnd;

	int N; 			// Number of molecules

    double width;
    double height;
    double acceleration;
    double max_x_acceleration;
	double volume;
	double eff_num;
	double mpv; 	// Most probable velocity
	double mfp; 	// Mean free path
	double dt;
	double coeff;
	double t;
	double T;
    double mass, diam, density;
    double wall_temperature;
	double *time_consumption;
	int collisions;
	int steps;

	Sorter *sorter;

	void printPositionsToFile(FILE *file);
    void initialize(Settings *settings_);
	void step();
	double rand_gauss(long *idum);
    System() { }
};
