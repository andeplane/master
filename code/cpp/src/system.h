#pragma once

class Molecule;
class Cell;
class Sorter;
class Wall;
class Grid;
class DSMC_IO;

#include <fstream>
#include "molecule.h"
#include "cell.h"
#include "sorter.h"
#include "random.h"
#include <vector>
#include <CIniFile.h>
#include <grid.h>
#include <settings.h>
#include <unitconverter.h>

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
    vector< vector< vector<Cell*> > > cells;

    DSMC_IO *io;
    Grid *world_grid;
    Grid *initial_world_grid;
    Settings *settings;
    UnitConverter * unit_converter;

    Random *rnd;

	int N; 			// Number of molecules

    double Lx;
    double Ly;
    double Lz;
    double acceleration;
    double max_x_acceleration;
	double volume;
	double eff_num;
	double mpv; 	// Most probable velocity
	double mfp; 	// Mean free path
	double dt;
	double t;
    double temperature;
    double mass, diam, density;
    double wall_temperature;
	double *time_consumption;

    double *positions;
    double *velocities;
    double *initial_positions;

	int collisions;
	int steps;

	Sorter *sorter;

    void initialize(Settings *settings_);
	void step();
    System() { }
};
